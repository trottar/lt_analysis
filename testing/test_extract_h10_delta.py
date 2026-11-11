import contextlib
import hashlib
import importlib
import io
import json
import os
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

REPO = Path(__file__).resolve().parents[1]
for folder in ("testing", "src/cuts", "src/binning", "src/utility"):
    sys.path.insert(0, str(REPO/folder))
import extract_h10_delta as target
from canonical_binning import find_canonical_bin
from pion_t_template_comparison import integrate_raw_templates

SETTINGS = ("Center_highe", "Center_lowe", "Left_lowe", "Left_highe", "Right_highe")
CHANNELS = {"pi_n":"neutron", "pi_delta":"delta", "pi_sidis":"sidis"}
T_EDGES = [.4, .5666666666666667, .7333333333333334, .9]
PHI_EDGES = list(np.linspace(-180,180,10))
MM_EDGES = [.7,.9,1.1,1.3]
SOURCE_FILES = ("src/cuts/apply_cuts.py", "src/cuts/data_coordinates.py",
                "src/binning/theta_cm.py", "src/utility/canonical_binning.py",
                "src/cuts/pion_t_template_comparison.py")


def tree_events(extra=False):
    records = []
    for t in (.45,.65,.85):
        for mm, weight, delta in ((.6,1.,-10.),(.95,2.,-7.),(1.15,1.5,20.),(1.4,.5,22.)):
            fields = dict.fromkeys(target.REQUIRED,0.)
            fields.update(missmass=mm,t=-t,phipq=0.,ssdelta=delta,Q2=4.4,W=2.74,
                          iter_weight=weight+(0.2 if extra else 0.))
            records.append(fields)
    return records


def row_arrays(events):
    vals = {k:[] for k in ("entry","weight","minus_t","phi_degrees","delta_bin",
                           "t_bin","phi_bin","mm_bin",*target.REQUIRED)}
    for i,e in enumerate(events):
        v = dict(e,entry=i,weight=e.get("iter_weight",1.),minus_t=-e["t"],
                 phi_degrees=e["phipq"]*180./np.pi)
        for name,coordinate,edges in (("delta_bin",v["ssdelta"],target.DELTA),
                                      ("t_bin",v["minus_t"],T_EDGES),
                                      ("phi_bin",v["phi_degrees"],PHI_EDGES)):
            v[name] = find_canonical_bin(coordinate,edges)
        v["mm_bin"] = np.searchsorted(MM_EDGES,v["missmass"],side="right")
        for k in vals:
            vals[k].append(v[k])
    return {k:np.asarray(v,dtype=int if k=="entry" or k.endswith("_bin") else float)
            for k,v in vals.items()}


def fixture_role(path, identity="Center_highe_pi_n"):
    source = dict(status="available",path=str(path),sha256=target.digest(path),
                  bytes=path.stat().st_size,tree_name="h10",tree_entries=12,
                  weight_branch_present=True,
                  weight_expression="getattr(evt, 'iter_weight', 1.0)",
                  coordinate="missmass; -t; phipq*180/pi",
                  hgcer_hole_points=[[30.,30.],[31.,30.],[31.,31.],[30.,30.]],
                  active_selection_state=dict(W="2p74",Q2="4p4",EPSSET="high",POL="+1",
                                              ParticleType="kaon",tmin=.4,tmax=.9,
                                              cut_poly=[[4.,2.],[5.,2.],[5.,3.],[4.,3.]],
                                              TMIN_RESOLUTION_THRESHOLD=1e-5))
    records = row_arrays(tree_events())
    raw = integrate_raw_templates(zip(records["missmass"],records["minus_t"],
                                      records["phi_degrees"],records["weight"]),
                                  T_EDGES,PHI_EDGES,MM_EDGES,source)
    role = dict(id=identity,source=source,t_edges=T_EDGES,phi_edges=PHI_EDGES,parents=[])
    for key in ("selected","accepted","excluded"):
        role[key] = raw[key]
    for p in raw["parents"]:
        d = dict(mm_edges=MM_EDGES,raw=p["raw"],children=p["raw_children"],
                 shift=.03,interpolation="linear",renormalize=True)
        L = target.transform(d)
        d["aligned"] = (L@d["raw"]["values"]).tolist()
        d["aligned_sumw2"] = (L**2@d["raw"]["variances"]).tolist()
        role["parents"].append(d)
    return role


class FakeTree:
    def __init__(self, events, missing=()):
        self.events = events
        self.missing = set(missing)

    def GetEntries(self):
        return len(self.events)

    def GetBranch(self, key):
        return key not in self.missing and (key in target.REQUIRED or key=="iter_weight")

    GetLeaf = GetBranch

    def GetListOfBranches(self):
        names=[k for k in (*target.REQUIRED,"iter_weight") if self.GetBranch(k)]
        return [SimpleNamespace(GetName=lambda n=n:n) for n in names]

    def __iter__(self):
        return iter(SimpleNamespace(**r) for r in self.events)


class FakeFile:
    def __init__(self, events, missing=()):
        self.tree = FakeTree(events,missing)
        self.closed = False

    def IsZombie(self):
        return False

    def Get(self,name):
        return self.tree if name=="h10" else None

    def Close(self):
        self.closed = True


class ProjectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.old = self.root/"old.root"
        self.old.write_bytes(b"historical-root")
        self.role = fixture_role(self.old)
        self.rows = row_arrays(tree_events())

    def test_historical_parent_closure_with_explicit_flows(self):
        receipt, arrays = target.check_role(self.rows,self.role)
        self.assertEqual(len(receipt),3)
        self.assertEqual(receipt[0]["delta_flow_count"],1)
        self.assertEqual(receipt[0]["mm_underflow"],1)
        self.assertEqual(receipt[0]["mm_overflow"],1)
        self.assertEqual(arrays["t1_count"].shape,(11,5))

    def test_changed_sample_fails_historical_histograms(self):
        with self.assertRaises(AssertionError):
            target.check_role(row_arrays(tree_events(extra=True)),self.role)

    def test_regenerated_reference_closes_changed_weights(self):
        admitted=dict(identity=target.input_identity(self.old),tree_entries=12,weight_branch_present=True)
        new, comparison = target.regenerated_role(row_arrays(tree_events(extra=True)),self.role,admitted)
        target.check_role(row_arrays(tree_events(extra=True)),new)
        self.assertFalse(comparison[0]["historical_equality_required"])
        self.assertTrue(any(abs(x)>0 for x in comparison[0]["aligned_values_difference"]))
        self.assertNotEqual(new["parents"][0]["raw"]["values"],self.role["parents"][0]["raw"]["values"])

    def test_regenerated_signed_moments(self):
        events = tree_events()
        events[3]["iter_weight"] = -.1
        rows = row_arrays(events)
        admitted=dict(identity=target.input_identity(self.old),tree_entries=12,weight_branch_present=True)
        new,_ = target.regenerated_role(rows,self.role,admitted)
        target.check_role(rows,new)
        self.assertGreater(new["selected"]["absolute"],new["selected"]["signed"])

    def test_zero_weight_parent_is_unavailable(self):
        rows = {k:v.copy() for k,v in self.rows.items()}
        rows["weight"][:] = 0
        admitted=dict(identity=target.input_identity(self.old),tree_entries=12,weight_branch_present=True)
        with self.assertRaisesRegex(RuntimeError,"parent unavailable"):
            target.regenerated_role(rows,self.role,admitted)

    def test_historical_file_hash_guard(self):
        self.assertEqual(target.find_input(self.role,[]),self.old.resolve())
        self.old.write_bytes(b"regenerated")
        with self.assertRaisesRegex(RuntimeError,"exact input unresolved"):
            target.find_input(self.role,[])

    def test_regenerated_explicit_root_never_uses_old_path(self):
        fresh = self.root/"fresh"/self.old.parent.name/self.old.name
        fresh.parent.mkdir(parents=True)
        fresh.write_bytes(b"new")
        self.assertEqual(target.find_input(self.role,[],self.root/"fresh"),fresh)
        with self.assertRaisesRegex(RuntimeError,"regenerated input missing"):
            target.find_input(self.role,[],self.root/"absent")

    def test_identity_detects_mutation_during_hash(self):
        digest = target.digest
        def mutate(path):
            value = digest(path)
            path.write_bytes(b"mutated")
            return value
        with patch.object(target,"digest",side_effect=mutate):
            with self.assertRaisesRegex(RuntimeError,"during hash"):
                target.input_identity(self.old)

    def test_regenerated_entries_and_weight_branch_are_actual(self):
        root = SimpleNamespace(TFile=SimpleNamespace(Open=lambda *a:FakeFile(tree_events()[:5],("iter_weight",))))
        f,t,optional=target.open_tree(root,self.role,self.old,True)
        self.assertEqual(t.GetEntries(),5)
        self.assertFalse(t.GetBranch("iter_weight"))
        f.Close()
        with self.assertRaisesRegex(RuntimeError,"entries changed"):
            target.open_tree(root,self.role,self.old,False)

    def test_missing_branch_closes_input(self):
        f=FakeFile(tree_events(),("ssdelta",))
        root=SimpleNamespace(TFile=SimpleNamespace(Open=lambda *a:f))
        with self.assertRaisesRegex(RuntimeError,"missing branches"):
            target.open_tree(root,self.role,self.old,True)
        self.assertTrue(f.closed)

    def test_empty_tree_rejected(self):
        root=SimpleNamespace(TFile=SimpleNamespace(Open=lambda *a:FakeFile([])))
        with self.assertRaisesRegex(RuntimeError,"empty tree"):
            target.open_tree(root,self.role,self.old,True)

    def test_exact_delta_endpoints(self):
        self.assertEqual([find_canonical_bin(x,target.DELTA) for x in target.DELTA],list(range(10))+[9])
        self.assertEqual(find_canonical_bin(20.001,target.DELTA),-1)
        self.assertEqual(find_canonical_bin(.9,T_EDGES),2)

    def main_fixture(self):
        root=self.root/"new"
        roles=[]
        for setting in SETTINGS:
            for channel,subdir in CHANNELS.items():
                name=f"Prod_Coin_Q4p4W2p74{setting.split('_')[0].lower()}_{setting.split('_')[1]}.root"
                old=self.root/"original"/subdir/name
                old.parent.mkdir(parents=True,exist_ok=True)
                old.write_bytes(b"historical "+setting.encode()+channel.encode())
                roles.append(fixture_role(old,setting+"_"+channel))
                actual=root/subdir/name
                actual.parent.mkdir(parents=True,exist_ok=True)
                actual.write_text(json.dumps(tree_events(extra=True)),encoding="utf-8")
        ref=dict(delta_edges=target.DELTA,roles=roles,primary_zip_sha256="b"*64,
                 source_hashes={p:hashlib.sha256((REPO/p).read_bytes().replace(b"\r\n",b"\n")).hexdigest() for p in SOURCE_FILES})
        reference=self.root/"reference.json"
        reference.write_text(json.dumps(ref),encoding="utf-8")
        output=self.root/"projection.zip"
        class Hole:
            def __init__(self,*a):pass
            def SetPoint(self,*a):pass
            def IsInside(self,*a):return False
        fake=SimpleNamespace(TFile=SimpleNamespace(Open=lambda p,*a:FakeFile(json.loads(Path(p).read_bytes()))),
                             TCutG=Hole,gROOT=SimpleNamespace(SetBatch=lambda *a:None,GetVersion=lambda:"fake"))
        args=["extract_h10_delta.py","--repo",str(REPO),"--reference",str(reference),
              "--regenerated-input-root",str(root),"--output",str(output)]
        return fake,args,output,root,reference

    def invoke(self,fake,args):
        import apply_cuts
        with patch.dict(sys.modules,ROOT=fake), patch.object(sys,"argv",args),              patch.object(apply_cuts,"apply_simc_sub_cuts",return_value=True),              contextlib.redirect_stdout(io.StringIO()):
            target.main()

    def test_main_regenerated_all_15_and_receipt(self):
        fake,args,output,_,_=self.main_fixture()
        self.invoke(fake,args)
        import zipfile
        with zipfile.ZipFile(output) as z:
            receipt=json.loads(z.read("receipt.json"))
            self.assertEqual(receipt["parent_channel_closures"],45)
            self.assertEqual(receipt["parent_closure_reference"],"current_inputs")
            self.assertFalse(receipt["historical_template_equivalent"])
            self.assertFalse(receipt["scientifically_accepted"])
            self.assertEqual(len(receipt["roles"]),15)
            self.assertIn("ssdelta",receipt["roles"][0]["input"]["available_branches"])
            self.assertNotEqual(receipt["roles"][0]["sha256"],receipt["roles"][0]["historical_source"]["sha256"])
            for name,meta in receipt["members"].items():
                self.assertEqual(hashlib.sha256(z.read(name)).hexdigest(),meta["sha256"])

    def test_main_no_overwrite(self):
        fake,args,output,_,_=self.main_fixture()
        output.write_bytes(b"independent-existing-file")
        with self.assertRaisesRegex(RuntimeError,"already exists"):
            self.invoke(fake,args)
        self.assertEqual(output.read_bytes(),b"independent-existing-file")

    def test_main_duplicate_physical_file(self):
        fake,args,output,root,_=self.main_fixture()
        paths=sorted(root.glob("*/*.root"))
        paths[1].unlink()
        os.link(paths[0],paths[1])
        with self.assertRaisesRegex(RuntimeError,"duplicate physical input"):
            self.invoke(fake,args)
        self.assertFalse(output.exists())

    def test_main_changed_input_before_publication(self):
        fake,args,output,root,_=self.main_fixture()
        original=target.check_role
        victim=next(root.glob("*/*.root"))
        calls=[]
        def mutate(rows,role):
            result=original(rows,role)
            calls.append(role["id"])
            if len(calls)==15:
                victim.write_bytes(b"replaced after earlier read")
            return result
        with patch.object(target,"check_role",side_effect=mutate):
            with self.assertRaisesRegex(RuntimeError,"before publication"):
                self.invoke(fake,args)
        self.assertFalse(output.exists())

    def test_main_source_pin_failure(self):
        fake,args,output,_,reference=self.main_fixture()
        ref=json.loads(reference.read_bytes())
        ref["source_hashes"]["src/cuts/apply_cuts.py"]="0"*64
        reference.write_text(json.dumps(ref),encoding="utf-8")
        with self.assertRaisesRegex(RuntimeError,"source bytes differ"):
            self.invoke(fake,args)
        self.assertFalse(output.exists())

    def test_main_missing_role_reference(self):
        fake,args,output,_,reference=self.main_fixture()
        ref=json.loads(reference.read_bytes())
        ref["roles"][0]["id"]=ref["roles"][1]["id"]
        reference.write_text(json.dumps(ref),encoding="utf-8")
        with self.assertRaisesRegex(RuntimeError,"role identities"):
            self.invoke(fake,args)
        self.assertFalse(output.exists())

    def test_main_conflicting_modes(self):
        fake,args,output,root,_=self.main_fixture()
        with self.assertRaisesRegex(RuntimeError,"one input admission mode"):
            self.invoke(fake,args+["--input-root",str(root)])
        self.assertFalse(output.exists())


    def test_main_historical_default_is_preserved(self):
        fake,args,output,_,reference=self.main_fixture()
        ref=json.loads(reference.read_bytes())
        for role in ref["roles"]:
            path=Path(role["source"]["path"])
            path.write_text(json.dumps(tree_events()),encoding="utf-8")
            role["source"]["sha256"]=target.digest(path)
            role["source"]["bytes"]=path.stat().st_size
        reference.write_text(json.dumps(ref),encoding="utf-8")
        offset=args.index("--regenerated-input-root")
        args=args[:offset]+args[offset+2:]
        self.invoke(fake,args)
        import zipfile
        with zipfile.ZipFile(output) as archive:
            receipt=json.loads(archive.read("receipt.json"))
        self.assertEqual(receipt["cohort_kind"],"historical")
        self.assertEqual(receipt["parent_closure_reference"],"historical_reference")
        self.assertTrue(receipt["historical_template_equivalent"])

    def test_main_regenerated_changed_entry_count(self):
        fake,args,output,root,_=self.main_fixture()
        path=next(root.glob("*/*.root"))
        events=json.loads(path.read_bytes())
        events.append(events[0])
        path.write_text(json.dumps(events),encoding="utf-8")
        self.invoke(fake,args)
        import zipfile
        with zipfile.ZipFile(output) as archive:
            receipt=json.loads(archive.read("receipt.json"))
        changed=[r for r in receipt["roles"] if r["input"]["tree_entries"]==13]
        self.assertEqual(len(changed),1)
        self.assertEqual(changed[0]["historical_source"]["tree_entries"],12)
        self.assertEqual(sum(p["raw_entries"] for p in changed[0]["parents"]),13)

    def test_nonfinite_selected_row_fails_before_publication(self):
        fake,args,output,root,_=self.main_fixture()
        path=next(root.glob("*/*.root"))
        events=json.loads(path.read_bytes())
        events[0]["ssdelta"]=float("nan")
        path.write_text(json.dumps(events),encoding="utf-8")
        with self.assertRaisesRegex(RuntimeError,"nonfinite selected row"):
            self.invoke(fake,args)
        self.assertFalse(output.exists())

    def test_publication_failure_does_not_replace_collision(self):
        output=self.root/"collision.zip"
        output.write_bytes(b"independent")
        with self.assertRaises(FileExistsError):
            target.publish_projection(output,{},{"ok":True})
        self.assertEqual(output.read_bytes(),b"independent")

    def test_publication_context_exit_failure_rolls_back(self):
        output=self.root/"out.zip"
        actual=target.tempfile.TemporaryDirectory
        class Context:
            def __init__(self,*args,**kwargs):
                self.inner=actual(*args,**kwargs)
            def __enter__(self):
                return self.inner.__enter__()
            def __exit__(self,*args):
                self.inner.__exit__(*args)
                raise RuntimeError("context-exit failure")
        with patch.object(target.tempfile,"TemporaryDirectory",Context):
            with self.assertRaisesRegex(RuntimeError,"context-exit failure"):
                target.publish_projection(output,{},{"ok":True})
        self.assertFalse(output.exists())

    def test_publication_context_exit_preserves_replacement(self):
        output=self.root/"out.zip"
        actual=target.tempfile.TemporaryDirectory
        class Context:
            def __init__(self,*args,**kwargs):
                self.inner=actual(*args,**kwargs)
            def __enter__(self):
                return self.inner.__enter__()
            def __exit__(self,*args):
                self.inner.__exit__(*args)
                replacement=output.with_suffix(".replacement")
                replacement.write_bytes(b"independent")
                os.replace(replacement,output)
                raise RuntimeError("first failure")
        with patch.object(target.tempfile,"TemporaryDirectory",Context):
            with self.assertRaisesRegex(RuntimeError,"first failure"):
                target.publish_projection(output,{},{"ok":True})
        self.assertEqual(output.read_bytes(),b"independent")


    def test_reference_changed_before_publication(self):
        fake,args,output,_,reference=self.main_fixture()
        original=target.check_role
        calls=[]
        def mutate(rows,role):
            result=original(rows,role)
            calls.append(role["id"])
            if len(calls)==15:
                reference.write_bytes(reference.read_bytes()+b" ")
            return result
        with patch.object(target,"check_role",side_effect=mutate):
            with self.assertRaisesRegex(RuntimeError,"reference changed"):
                self.invoke(fake,args)
        self.assertFalse(output.exists())


if __name__=="__main__":
    unittest.main()
