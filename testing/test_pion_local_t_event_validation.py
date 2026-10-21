"""Deterministic scalar/synthetic checks; no farm or ROOT integration claims."""
from copy import deepcopy
import ast
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
import zipfile
from unittest.mock import patch

from testing import pion_local_t_event_validation_capture as cap
from testing import pion_local_t_event_validation_replay as replay

ROOT = Path(__file__).resolve().parents[1]
FP = "b"*64
EDGES = [.7,1.05,1.10,1.16,1.3]


def histogram(values, variance=None):
    variance = variance or [0.]*len(values)
    errors = [math.sqrt(x) for x in variance]
    packed = [float(len(EDGES)-1),EDGES[0],EDGES[-1],*EDGES]
    packed += [v for pair in zip(values,errors) for v in pair]
    return {"status":"available", "edges":EDGES[:], "values":values[:], "errors":errors,
            "stored_sumw2":variance[:], "content_error_fingerprint":hashlib.sha256(b"".join(struct.pack("!d",v) for v in packed)).hexdigest()}


def fixture(setting="Center"):
    coeff = cap.coefficients(2.,.5,2.)
    events=[]; child={s:{k:[] for k in ("entry_index","adj_MM","adj_t","phi_degrees","allcuts","nommcuts")} for s in cap.SOURCES}
    kaon={s:{k:[] for k in ("source_label","source_entry_index","adj_MM","adj_t","phi_shift_deg","t_index","phi_index","allcuts","nommcuts","proton_cleaning_factor","proton_rf_accept")} for s in cap.SOURCES}
    parents=[]; crows=[]
    for ti in range(3):
        for s in cap.SOURCES:
            for j,(mm,allcuts) in enumerate(((1.,False),(1.12,True))):
                entry=ti*2+j
                r=dict(source_label=s,entry_index=entry,coefficient=coeff[s],adj_MM=mm,adj_t=ti+.5,
                       phi_degrees=30.,t_index=ti,phi_index=0,allcuts=allcuts,nommcuts=True,coordinate_fingerprint=FP)
                events.append(r)
                for key in child[s]: child[s][key].append(r[key])
                values=(s,entry,mm,ti+.5,30.,ti,0,allcuts,True,.8,True)
                for key,value in zip(kaon[s],values): kaon[s][key].append(value)
        parent=dict(t_bin_index=ti,t_edges=[float(ti),float(ti+1)],pion_parent_id=setting+str(ti),canonical_interval_pair_id="pair",
                    canonical_interval_pair_hash=FP,coordinate_fingerprint=FP,baseline_parent_fingerprint=FP)
        h=histogram([0.,.75,0.,.75,0.,0.],[0.,5.3125,0.,5.3125,0.,0.])
        model=histogram([0.,.75,.1,.75,.1,0.]); k=histogram([0.,.075,.01,.075,.01,0.])
        template=histogram([0.,.25,.25,.25,.25,0.])
        fit={"H_pion_control_input":h,"H_kaon_nosub_input":h,"H_pion_fit_total_refined":model,
             "H_kaon_pion_bg_fit_total_refined":k,"pion_component_alignment":{"components":{"pi_n":{},"pi_sidis":{},"pi_delta":{}}},
             "resolved_subtraction_config":{"synthetic":True}}
        for n in replay.c.COMPONENTS:fit["H_simc_shape_"+n]=template
        p={**parent,"status":"available","authoritative_parent_unchanged":True,"production_applied":False,"rejection_reasons":[],
           "alignment_mode":"fixed_parent_alignment_template_shape_test","baseline_A":deepcopy(fit),"candidate_B":deepcopy(fit),
           "frozen_inputs":{"pion":deepcopy(h),"kaon":deepcopy(h)},"templates":{}}
        for n in replay.c.COMPONENTS:
            p["templates"][n]={"status":"genuinely_local","normalization_count":1,"canonical_selection_closure":True,
                "weighted_child_to_parent_closure":True,"production_global_selection_closure":True,"t_edges":parent["t_edges"],
                "provenance":{"status":"available","sha256":FP,"tree_name":"synthetic","missing_branches":[]},
                "histogram":deepcopy(template),"normalized_values":template["values"],"normalized_variances":template["stored_sumw2"],"normalized_integral":1.}
        transfer={"weights":[0.,.1,.1,.1,.1,0.]}
        crows.append({"setting":setting,"t_bin":ti+1,"t_edges":parent["t_edges"],"A_B_unchanged":p,"fixed_K_B":k["values"],
                      "A_transfer":deepcopy(transfer),"B_transfer":deepcopy(transfer),"C":[{"lambda":lam,"status":"available","selected":False,
                        "model":model["values"],"transfer":deepcopy(transfer)} for lam in replay.c.LAMBDAS]})
        parents.append({**parent,"pion_input":deepcopy(h),"kaon_input":deepcopy(h),
                        "production_weight_configuration":{"weight_denominator_floor":1e-12,"weight_clip_bounds":[None,None]}})
    hist={"normfac_data":2.,"normfac_dummy":.5,"nWindows":2.,"phi_setting":setting,
          "_authoritative_pion_control_source_cache":{"records":events,"coordinate_fingerprint":FP,"kaon_data_coordinate":{"frame":"synthetic"},
              "physical_pion_control_mask_fingerprint":FP,"source_accounting":{s:{"records_inside_canonical_t":6} for s in cap.SOURCES}},
          "_yield_data_event_cache":kaon,"_yield_sub_event_cache":child}
    identities={"data":{"sha256":"a"*64,"bytes":10},"dummy":{"sha256":"d"*64,"bytes":12}}
    packet=cap.build_packet(hist,{},[0.,1.,2.,3.],[0.,60.,120.],source_head="a"*40,input_identity=identities,parent_records=parents,source_files={})
    return packet,crows,hist


def comparison_fixture():
    packets={}; rows=[]; metadata={}
    for s in replay.c.SETTINGS:
        packet,r,_=fixture(s);packets[s]=packet;rows.extend(r)
        metadata[s]={"canonical_t_edges":packet["t_edges"],"canonical_phi_edges":packet["phi_edges"],
                     "original_analysis_input_identity":packet["input_identity"],"source_files":{}}
    candidate={"schema_version":replay.c.SCHEMA,"status":"available","inputs_unchanged":True,"diagnostic_only":True,
               "production_applied":False,"scientifically_accepted":False,"source_HEAD_of_input":replay.c.SOURCE,
               "policy":deepcopy(replay.c.POLICY),"parents":rows,"input_metadata":metadata,
               "provenance":{s:{"sha256":replay.ACCEPTED_AB[s]} for s in replay.c.SETTINGS}}
    return packets,candidate


class BackingRootFile:
    """Only the existing TFile identity interface; never reads a ROOT tree."""
    def __init__(self, path, *, opened=True, zombie=False):
        self.path, self.opened, self.zombie = str(path), opened, zombie
    def GetName(self): return self.path
    def GetFile(self): return self
    def IsOpen(self): return self.opened
    def IsZombie(self): return self.zombie


class BackedTree:
    def __init__(self, name, backing): self.name, self.backing = name, backing
    def GetName(self): return self.name
    def GetCurrentFile(self): return self.backing
    def __iter__(self): raise AssertionError("observer must not iterate a tree")


class CaptureHistogram:
    """Read-only TH1 interface used by the actual retained-histogram serializer."""
    def __init__(self, stored): self.stored = deepcopy(stored)
    def GetNbinsX(self): return len(self.stored["edges"])-1
    def GetXaxis(self): return self
    def GetXmin(self): return self.stored["edges"][0]
    def GetXmax(self): return self.stored["edges"][-1]
    def GetBinLowEdge(self, i): return self.stored["edges"][i-1]
    def GetBinContent(self, i): return self.stored["values"][i]
    def GetBinError(self, i): return self.stored["errors"][i]
    def GetSumw2N(self): return len(self.stored["stored_sumw2"])
    def GetSumw2(self): return self.stored["stored_sumw2"]


class EnabledObserverTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.outdir = self.directory / "packets"
        self.outdir.mkdir()
        self.paths = {role:self.directory/(role+".root") for role in ("data","dummy")}
        for role,path in self.paths.items(): path.write_bytes(("synthetic backing "+role).encode())
        self.identities = {role:cap.file_identity(path) for role,path in self.paths.items()}
        self.inp = {"Q2":"4p4","W":"2p74","ParticleType":"kaon","EPSSET":"low",
                    "particle_subtraction_weight_denom_floor":2.e-8,
                    "particle_subtraction_weight_clip_min":.02,"particle_subtraction_weight_clip_max":7.5}
        # Resolve real source modules, not mocked production policy functions.
        path_patch = patch.object(sys,"path",[str(ROOT/"src/cuts"),str(ROOT/"src/utility"),*sys.path])
        path_patch.start(); self.addCleanup(path_patch.stop)
        env_patch = patch.dict(os.environ,{cap.ENABLE_ENV:"1"})
        env_patch.start(); self.addCleanup(env_patch.stop)
        pins = patch.object(cap,"ACCEPTED_INPUT_IDENTITIES",
                            {s:deepcopy(self.identities) for s in ("Center","Left")})
        pins.start(); self.addCleanup(pins.stop)

    def runtime_hist(self, setting="Center", *, trees=True):
        packet,_,hist = fixture(setting)
        parents=[]
        for p in packet["parents"]:
            parent = {k:deepcopy(v) for k,v in p.items() if k not in (
                "pion_input","kaon_input","production_weight_configuration","baseline_parent_fingerprint")}
            parent.update(H_pion_control=CaptureHistogram(p["pion_input"]),
                          H_proton_cleaned_final_rf=CaptureHistogram(p["kaon_input"]),
                          fit_result={"analysis_scope":"t_bin"+str(p["t_bin_index"]+1)},
                          parent_fit_configuration_hash=FP,
                          final_diagnostic_application_result={"diagnostics":{
                              "denom_floor":2.e-8,"clip_min":.02,"clip_max":7.5}})
            self.assertNotIn("parent_fit_configuration",parent)
            parents.append(parent)
        hist["_pion_t_bin_parent_results"] = tuple(parents)
        for role,key in (("data","InFile_DATA"),("dummy","InFile_DUMMY")):
            backing=BackingRootFile(self.paths[role])
            hist[key]=BackedTree("T_coin_prompt_kaon_noRF",backing) if trees else backing
        return hist,packet

    def observe(self,hist,packet):
        cap.observe(hist,self.inp,packet["t_edges"],packet["phi_edges"],self.outdir)
        path=self.outdir/(hist["phi_setting"]+"_kaon_Q4p4W2p74_lowe_pion-local-t-events.json")
        return replay.strict_json(path.read_bytes())

    def test_enabled_observer_uses_backing_files_and_real_resolved_policy(self):
        from pion_background_chain_audit import digest as source_digest
        for setting,trees in (("Center",True),("Left",False)):
            with self.subTest(setting=setting,trees=trees):
                hist,packet=self.runtime_hist(setting,trees=trees)
                before=source_digest(hist); inp_before=deepcopy(self.inp)
                sources={key:hist[key] for key in ("InFile_DATA","InFile_DUMMY")}
                result=self.observe(hist,packet)
                self.assertEqual(result["input_identity"],self.identities)
                self.assertEqual(result["source_input_paths"],{r:str(p.resolve()) for r,p in self.paths.items()})
                self.assertEqual(cap.validate_packet(result),result)
                self.assertEqual(result["parents"][0]["production_weight_configuration"],
                                 {"weight_denominator_floor":2.e-8,"weight_clip_bounds":[.02,7.5]})
                self.assertEqual(source_digest(hist),before)
                self.assertEqual(self.inp,inp_before)
                for key,source in sources.items():
                    self.assertIs(hist[key],source)
                    backing=source.GetCurrentFile() if trees else source.GetFile()
                    self.assertTrue(backing.IsOpen())
                    self.assertFalse(backing.IsZombie())
                self.assertTrue(all(self.paths[r].read_bytes()==("synthetic backing "+r).encode() for r in self.paths))

    def test_enabled_observer_rejects_missing_backing_provenance(self):
        # A path-looking tree name is still not file provenance.
        for role,key in (("data","InFile_DATA"),("dummy","InFile_DUMMY")):
            for backing in (None,BackingRootFile(self.paths[role],opened=False),
                            BackingRootFile(self.paths[role],zombie=True)):
                hist,packet=self.runtime_hist()
                hist[key]=BackedTree(str(self.paths[role]),backing)
                with self.subTest(role=role,backing=backing),self.assertRaisesRegex(ValueError,"missing_backing_file_provenance:"+role):
                    self.observe(hist,packet)
                self.assertFalse(list(self.outdir.iterdir()))
        class NameOnly:
            def GetName(inner): return str(self.paths["data"])
        hist,packet=self.runtime_hist();hist["InFile_DATA"]=NameOnly()
        with self.assertRaisesRegex(ValueError,"missing_backing_file_provenance:data"):
            self.observe(hist,packet)
        self.assertFalse(list(self.outdir.iterdir()))

    def test_enabled_observer_rejects_input_mismatch_before_parent_serialization(self):
        for role in ("data","dummy"):
            hist,packet=self.runtime_hist()
            self.paths[role].write_bytes(b"different input")
            with patch("pion_background_chain_audit.histogram",side_effect=AssertionError("parent serialization was too early")) as serializer:
                with self.assertRaisesRegex(ValueError,"accepted_AB_input_identity_mismatch:Center:"+role):
                    self.observe(hist,packet)
                serializer.assert_not_called()
            self.paths[role].write_bytes(("synthetic backing "+role).encode())
            self.assertFalse(list(self.outdir.iterdir()))

    def test_enabled_observer_rejects_configuration_disagreement(self):
        from pion_background_chain_audit import digest as source_digest
        for key,value in (("denom_floor",1.e-4),("clip_min",None),("clip_max",8.)):
            hist,packet=self.runtime_hist()
            hist["_pion_t_bin_parent_results"][1]["final_diagnostic_application_result"]["diagnostics"][key]=value
            before=source_digest(hist)
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,"production_weight_configuration_disagreement"):
                self.observe(hist,packet)
            self.assertEqual(source_digest(hist),before)
            self.assertFalse(list(self.outdir.iterdir()))

    def test_enabled_observer_captures_policy_without_optional_parent_diagnostics(self):
        hist,packet=self.runtime_hist()
        for p in hist["_pion_t_bin_parent_results"]: p["final_diagnostic_application_result"]=None
        result=self.observe(hist,packet)
        self.assertEqual(result["parents"][2]["production_weight_configuration"],
                         {"weight_denominator_floor":2.e-8,"weight_clip_bounds":[.02,7.5]})


class EventTests(unittest.TestCase):
    def test_capture_copies_original_caches(self):
        _,_,hist=fixture();before=deepcopy(hist)
        packet,_,_=fixture()
        self.assertEqual(hist,before)
        self.assertEqual(cap.validate_packet(packet),packet)

    def test_signed_quadratic_sumw2_and_cut_separation(self):
        packet,_,_=fixture();rows=packet["events"]["pion"]
        a=replay.spectrum([r for r in rows if r["t_index"]==0],EDGES)
        self.assertEqual(a["domains"]["normal_full_MM"]["signed"],1.5)
        self.assertEqual(a["domains"]["normal_full_MM"]["sumw2"],10.625)
        self.assertEqual(a["domains"]["signal_1.10_1.16"]["signed"],.75)
        cut=replay.spectrum([r for r in rows if r["t_index"]==0 and r["allcuts"]],EDGES)
        self.assertEqual(cut["domains"]["normal_full_MM"]["sumw2"],5.3125)

    def test_signed_cancellation_is_not_zero_variance(self):
        self.assertEqual(replay.moments([2.,-2.])["signed"],0.)
        self.assertEqual(replay.moments([2.,-2.])["sumw2"],8.)

    def test_declared_synthetic_truth_transfer(self):
        # Analytic synthetic truth only: transfer = .1, no physics bias claim.
        packet,_,_=fixture();rows=[r for r in packet["events"]["pion"] if r["t_index"]==0]
        a=replay.spectrum(rows,EDGES,[0.,.1,.1,.1,.1,0.])
        self.assertAlmostEqual(a["domains"]["normal_full_MM"]["signed"],.15)
        self.assertAlmostEqual(a["domains"]["normal_full_MM"]["sumw2"],.10625)

    def test_canonical_endpoint_and_flow_policy(self):
        self.assertEqual(cap.bin_index(3.,[0.,1.,2.,3.],final_endpoint=True),2)
        self.assertIsNone(cap.bin_index(3.,[0.,1.,2.,3.]))
        self.assertEqual(cap.bin_index(1.,[0.,1.,2.,3.]),1)
        packet,_,_=fixture();r=packet["events"]["pion"][0]
        rows=[{**r,"MM":.6},{**r,"MM":1.3},{**r,"MM":.7}]
        a=replay.spectrum(rows,EDGES,[0.,1.,1.,1.,1.,0.])
        self.assertEqual(a["signed"],[0.,2.,0.,0.,0.,0.])

    def test_parent_child_and_source_closure_with_outside_phi(self):
        packet,_,_=fixture();rows=deepcopy(packet["events"]["pion"])
        rows[0]["phi_index"]=None;rows[1]["phi_index"]=1
        result=replay.split_spectra(rows,EDGES)
        self.assertEqual(result["closure"],"PASS")
        self.assertEqual(result["outside_phi"]["count"][1],1)

    def test_packet_negative_identity_and_selection_gates(self):
        packet,_,_=fixture()
        cases=[(lambda p:p["events"]["pion"].append(p["events"]["pion"][0]),"duplicate"),
               (lambda p:p["events"]["pion"][0].update(coordinate_fingerprint="wrong"),"coordinate"),
               (lambda p:p["events"]["pion"][0].update(post_proton_factor=.5),"proton_on_pion"),
               (lambda p:p["events"]["pion"][0].update(source_sha256="wrong"),"source_file"),
               (lambda p:p["events"]["pion"][0].update(t_index=2),"event_bin"),
               (lambda p:p["events"]["pion"][0].update(MM=float("nan")),"nonfinite"),
               (lambda p:p["events"]["pion"][0].update(fold=99),"fold_identity"),
               (lambda p:p["events"]["pion"].pop(),"missing_source_rows"),
               (lambda p:p.update(setting="Right"),"packet_setting")]
        for mutate,reason in cases:
            p=deepcopy(packet);mutate(p)
            with self.subTest(reason=reason),self.assertRaisesRegex(ValueError,reason):cap.validate_packet(p)

    def test_missing_allcuts_and_changed_hash_rejected(self):
        packet,_,_=fixture()
        for r in packet["events"]["pion"]:r["allcuts"]=False
        packet["event_hashes"]["pion"]=cap.digest(packet["events"]["pion"])
        with self.assertRaisesRegex(ValueError,"missing_allcuts"):cap.validate_packet(packet)
        packet,_,_=fixture();packet["event_hashes"]["kaon"]="wrong"
        with self.assertRaisesRegex(ValueError,"hash_mismatch"):cap.validate_packet(packet)
        packet,_,_=fixture();packet["events"]["kaon"] = deepcopy(packet["events"]["pion"])
        with self.assertRaisesRegex(ValueError,"cross_population_leakage"):cap.validate_packet(packet)

    def test_fold_stability_balance_and_no_false_holdout(self):
        packet,_,_=fixture();a=replay.fold_report(packet["events"])
        reverse={p:list(reversed(rows)) for p,rows in packet["events"].items()}
        self.assertEqual(a,replay.fold_report(reverse))
        self.assertEqual(a["independent_model_evaluation"],"unavailable")
        for population in a["populations"].values():
            for f in population.values():
                self.assertEqual(sum(v["count"] for side in f.values() for v in side.values()),24)

    def test_full_synthetic_replay_and_immutable_inputs(self):
        packets,candidate=comparison_fixture();before=deepcopy((packets,candidate))
        result=replay.replay(packets,candidate,expected_head="a"*40)
        self.assertEqual((packets,candidate),before)
        self.assertEqual(len(result["parents"]),6)
        self.assertTrue(all(len(p["models"])==9 for p in result["parents"]))
        self.assertEqual(result["limits"]["bias"],"bias_truth_unavailable")
        self.assertEqual(result["limits"]["MSE"],"MSE_selection_unavailable")
        self.assertIsNone(result["selected_lambda"])

    def test_replay_provenance_and_model_negative_gates(self):
        cases=[(lambda ps,c:ps["Center"]["parents"][0].update(baseline_parent_fingerprint="wrong"),"parent_model_identity"),
               (lambda ps,c:c["input_metadata"]["Left"].update(canonical_t_edges=[0.,1.,2.,4.]),"canonical_edges"),
               (lambda ps,c:c["provenance"]["Center"].update(sha256="wrong"),"original_AB"),
               (lambda ps,c:c["parents"][0]["C"].pop(),"seven_predeclared"),
               (lambda ps,c:c["parents"][0]["fixed_K_B"].__setitem__(1,5.),"fixed_protected"),
               (lambda ps,c:c["parents"][0]["C"][0]["transfer"]["weights"].__setitem__(1,10.),"weight_ratio"),
               (lambda ps,c:c["parents"][0]["C"][0]["model"].__setitem__(1,0.),"unsupported_denominator")]
        for mutate,reason in cases:
            ps,c=comparison_fixture();mutate(ps,c)
            with self.subTest(reason=reason),self.assertRaisesRegex(ValueError,reason):replay.replay(ps,c,expected_head="a"*40)

    def test_zero_numerator_and_unsupported_event(self):
        packet,_,_=fixture();rows=packet["events"]["pion"]
        a=replay.spectrum(rows,EDGES,[0.]*6)
        self.assertEqual(a["domains"]["normal_full_MM"]["signed"],0.)
        with self.assertRaisesRegex(ValueError,"unsupported_denominator"):
            replay.spectrum(rows,EDGES,[0.,None,0.,0.,0.,0.])

    def test_default_off_and_import_isolation(self):
        with patch.dict(os.environ,{},clear=True):self.assertIsNone(cap.observe({}, {}, [], [], "unused"))
        script="import sys;from testing import pion_local_t_event_validation_capture;assert 'ROOT' not in sys.modules;assert 'numpy' not in sys.modules"
        subprocess.run([sys.executable,"-B","-c",script],cwd=ROOT,check=True)
        tree=ast.parse((ROOT/"src/binning/calculate_yield.py").read_text())
        function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="calculate_yield_data")
        guard=next(n for n in function.body if isinstance(n,ast.If) and cap.ENABLE_ENV in ast.unparse(n.test))
        self.assertIn("observer.observe",ast.unparse(guard))
        self.assertFalse(any(isinstance(n,(ast.Import,ast.ImportFrom)) and "pion_local_t_event_validation" in ast.unparse(n)
                             for n in tree.body))
        original=ast.parse(subprocess.check_output(["git","show","HEAD:src/binning/calculate_yield.py"],cwd=ROOT).decode())
        function.body.remove(guard)
        self.assertEqual(ast.dump(tree),ast.dump(original))

    def test_post_proton_factor_is_quadratic_and_rf_zero_is_preserved(self):
        packet,_,_=fixture();row=deepcopy(packet["events"]["kaon"][0])
        row.update(post_proton_factor=0.,proton_rf_accept=False)
        result=replay.spectrum([row],EDGES)
        self.assertEqual(result["domains"]["normal_full_MM"]["sumw2"],0.)
        row["post_proton_factor"]=.25
        result=replay.spectrum([row],EDGES)
        self.assertEqual(result["domains"]["normal_full_MM"]["signed"],.5)
        self.assertEqual(result["domains"]["normal_full_MM"]["sumw2"],.25)

    def test_capture_rejects_cache_disagreement_without_mutation(self):
        _,_,hist=fixture();hist["_yield_sub_event_cache"]["prompt"]["adj_MM"][0]=.9
        before=deepcopy(hist)
        packet,_,_=fixture()
        with self.assertRaisesRegex(ValueError,"child_cache_row_mismatch"):
            cap.build_packet(hist,{},packet["t_edges"],packet["phi_edges"],source_head=packet["source_HEAD"],
                input_identity=packet["input_identity"],parent_records=packet["parents"],source_files={})
        self.assertEqual(hist,before)

    def test_production_convention_floor_reported_without_changing_C(self):
        packets,candidate=comparison_fixture()
        for packet in packets.values():
            for p in packet["parents"]:p["production_weight_configuration"]["weight_denominator_floor"]=1.
        result=replay.replay(packets,candidate,expected_head="a"*40)
        model=result["parents"][0]["models"]["C_0"]
        self.assertTrue(model["convention_differs"])
        self.assertGreater(model["application"]["allcuts"]["parent"]["domains"]["signal_1.10_1.16"]["sumw2"],0)
        self.assertEqual(model["production_convention_application"]["allcuts"]["parent"]["domains"]["signal_1.10_1.16"]["sumw2"],0)

    def test_pdf_plot_branch_synthetic_layout_only(self):
        packets,candidate=comparison_fixture()
        result=replay.replay(packets,candidate,expected_head="a"*40)
        # No PDF is created here: inspect plotting commands and page layouts.
        from matplotlib.backends import backend_pdf
        pages=[]
        class Pages:
            def __init__(self,*args,**kwargs):pass
            def __enter__(self):return self
            def __exit__(self,*args):pass
            def savefig(self,fig):
                fig.canvas.draw()
                pages.append((len(fig.axes),fig._suptitle.get_text()))
        with patch.object(backend_pdf,"PdfPages",Pages):replay.render_pdf(result,"unused.pdf")
        self.assertEqual(len(pages),12)
        self.assertTrue(all(n==4 for n,_ in pages))

    def test_strict_serialization(self):
        for raw in ('{"x":NaN}','{"x":1e999}','{"x":1,"x":2}'):
            with self.assertRaises(ValueError):replay.strict_json(raw)
        packet,_,_=fixture();self.assertEqual(cap.encoded(packet),cap.encoded(replay.strict_json(cap.encoded(packet))))

    @unittest.skipUnless(os.environ.get("KAONLT_ACCEPTED_C_JSON"),"original saved C JSON not supplied")
    def test_supplied_saved_C_schema_without_refitting(self):
        replay.validate_candidate(replay.strict_json(Path(os.environ["KAONLT_ACCEPTED_C_JSON"]).read_bytes()))

    @unittest.skipUnless(os.environ.get("KAONLT_ACCEPTED_AB_ZIP"),"accepted A/B archive not supplied")
    def test_capture_input_pins_match_original_accepted_AB_artifacts(self):
        path=Path(os.environ["KAONLT_ACCEPTED_AB_ZIP"])
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         "8bd71f9ab17fa456a86911b8b43eed9c2a1a9da3101cbbff776c56f66fa901ee")
        with zipfile.ZipFile(path) as archive:
            for setting in ("Center","Left"):
                raw=archive.read("global/"+setting+"_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json")
                self.assertEqual(hashlib.sha256(raw).hexdigest(),cap.ACCEPTED_AB_JSON_SHA256[setting])
                doc=replay.strict_json(raw)
                self.assertEqual(doc["source_HEAD"],replay.c.SOURCE)
                inputs={r:{k:v[k] for k in ("sha256","bytes")} for r,v in doc["original_analysis_input_identity"].items()}
                cap.verify_accepted_input_identity(inputs,setting)


if __name__ == "__main__":unittest.main()
