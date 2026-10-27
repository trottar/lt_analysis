"""ROOT-free replay of captured signed events using the accepted saved C sweep.

No fit, lambda selection, cross-section extraction or production application.
"""
from __future__ import annotations

import argparse
import bisect
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

try:
    from . import pion_local_t_event_validation_capture as capture
    from . import pion_local_t_transfer_variance_comparison as c
    from . import pion_local_t_statistical_validation as statistics
except ImportError:
    import pion_local_t_event_validation_capture as capture
    import pion_local_t_transfer_variance_comparison as c
    import pion_local_t_statistical_validation as statistics

require, digest = capture.require, capture.digest
DOMAINS = {"normal_full_MM": None, "full_0.7_1.3": (.7,1.3), "signal_1.10_1.16": (1.10,1.16),
           "visual_1.105_1.135": (1.105,1.135), "troublesome_0.94_1.03": (.94,1.03)}
ACCEPTED_AB = {"Center": "7e34ec134ef28378ffebf1686c47915a266478531d9910bbdd211f5392cb1f24",
               "Left": "f5c8cf2483ce77d74dc2e2a497d8b7c5f1b0aa3c96b4d76c3964d813e675b812"}
LIMITATIONS = {
    "independent_fitted_estimator_evaluation": {"status": "unavailable",
        "first_missing_step": "training-only staged/protected kaon-side fits and priors; saved K_B and stages use full-data events",
        "dependencies": ["K_B", "staged amplitudes", "staged priors", "protected signal model", "local SIMC templates"],
        "local_SIMC": "independent of data folds, but no retraining/independent template-uncertainty ensemble supplied",
        "fold_claim": "disjoint record partitions only; no held-out model scores"},
    "bias": "bias_truth_unavailable", "MSE": "MSE_selection_unavailable",
    "truth_audit": "Existing pion-control SIMC shapes describe selected pion channels, not independently truth-tagged kaon-PID leakage with source normalization. No leakage truth is supplied.",
    "uncertainty": "fixed-weight Sumw2 only; no fitted-parameter covariance or repeated-training variance",
    "yield": "allcuts event window integral only; not authoritative extracted yield/error or absolute SIMC comparison",
    "cross_population": "tree-local entry identities cannot prove physical-event independence between different PID trees",
    "selection": "none; all seven saved strengths retained"}


def strict_json(raw):
    def pairs(items):
        result = {}
        for key,value in items:
            require(key not in result, "duplicate_json_key")
            result[key] = value
        return result
    def invalid(value):
        raise ValueError("nonfinite_json:" + value)
    result = json.loads(raw, object_pairs_hook=pairs, parse_constant=invalid)
    capture.encoded(result)  # Includes overflow such as 1e999.
    return result


def moments(weights, *, include_tails=False):
    positive = math.fsum(w for w in weights if w > 0)
    negative = math.fsum(w for w in weights if w < 0)
    signed = math.fsum(weights)
    absolute = math.fsum(abs(w) for w in weights)
    sumw2 = math.fsum(w*w for w in weights)
    result = {"count": len(weights), "signed": signed, "absolute": absolute, "positive": positive,
              "negative": negative, "sumw2": sumw2, "fixed_weight_error": math.sqrt(sumw2),
              "effective_absolute_support": absolute**2/sumw2 if sumw2 else 0.,
              "effective_signed_support": signed**2/sumw2 if sumw2 else 0.,
              "maximum_absolute_contribution": max(map(abs, weights), default=0.)}
    if not include_tails:
        return result
    ordered = sorted(abs(w) for w in weights)
    def quantile(fraction):
        if not ordered:
            return 0.
        index = (len(ordered)-1)*fraction
        lower = int(index)
        return ordered[lower] + (ordered[min(lower+1, len(ordered)-1)]-ordered[lower])*(index-lower)
    largest = ordered[-max(1, math.ceil(.01*len(ordered))):] if ordered else []
    result.update(absolute_event_contribution_quantiles={str(p): quantile(p) for p in (.5, .9, .95, .99, .999)},
            largest_one_percent_absolute_fraction=math.fsum(largest)/absolute if absolute else 0.,
            largest_one_percent_Sumw2_fraction=math.fsum(w*w for w in largest)/sumw2 if sumw2 else 0.,
            tail_policy="absolute signed-source contributions; linear unweighted event quantiles; ceil(1% of records)")
    return result


def spectrum(rows, edges, transfer=None, *, include_tails=False):
    """ROOT bin indexing; flows explicitly retained, transfer flow weights zero."""
    bins = [[] for _ in range(len(edges)+1)]
    for r in rows:
        i = bisect.bisect_right(edges, r["MM"])
        w = r["coefficient"] * r["post_proton_factor"]
        if transfer is not None:
            require(transfer[i] is not None or w == 0, "unsupported_denominator_on_event:MM="+str(r["MM"])+":bin="+str(i))
            w *= transfer[i] if transfer[i] is not None else 0.
        require(math.isfinite(w), "nonfinite_application")
        bins[i].append(w)
    cells = [moments(weights, include_tails=include_tails) for weights in bins]
    result = {k:[cell[k] for cell in cells] for k in cells[0]}
    result["domains"] = {}
    centers = [(a+b)/2 for a,b in zip(edges,edges[1:])]
    for name, window in DOMAINS.items():
        ix = [i+1 for i,x in enumerate(centers) if window is None or window[0] <= x <= window[1]]
        result["domains"][name] = moments([w for i in ix for w in bins[i]], include_tails=include_tails)
    result["window_policy"] = "inclusive normal-bin centers; flow bins excluded from window totals"
    return result


def split_spectra(rows, edges, weights=None, *, include_tails=False):
    parent = spectrum(rows, edges, weights, include_tails=include_tails)
    sources = {s: spectrum([r for r in rows if r["source_label"] == s], edges, weights, include_tails=include_tails) for s in capture.SOURCES}
    children = {str(pi):spectrum([r for r in rows if r["phi_index"] == pi], edges, weights, include_tails=include_tails)
                for pi in sorted({r["phi_index"] for r in rows if r["phi_index"] is not None})}
    outside = spectrum([r for r in rows if r["phi_index"] is None], edges, weights, include_tails=include_tails)
    for parts in (list(sources.values()), [*children.values(),outside]):
        for key in ("signed", "absolute", "positive", "negative", "sumw2", "count"):
            require(all(math.isclose(v, math.fsum(p[key][i] for p in parts), rel_tol=1e-11, abs_tol=1e-12)
                        for i,v in enumerate(parent[key])), "parent_child_or_source_closure")
    return {"parent": parent, "sources": sources, "children": children, "outside_phi": outside, "closure": "PASS"}


def fold_report(events):
    result = {"policy": deepcopy(capture.FOLD_POLICY), "independent_model_evaluation": "unavailable", "populations": {}}
    for population,rows in events.items():
        folds = {}
        for f in range(capture.FOLD_POLICY["count"]):
            test = [r for r in rows if r["fold"] == f]
            train = [r for r in rows if r["fold"] != f]
            identity = lambda r: (r["source_sha256"],r["source_label"],r["entry_index"])
            require(not set(map(identity,train)) & set(map(identity,test)), "fold_leakage")
            folds[str(f)] = {role: {s:moments([r["coefficient"]*r["post_proton_factor"] for r in subset if r["source_label"] == s])
                                   for s in capture.SOURCES} for role,subset in (("training",train),("evaluation",test))}
        result["populations"][population] = folds
    return result


def validate_candidate(candidate, *, matched_authority=None):
    require(candidate.get("schema_version") == c.SCHEMA and candidate.get("status") == "available"
            and candidate.get("inputs_unchanged") is True and candidate.get("diagnostic_only") is True
            and candidate.get("production_applied") is False and candidate.get("scientifically_accepted") is False,
            "accepted_C_schema_or_authority")
    require(candidate.get("source_HEAD_of_input") == (c.SOURCE if matched_authority is None else matched_authority['source_HEAD'])
            and digest(candidate.get("policy")) == digest(c.POLICY), "accepted_C_policy")
    require(len(candidate.get("parents",[])) == 6
            and {(p["setting"],p["t_bin"]) for p in candidate["parents"]} == {(s,t) for s in c.SETTINGS for t in (1,2,3)}, "six_C_parents")
    if matched_authority is None:
        require('cohort_authority' not in candidate, 'matched_candidate_requires_explicit_mode')
        for s in c.SETTINGS:
            require(candidate["provenance"][s]["sha256"] == ACCEPTED_AB[s], "original_AB_identity")
    else:
        require(candidate.get('cohort_authority') == matched_authority, 'matched_C_authority_identity')
    for row in candidate["parents"]:
        c.validate_parent(row["A_B_unchanged"])
        require([v["lambda"] for v in row["C"]] == list(c.LAMBDAS), "seven_predeclared_strengths")
        require(all(v.get("status") == "available" and v.get("selected") is False for v in row["C"]), "unavailable_or_selected_C")
        if matched_authority is not None:
            for fit in row['C']:
                expected = c.matched_fit_identity(fit, matched_authority, candidate['provenance'], row['setting'],
                                                   row['A_B_unchanged'], candidate['implementation_sha256'])
                require(fit.get('fit_identity') == expected and fit.get('model_identity') ==
                        digest([expected, fit['model'], row['fixed_K_B']]), 'matched_fit_or_model_identity')
        edges = row["A_B_unchanged"]["frozen_inputs"]["pion"]["edges"]
        k = row["fixed_K_B"]
        require(k == row["A_B_unchanged"]["candidate_B"]["H_kaon_pion_bg_fit_total_refined"]["values"], "fixed_protected_K_B")
        for model, numerator, transfer in [(v["model"],k,v["transfer"]) for v in row["C"]] + [
            (row["A_B_unchanged"]["baseline_A"]["H_pion_fit_total_refined"]["values"],
             row["A_B_unchanged"]["baseline_A"]["H_kaon_pion_bg_fit_total_refined"]["values"],row["A_transfer"]),
            (row["A_B_unchanged"]["candidate_B"]["H_pion_fit_total_refined"]["values"],k,row["B_transfer"])]:
            require(len(model) == len(k) == len(transfer["weights"]) == len(edges)+1, "model_bin_layout")
            for i in range(1,len(edges)):
                require(all(math.isfinite(x) and x >= 0 for x in (model[i],numerator[i])), "model_nonfinite_or_negative")
                require(numerator[i] == 0 or model[i] > 0,
                        f'unsupported_denominator:{row["setting"]}/t{row["t_bin"]}:bin={i}:MM={edges[i-1:i+1]}')
                expected = numerator[i]/model[i] if model[i] > 0 else 0.
                actual = transfer["weights"][i]
                require(actual is not None and math.isfinite(actual) and math.isclose(expected,actual,rel_tol=1e-12,abs_tol=1e-14), "saved_weight_ratio_mismatch")
    return candidate


def replay(packets, candidate, *, expected_head, matched_authority=None, include_statistics=False):
    before = digest([packets,candidate])
    detailed = matched_authority is not None
    validate_candidate(candidate, matched_authority=matched_authority)
    require(set(packets) == set(c.SETTINGS), "both_settings_required")
    output = {"schema_version": "pion_local_t_event_validation_replay/v1", "diagnostic_only": True,
              "production_applied": False, "scientifically_accepted": False, "selected_lambda": None,
              "candidate_sha256_of_canonical_json": digest(candidate), "limits": deepcopy(LIMITATIONS), "parents": [], "folds": {}}
    for s,packet in packets.items():
        capture.validate_packet(packet)
        require(packet["source_HEAD"] == expected_head, "capture_source_HEAD_mismatch")
        require(packet["setting"] == s, "wrong_setting")
        if matched_authority is not None:
            auth = packet.get('cohort_authority', {})
            expected_inputs = matched_authority['inputs'][s]
            if matched_authority.get('private_input_policy') == 'frozen_shiftprep_private_copy':
                require(auth.get('original_input_identity') == expected_inputs, 'matched_replay_original_inputs')
                expected_inputs = candidate['provenance'][s]['runtime_input_identity']
                require(auth.get('runtime_input_identity') == expected_inputs, 'matched_replay_runtime_inputs')
            require(auth.get('mode') == 'matched_cohort' and auth.get('manifest_sha256') == digest(matched_authority)
                    and auth.get('attempt_id') == matched_authority['attempt_id']
                    and packet['input_identity'] == expected_inputs, 'matched_replay_authority')
        metadata = candidate["input_metadata"][s]
        require(packet["t_edges"] == metadata["canonical_t_edges"] and packet["phi_edges"] == metadata["canonical_phi_edges"], "canonical_edges_mismatch")
        for role in ("data", "dummy"):
            require(all(packet["input_identity"][role][k] == metadata["original_analysis_input_identity"][role][k]
                        for k in ("sha256","bytes")), "original_source_data_mismatch")
        for path,identity in metadata["source_files"].items():
            require(all(packet["source_files"][path][k] == identity[k] for k in ("sha256","bytes")), "frozen_model_source_mismatch:"+path)
        output["folds"][s] = fold_report(packet["events"])
        for row in [r for r in candidate["parents"] if r["setting"] == s]:
            ti = row["t_bin"]-1
            p = row["A_B_unchanged"]
            parent = packet["parents"][ti]
            for key in ("t_bin_index", "t_edges", "pion_parent_id", "canonical_interval_pair_id", "canonical_interval_pair_hash", "coordinate_fingerprint", "baseline_parent_fingerprint"):
                require(parent[key] == p[key], "parent_model_identity_mismatch:"+key)
            require(packet["coordinate_fingerprint"] == p["coordinate_fingerprint"], "coordinate_mismatch")
            for side, key in (("pion","pion_input"),("kaon","kaon_input")):
                require(parent[key] == p["frozen_inputs"][side], "frozen_parent_input_mismatch")
            edges = p["frozen_inputs"]["pion"]["edges"]
            record = {"setting": s, "t_bin": ti+1, "edges": edges, "populations": {}, "models": {}, "support_failures": []}
            if detailed:
                record['pion_fit_models'] = {'A': p['baseline_A']['H_pion_fit_total_refined']['values'],
                                             'B': p['candidate_B']['H_pion_fit_total_refined']['values']}
                record['pion_fit_models'].update({'C_' + format(v['lambda'], 'g'): v['model'] for v in row['C']})
                record['kaon_fit_models'] = {'A': p['baseline_A']['H_kaon_pion_bg_fit_total_refined']['values'],
                                             'B_C': row['fixed_K_B']}
                fit_bins = c.problem(p)['ix'].tolist() if matched_authority is not None else []
                if fit_bins:
                    observed_fit = p['frozen_inputs']['pion']
                    record['common_fit_bins'] = fit_bins
                    record['descriptive_data_chi2'] = {
                        name: math.fsum((observed_fit['values'][i] - model[i])**2 / observed_fit['stored_sumw2'][i]
                                        for i in fit_bins) for name, model in record['pion_fit_models'].items()}
            for population in ("pion","kaon"):
                rows = [r for r in packet["events"][population] if r["t_index"] == ti]
                record["populations"][population] = {cut:split_spectra([r for r in rows if r[cut]],edges,include_tails=detailed) for cut in ("nommcuts","allcuts")}
            # The pion parent uses all phi, including explicit outside-acceptance rows.
            control = record["populations"]["pion"]["nommcuts"]["parent"]
            reference = p["frozen_inputs"]["pion"]
            for key,stored in (("signed","values"),("sumw2","stored_sumw2")):
                require(all(math.isclose(x,y,rel_tol=1e-9,abs_tol=1e-12) for x,y in zip(control[key],reference[stored])), "pion_input_selection_closure")
            variants = [("A",row["A_transfer"],p["baseline_A"]),("B",row["B_transfer"],p["candidate_B"])]
            variants += [("C_"+format(v["lambda"],"g"),v["transfer"],v) for v in row["C"]]
            for name,transfer,fit in variants:
                # Saved ratio is diagnostic; actual production convention may differ.
                ratio = [0., *transfer["weights"][1:-1], 0.]
                model = fit["model"] if name.startswith("C_") else fit["H_pion_fit_total_refined"]["values"]
                numerator = row["fixed_K_B"] if name != "A" else p["baseline_A"]["H_kaon_pion_bg_fit_total_refined"]["values"]
                cfg = parent["production_weight_configuration"]
                require(all(math.isfinite(x) for x in ratio), "unsupported_positive_transfer")
                convention = cfg.get("weight_denominator_floor")
                require(convention is not None and math.isfinite(convention) and convention >= 0, "production_weight_floor_missing")
                physical = [0.] + [0. if model[i] <= convention else numerator[i]/model[i] for i in range(1,len(edges))] + [0.]
                # Configured clipping is preserved, never introduced for C risk.
                require(len(cfg.get("weight_clip_bounds", [])) == 2, "production_clip_policy_missing")
                for bound,operation in zip(cfg["weight_clip_bounds"],(max,min)):
                    if bound is not None:
                        require(math.isfinite(bound), "nonfinite_clip_policy")
                        physical[1:-1] = [operation(bound,v) for v in physical[1:-1]]
                rows = [r for r in packet["events"]["pion"] if r["t_index"] == ti]
                record["models"][name] = {"maximum_ratio": max(ratio), "weights": ratio,
                    "production_convention_weights": physical,
                    "convention_differs": any(not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14) for a,b in zip(ratio,physical)),
                    "numerical_fit": "saved full-sample model; no fit performed",
                    "fit_quality": fit.get("data_chi2_ndf", "stored A/B diagnostics in original C input"),
                    "shape_deterioration": fit.get("shape_fit_worsened", "see original A/B diagnostics"),
                    "application": {cut:split_spectra([r for r in rows if r[cut]],edges,ratio,include_tails=detailed) for cut in ("nommcuts","allcuts")},
                    "production_convention_application": {cut:split_spectra([r for r in rows if r[cut]],edges,physical,include_tails=detailed) for cut in ("nommcuts","allcuts")}}
            if detailed:
                for name, model in record['models'].items():
                    reference = record['models']['B']['application']
                    model['paired_difference_from_B'] = {}
                    for label, key in (('ratio', 'weights'), ('production_convention', 'production_convention_weights')):
                        delta = [a-b for a, b in zip(model[key], record['models']['B'][key])]
                        model['paired_difference_from_B'][label] = {
                            cut: split_spectra([r for r in rows if r[cut]], edges, delta, include_tails=True)
                            for cut in ('nommcuts', 'allcuts')}
                    model['paired_difference_uncertainty'] = 'same-event difference weights squared; fitted-parameter uncertainty unavailable'
                    if matched_authority is not None:
                        model['numerical_fit'] = ('fresh same-cohort C sweep upstream; this replay performs no fit'
                                                  if name.startswith('C_') else 'fresh same-run A/B upstream; this replay performs no fit')
                    model['descriptive_flags'] = {
                        'fit_degradation': model['shape_deterioration'],
                        'variance_amplification_domains': {cut: [domain for domain in DOMAINS
                            if model['application'][cut]['parent']['domains'][domain]['sumw2'] >
                               reference[cut]['parent']['domains'][domain]['sumw2'] and not math.isclose(
                                   model['application'][cut]['parent']['domains'][domain]['sumw2'],
                                   reference[cut]['parent']['domains'][domain]['sumw2'], rel_tol=1e-10, abs_tol=1e-14)]
                            for cut in ('nommcuts', 'allcuts')},
                        'interpretation': 'descriptive full-sample comparisons only; no accepted or selected strength'}
            if detailed or include_statistics:
                record["training_dependency_audit"] = statistics.training_dependency_audit(p)
                record["conditional_fixed_weight_covariance"] = statistics.event_replay_statistics(record, DOMAINS)
            output["parents"].append(record)
    require(digest([packets,candidate]) == before, "replay_input_mutation")
    output["inputs_unchanged"] = True
    if matched_authority is not None:
        output['cohort_authority_sha256'] = digest(matched_authority)
    return output


def render_pdf(result, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(path, metadata={"Title":"Detached signed-event replay; NO selected strength", "CreationDate":None}) as pdf:
        for p in result["parents"]:
            fig, axes = plt.subplots(2,2,figsize=(12,9),layout="constrained")
            fig.suptitle(f'{p["setting"]}/lowe t{p["t_bin"]} — full-sample signed-event replay\nNO held-out validation; bias/MSE unavailable; NO selected strength',fontsize=12)
            x = [(a+b)/2 for a,b in zip(p["edges"],p["edges"][1:])]
            for row,cut in enumerate(("nommcuts","allcuts")):
                for col,window in enumerate(((p["edges"][0],p["edges"][-1]),(1.10,1.16))):
                    ax = axes[row,col]
                    observed = p["populations"]["kaon"][cut]["parent"]
                    ax.errorbar(x,observed["signed"][1:-1],observed["fixed_weight_error"][1:-1],fmt="k.",ms=2,label="kaon post-proton cache")
                    for name,model in p["models"].items():
                        a = model["application"][cut]["parent"]
                        v = a["domains"]["signal_1.10_1.16"]["sumw2"]
                        ax.plot(x,a["signed"][1:-1],label=f'{name}; window Sumw2={v:.3g}',lw=1)
                    ax.set(xlim=window,title=cut+(" full MM" if col == 0 else " Lambda"),xlabel="shifted MM (GeV)",ylabel="signed transferred bin content")
                    visible = [float(y) for line in ax.lines for xx,y in zip(line.get_xdata(),line.get_ydata()) if window[0] <= xx <= window[1]]
                    visible += [y+sign*e for xx,y,e in zip(x,observed["signed"][1:-1],observed["fixed_weight_error"][1:-1])
                                if window[0] <= xx <= window[1] for sign in (-1,1)]
                    if visible:
                        lo,hi=min(0.,min(visible)),max(0.,max(visible)); pad=max((hi-lo)*.08,1e-12);ax.set_ylim(lo-pad,hi+pad)
            axes[0,0].legend(fontsize=6,ncol=2)
            pdf.savefig(fig);plt.close(fig)
            fig, axes = plt.subplots(2,2,figsize=(12,9),layout="constrained")
            fig.suptitle(f'{p["setting"]}/lowe t{p["t_bin"]} — detached support and fixed-weight variance\nNo selected strength; no physical bias truth',fontsize=12)
            for ax,cut,window in ((axes[0,0],"nommcuts",(.94,1.03)),(axes[0,1],"allcuts",(1.105,1.135))):
                for name,model in p["models"].items():
                    ax.plot(x,model["application"][cut]["parent"]["signed"][1:-1],label=name)
                ax.set(xlim=window,title=cut,xlabel="shifted MM (GeV)",ylabel="signed transferred bin content")
                visible=[float(y) for line in ax.lines for xx,y in zip(line.get_xdata(),line.get_ydata()) if window[0]<=xx<=window[1]]
                if visible:
                    lo,hi=min(0.,min(visible)),max(0.,max(visible));pad=max((hi-lo)*.08,1e-12);ax.set_ylim(lo-pad,hi+pad)
            names=list(p["models"])
            for cut,style in (("nommcuts","o-"),("allcuts","x--")):
                for domain in ("normal_full_MM","signal_1.10_1.16"):
                    axes[1,0].plot(range(len(names)),[p["models"][n]["application"][cut]["parent"]["domains"][domain]["sumw2"] for n in names],style,label=cut+" "+domain)
            axes[1,0].set_xticks(range(len(names)),names,rotation=45,fontsize=7)
            axes[1,0].set_yscale("symlog",linthresh=1e-12);axes[1,0].set_ylabel("fixed-weight Sumw2 (not training variance)")
            axes[1,0].legend(fontsize=6)
            lines=["All source/child sum closures: PASS", "Unsupported positive denominator: none admitted", "Phi children: "+str(list(p["populations"]["pion"]["nommcuts"]["children"])),
                   "All strengths from saved full-sample fits:"]
            for name,model in p["models"].items():
                lines.append(f'{name}: chi2/ndf={model["fit_quality"] if isinstance(model["fit_quality"],(int,float)) else "original A/B"}; max ratio={model["maximum_ratio"]:.4g}; convention differs={model["convention_differs"]}')
            lines += ["Held-out fitted estimator: unavailable", "First missing step: training-only staged/protected fits/priors", "bias_truth_unavailable; MSE_selection_unavailable", "No covariance or repeated-training variance inferred"]
            axes[1,1].axis("off");axes[1,1].text(0,1,"\n".join(lines),va="top",fontsize=8)
            axes[0,0].legend(fontsize=6,ncol=3)
            pdf.savefig(fig);plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--center-packet",type=Path,required=True)
    parser.add_argument("--left-packet",type=Path,required=True)
    parser.add_argument("--candidate-json",type=Path,required=True)
    parser.add_argument("--candidate-sha256",required=True,help="independently reviewed original c_review_candidate_c.json hash")
    parser.add_argument("--capture-source-commit",required=True,help="independently reviewed published capture HEAD")
    parser.add_argument("--statistical-audit",action="store_true",help="Add conditional same-event covariance; never held-out model validation")
    parser.add_argument("--output-json",type=Path,required=True)
    parser.add_argument("--output-pdf",type=Path,required=True)
    args = parser.parse_args(argv)
    inputs = (args.center_packet,args.left_packet,args.candidate_json)
    outputs = (args.output_json,args.output_pdf)
    require(len({p.resolve() for p in (*inputs,*outputs)}) == 5, "path_collision")
    require(all(not p.exists() and not p.is_symlink() for p in outputs), "output_exists")
    require(args.candidate_json.name == "c_review_candidate_c.json", "original_candidate_basename")
    raw = args.candidate_json.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == args.candidate_sha256, "accepted_candidate_hash_mismatch")
    candidate = strict_json(raw)
    packets = {s:strict_json(p.read_bytes()) for s,p in zip(c.SETTINGS,inputs[:2])}
    result = replay(packets,candidate,expected_head=args.capture_source_commit,include_statistics=args.statistical_audit)
    result["input_files"] = {str(p):capture.file_identity(p) for p in inputs}
    render_pdf(result,args.output_pdf)
    result["pdf"] = capture.file_identity(args.output_pdf)
    with args.output_json.open("xb") as stream:
        stream.write(capture.encoded(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
