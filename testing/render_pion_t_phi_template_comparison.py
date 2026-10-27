"""Offline plots of captured detached child fits; never fits or loads ROOT."""
from __future__ import annotations

import argparse
from bisect import bisect_right
import json
import math
from pathlib import Path


def window_sum(values, edges, low, high):
    """Same inclusive normal-bin convention as utility.integrate_hist_range.

    Captured histogram sums are diagnostics, never extracted kaon yields.
    """
    first=max(1,bisect_right(edges,low))
    last=min(len(edges)-1,bisect_right(edges,high))
    return math.fsum(values[first:last+1])


def render(document, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np
    if document.get("schema_version") != "pion_t_phi_template_comparison/v1":
        raise ValueError("unsupported child comparison schema")
    if document.get("status") != "captured_pending_independent_validation" or not document.get("production_parents_unchanged"):
        raise ValueError("failed capture is not scientific evidence")
    label = document.get("data_class", "genuine selected runtime observations; acceptance pending")
    setting = document["phi_setting"]
    window=document["signal_window"]
    colors = dict(Global="tab:blue", Local_t="tab:orange", Local_t_phi="tab:green")
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11})
    output = Path(output)
    if output.exists():
        raise ValueError("PDF destination exists")
    with PdfPages(output) as pdf:
        for cell in document["cells"]:
            fig, axes = plt.subplots(3, 2, figsize=(13, 12), constrained_layout=True)
            fit = cell.get("fit") or {}
            def curve(ax, h, title, color="tab:green", data=False):
                if not isinstance(h, dict) or h.get("status") != "available":
                    return
                edges = h["edges"]
                x = (np.array(edges[:-1])+np.array(edges[1:]))/2
                if data:
                    ax.errorbar(x,h["values"][1:-1],yerr=h["errors"][1:-1],fmt=".",ms=3,color="black",label=title)
                else:
                    ax.stairs(h["values"][1:-1],edges,label=title,color=color)
            pion,kaon = cell["inputs"]["pion"],cell["inputs"]["kaon"]
            curve(axes[0,0],pion,"observed pion control",data=True)
            curve(axes[0,0],fit.get("H_pion_fit_total_refined") or fit.get("H_pion_fit_total"),"child control fit")
            for ax in (axes[0,1],axes[1,1],axes[2,1]):
                curve(ax,kaon,"observed kaon, post proton/RF",data=True)
                curve(ax,fit.get("H_kaon_fit_total_refined") or fit.get("H_kaon_fit_total"),"child full-fit candidate")
                for c,col in zip(("pi_n","pi_delta","pi_sidis"),("tab:purple","tab:red","tab:brown")):
                    curve(ax,fit.get("H_kaon_fit_"+c+"_scaled_refined") or fit.get("H_kaon_fit_"+c+"_scaled"),c,col)
            for side,col in (("pion","tab:green"),("kaon","tab:purple")):
                curve(axes[1,0],fit.get("H_fit_residual_"+side+"_refined"),side+" child fit residual",col)
            axes[1,0].axhline(0,color="grey",lw=.6)
            unavailable=[]
            for name,app in cell["applications"].items():
                if app["status"] != "available":
                    unavailable.append(name+": "+app.get("reason", "unavailable"))
                    continue
                spec = app["cuts"]["nommcuts"]
                x=(np.array(pion["edges"][:-1])+np.array(pion["edges"][1:]))/2
                axes[2,0].stairs(spec["values"][1:-1],pion["edges"],color=colors[name],label=name+" event-applied background")
                axes[2,0].fill_between(x,np.array(spec["values"][1:-1])-np.sqrt(spec["variances"][1:-1]),
                    np.array(spec["values"][1:-1])+np.sqrt(spec["variances"][1:-1]),color=colors[name],alpha=.12)
            axes[1,1].set_xlim(1.08,1.16)
            axes[2,1].set_xlim(1.18,1.5)
            for ax,lo,hi in ((axes[1,1],1.08,1.16),(axes[2,1],1.18,1.5)):
                yy=[]
                for h in [kaon,fit.get("H_kaon_fit_total_refined") or fit.get("H_kaon_fit_total")]:
                    if isinstance(h,dict) and h.get("status")=="available":
                        for i,(left,right) in enumerate(zip(h["edges"],h["edges"][1:]),1):
                            if right>=lo and left<=hi:
                                yy.extend([h["values"][i]-h["errors"][i],h["values"][i]+h["errors"][i]])
                if yy:
                    bottom=min(0.,min(yy)); top=max(yy)
                    ax.set_ylim(bottom,top+max((top-bottom)*.12,1e-12))
            axes[1,1].axvspan(*window,color="grey",alpha=.08,label="configured signal window")
            titles=("Pion control: independent child fit", "Kaon full MM: candidate model and pion components",
                    "Stored child fit residuals", "Lambda window: fitted candidate, protection status retained",
                    "Three-way signed transferred background; bands = fixed-weight Sumw2", "High-MM shoulder (Center t3 focus)")
            for ax,title in zip(axes.flat,titles):
                ax.set_title(title)
                ax.set_xlabel("Missing mass (GeV)")
                ax.set_ylabel("Signed normalized content / bin")
                ax.grid(alpha=.15)
                handles,_=ax.get_legend_handles_labels()
                if handles:
                    ax.legend(fontsize=7)
            state="; ".join(c+":"+s for c,s in cell["template_scopes"].items())
            fig.suptitle(f"{setting} t{cell['t_index']+1} phi{cell['phi_index']+1}: {cell['status']}\n{label}\n{state}",fontsize=12)
            estimates=[]
            for name,app in cell["applications"].items():
                if app["status"]=="available":
                    spec=app["cuts"]["nommcuts"]
                    estimates.append(f"{name}={window_sum(spec['values'],pion['edges'],*window):.4g} +/- {math.sqrt(window_sum(spec['variances'],pion['edges'],*window)):.3g}")
            protected=((fit.get("diagnostics") or {}).get("kaon") or {}).get("pi_delta_signal_protected_fit") or {}
            note="A/B parent control curves have no child projection; not overlaid as child fits. No production subtraction."
            note+=f"\nSignal-window background ({window[0]:g}-{window[1]:g}, inclusive ROOT bins; fixed-weight error): "+"; ".join(estimates)
            note+=f"\nProtected pi_delta: accepted={protected.get('accepted','unavailable')}, status={protected.get('status','unavailable')}; kaon fit status={fit.get('fit_status_kaon','unavailable')}"
            if unavailable:
                note += "\nUnavailable applications:\n"+"\n".join(s[:110] for s in unavailable)
            if cell.get("reason"):
                note += "\nChild fit unavailable: "+cell["reason"][:110]
            fig.supxlabel(note,fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)
        # Numerical support and fit-quality maps show missing cells as masked,
        # never as successful zero-valued fits.
        for group in ("support / fit quality",):
            fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
            keys = ("pion neff","kaon neff","local component count","chi2_ndf_pion","chi2_ndf_kaon","A_sidis")
            for ax,key in zip(axes.flat,keys):
                values=np.full((3,9),np.nan)
                for cell in document["cells"]:
                    fit=cell.get("fit") or {}
                    if key.endswith("neff"):
                        h=cell["inputs"][key.split()[0]]
                        den=sum(h["stored_sumw2"][1:-1]); num=sum(h["values"][1:-1])
                        value=num*num/den if den>0 else None
                    elif key=="local component count":
                        value=sum(s=="local_t_phi" for s in cell["template_scopes"].values())
                    else:
                        value=fit.get(key)
                    if value is not None:
                        values[cell["t_index"],cell["phi_index"]]=value
                im=ax.imshow(np.ma.masked_invalid(values),aspect="auto",origin="lower")
                ax.set(title=key,xlabel="phi index (one based)",ylabel="t index (one based)")
                ax.set_xticks(range(9),range(1,10)); ax.set_yticks(range(3),range(1,4))
                fig.colorbar(im,ax=ax,shrink=.75)
            fig.suptitle(f"{setting}: conditional diagnostics — {label}\nNeff uses signed sum squared / Sumw2; no goodness-of-fit calibration")
            pdf.savefig(fig); plt.close(fig)
        fig,axes=plt.subplots(1,3,figsize=(13,4.8),constrained_layout=True)
        for ax,(a,b,domain) in zip(axes,(("Local_t","Global","full MM"),("Local_t_phi","Local_t","signal window"),("Local_t_phi","Global","signal window"))):
            values=np.full((3,9),np.nan)
            for cell in document["cells"]:
                aa=cell["applications"].get(a,{}); bb=cell["applications"].get(b,{})
                if aa.get("status")==bb.get("status")=="available":
                    difference=[x-y for x,y in zip(aa["cuts"]["nommcuts"]["values"],bb["cuts"]["nommcuts"]["values"])]
                    value=sum(difference[1:-1]) if domain=="full MM" else window_sum(difference,cell["inputs"]["pion"]["edges"],*window)
                    values[cell["t_index"],cell["phi_index"]]=value
            im=ax.imshow(np.ma.masked_invalid(values),aspect="auto",origin="lower",cmap="coolwarm")
            ax.set(title=f"{a} - {b}\n{domain} signed background",xlabel="phi index",ylabel="t index")
            ax.set_xticks(range(9),range(1,10)); ax.set_yticks(range(3),range(1,4)); fig.colorbar(im,ax=ax,shrink=.75)
        fig.suptitle(f"{setting}: same-child transferred-background differences\n{label}")
        pdf.savefig(fig); plt.close(fig)
        fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
        for ax,key in zip(axes.flat,("B_n","B_delta","B_sidis","A_n","A_delta","A_sidis")):
            for ti in range(3):
                cells=[c for c in document["cells"] if c["t_index"]==ti and (c.get("fit") or {}).get(key) is not None]
                ax.plot([c["phi_index"]+1 for c in cells],[(c["fit"])[key] for c in cells],"o-",label=f"t{ti+1}")
            ax.set(title=key,xlabel="phi index",ylabel="fitted amplitude"); ax.legend()
        fig.suptitle(f"{setting}: independent child amplitudes; {label}\nNo parent prior measurement, physical truth or calibrated errors implied")
        pdf.savefig(fig); plt.close(fig)
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    render(json.loads(args.input.read_text()),args.output)


if __name__ == "__main__":
    main()
