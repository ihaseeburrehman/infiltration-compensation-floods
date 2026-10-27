"""Performance diagram: peak ratio (log axis) against NSE, one point per forcing and configuration, with an arrow from the
run without infiltration to the run with Green-Ampt. Left: Alzette gauge means (station section), circles LISFLOOD-FP,
squares HiPIMS, rescaled forecasts as smaller markers with dashed arrows. Right: Ernz, Larochette. NSE below the axis limit is
drawn at the limit. usage: python make_perf_diagram.py <paper_v4>"""
import sys, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P = sys.argv[1].rstrip("/")
import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import figstyle as FS
C = ["R_1km_5min", "F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12", "VM_graphcast", "VM_aifs", "VM_fuxi", "VM_wrf12noda", "VM_wrf12"]
COL = {c: FS.col(c) for c in C}
LAB = {"R_1km_5min": "Radar", "F_graphcast": "GraphCast", "F_aifs": "AIFS", "F_fuxi": "FuXi", "F_wrf12noda": "WRF before DA", "F_wrf12": "WRF after DA"}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1.25, 1]))
def draw(ax, pts, title, XMIN):
    """pts: list of (case, marker, x_none, y_none, x_ga, y_ga)"""
    for c, mk, x0, y0, x1, y1 in pts:
        vm = c.startswith("VM_") or c.endswith("_ksdeep"); ms = 4.5 if vm else 6.5; ls = "--" if vm else "-"; c = c.replace("_ksdeep", "")
        X0, X1 = max(x0, XMIN), max(x1, XMIN)
        ax.annotate("", xy=(X1, y1), xytext=(X0, y0), arrowprops=dict(arrowstyle="-|>", color=COL[c], lw=0.8, ls=ls, shrinkA=3, shrinkB=3))
        ax.plot(X0, y0, mk, mfc="white", mec=COL[c], ms=ms, mew=1.0, zorder=3); ax.plot(X1, y1, mk, mfc=COL[c], mec=COL[c], ms=ms, zorder=3)
        for X, x in ((X0, x0), (X1, x1)):
            if x < XMIN: ax.plot(X - 0.04, (y0 if X == X0 else y1), "<", color=COL[c], ms=3, zorder=4)
    ax.axhline(1, color="grey", lw=0.6); ax.axvline(0, color="grey", lw=0.6, ls=":")
    ax.set_yscale("log"); ax.set_xlabel("NSE"); ax.set_ylabel("Simulated / observed peak"); ax.set_title(title, loc="left", fontsize=8)
    ax.set_xlim(XMIN - 0.1, 1.0)
    # Alzette
s = pd.read_csv(P + "/data/scores.csv"); s = s[(s.width == "full") & s.case.isin(C) & s.infiltration.isin(["none", "GA"])]
gm = s.groupby(["case", "model", "infiltration"])[["peak_ratio", "nse"]].mean()
pts = []
for c in C:
    for model, mk in (("LISFLOOD-FP", "o"), ("HiPIMS", "s")):
        a, b = gm.loc[(c, model, "none")], gm.loc[(c, model, "GA")]; pts.append((c, mk, a.nse, a.peak_ratio, b.nse, b.peak_ratio))
draw(axs[0], pts, "(a) Alzette 2021, mean of four gauges", -0.6); axs[0].set_ylim(0.2, 3)
# Ernz
e = pd.read_csv(P + "/../ernz2018/results/scores_ernz.csv").set_index(["case", "infiltration"]); pts = []
for c in C[:6]:
    if (c, "noinf") not in e.index: continue
    a, b = e.loc[(c, "noinf")], e.loc[(c, "ga")]; pts.append((c, "o", a.nse, a.peak_ratio, b.nse, b.peak_ratio))
    if (c + "_ksdeep", "ga") in e.index: d = e.loc[(c + "_ksdeep", "ga")]; pts.append((c + "_ksdeep", "o", b.nse, b.peak_ratio, d.nse, d.peak_ratio))
draw(axs[1], pts, "(b) Ernz 2018, Larochette", -2.0); axs[1].set_ylim(0.03, 12)
hd = [Line2D([], [], color=COL[c], marker="o", lw=0, ms=5, label=LAB[c]) for c in C[:6]]
hd += [Line2D([], [], color="grey", marker="o", mfc="white", lw=0, ms=5, label="no infiltration"), Line2D([], [], color="grey", marker="o", lw=0, ms=5, label="Green–Ampt"),
       Line2D([], [], color="grey", marker="s", mfc="white", lw=0, ms=5, label="HiPIMS (squares)"), Line2D([], [], color="grey", ls="--", lw=0.8, label="rescaled forecast (a), $K_s$ of 30 cm layer (b) (small)")]
fig.legend(handles=hd, loc="lower center", ncol=5, frameon=False, fontsize=7); fig.tight_layout(rect=(0, 0.12, 1, 1))
fig.savefig(P + "/figures/fig11_perf_diagram.pdf"); fig.savefig(P + "/figures/fig11_perf_diagram.png", dpi=300); print("saved"); print(gm.round(2).to_string())
