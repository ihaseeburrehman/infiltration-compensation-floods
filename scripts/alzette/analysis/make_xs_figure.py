"""Minimal cross section figure for the appendix: water depth (filled) and unit discharge q.n (line) along each gauge
section at the hour of the simulated peak, radar forcing, LISFLOOD-FP, no infiltration, with the mapped channel width and
the station section marked. usage: python make_xs_figure.py <paper_v4 dir>"""
import sys, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P = sys.argv[1].rstrip("/")
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False})
TA = json.load(open(P + "/scripts/extraction/transects_fullvalley.json")); TE = json.load(open(P + "/../ernz2018/scripts/transects_ernz.json"))
ST = {"Livange": (-200, 200), "Hesperange": (-200, 200), "Pfaffenthal": (-200, 200), "Steinsel": (-400, 300), "Larochette": (-300, 300), "Mullerthal": (-200, 200)}   # station section bounds (m), as in score_all.py and xs_ernz.py
rows = [("Livange", P + "/data/segflux/L84_R_1km_5min_noinf.npz", TA), ("Hesperange", P + "/data/segflux/L84_R_1km_5min_noinf.npz", TA),
        ("Pfaffenthal", P + "/data/segflux/L84_R_1km_5min_noinf.npz", TA), ("Steinsel", P + "/data/segflux/L84_R_1km_5min_noinf.npz", TA),
        ("Larochette", P + "/../ernz2018/results/segflux/R_1km_5min_noinf.npz", TE), ("Mullerthal", P + "/../ernz2018/results/segflux/R_1km_5min_noinf.npz", TE)]
fig, axs = plt.subplots(2, 3, figsize=(7.2, 3.9)); axs = axs.ravel()
for ax, (g, f, T) in zip(axs, rows):
    z = np.load(f); t, L, qn, wd = z[g + "_t"], z[g + "_L"], z[g + "_qn"], z[g + "_wd"]; Q = (qn * L).sum(1); i = int(Q.argmax())
    a2 = ax.twinx(); a2.spines["top"].set_visible(False)
    ax.fill_between(t, 0, wd[i], color="#9ecae1", lw=0, label="depth"); a2.plot(t, qn[i], color="#000000", lw=0.9, label="q·n")
    w0 = T[g]["W0"]; sa, sb = ST[g]
    for x in (-w0 / 2, w0 / 2): ax.axvline(x, color="#3f8fd2", lw=0.9, ls=":")
    for x in (sa, sb): ax.axvline(x, color="#000000", lw=0.8, ls="--")
    lo, hi = min(q[2] for q in T[g]["segs"]), max(q[3] for q in T[g]["segs"]); ax.set_xlim(lo, hi)
    ax.set_ylim(0, max(1.05 * wd[i].max(), 0.5)); a2.set_ylim(0, 1.15 * qn[i].max())
    frac = 100 * np.abs(qn[i] * L)[(t >= sa) & (t <= sb)].sum() / np.abs(qn[i] * L).sum(); fs = 100 * np.abs(qn[i] * L)[np.abs(t) <= w0 / 2].sum() / np.abs(qn[i] * L).sum()
    ax.set_title(f"({'abcdef'[list(axs).index(ax)]}) {g} (peak {Q[i]:.0f} m$^3$ s$^{{-1}}$)", loc="left"); print(f"{g:12s} mapped {w0:.0f} m carries {fs:.0f} %, station section {sa:.0f} to {sb:.0f} m carries {frac:.0f} % of the peak {Q[i]:.0f} m3/s at hour {z['steps'][i]}")
    ax.set_ylabel("Depth (m)"); a2.set_ylabel("q·n (m$^2$ s$^{-1}$)")
    if g in ("Steinsel", "Larochette", "Mullerthal"): ax.set_xlabel("Distance from the gauge (m)")
    ax.tick_params(labelsize=7); a2.tick_params(labelsize=7)
from matplotlib.lines import Line2D; from matplotlib.patches import Patch
hd = [Patch(color="#9ecae1", label="Water depth at the simulated peak"), Line2D([], [], color="#000000", lw=0.9, label="Unit discharge normal to the section"),
      Line2D([], [], color="#3f8fd2", lw=0.9, ls=":", label="Mapped channel width"), Line2D([], [], color="#000000", lw=0.8, ls="--", label="Station section")]
fig.legend(handles=hd, loc="lower center", ncol=4, frameon=False, fontsize=7); fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(P + "/figures/figB1_cross_sections.pdf"); print("saved", P + "/figures/figB1_cross_sections.pdf")
