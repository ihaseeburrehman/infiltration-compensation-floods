"""Figures for the paper from data/*.csv. Every forcing has one fixed colour in every figure.
usage: python make_figures.py <paper_v4 dir>"""
import sys, os, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import figstyle as FS
P = sys.argv[1].rstrip("/"); D, F = P + "/data/", P + "/figures/"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e3e3e3", "grid.linewidth": 0.5, "lines.linewidth": 1.1})
G = ["Livange", "Hesperange", "Pfaffenthal", "Steinsel"]
# one colour per forcing, used in every figure: radar black, then Paul Tol's colour blind safe "muted" set
C = ["R_1km_5min", "F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12", "VM_graphcast", "VM_aifs", "VM_fuxi", "VM_wrf12noda", "VM_wrf12"]
COL = {c: FS.col(c) for c in C}
LS = {c: ("--" if c.startswith("VM_") else "-") for c in C}   # rescaled forecasts dashed, in the colour of their forecast
LAB = {"R_1km_5min": "Radar", "F_graphcast": "GraphCast", "F_aifs": "AIFS", "F_fuxi": "FuXi", "F_wrf12noda": "WRF before DA", "F_wrf12": "WRF after DA",
       "VM_fuxi": "FuXi VM", "VM_graphcast": "GraphCast VM", "VM_aifs": "AIFS VM", "VM_wrf12noda": "WRF before DA VM", "VM_wrf12": "WRF after DA VM"}
LW = {c: (1.6 if c == "R_1km_5min" else 1.0) for c in C}
INF = {"none": "no infiltration", "GA": "Green–Ampt infiltration"}
sc = pd.read_csv(D + "scores.csv"); hy = pd.read_csv(D + "hydrographs.csv"); ob = pd.read_csv(D + "obs_utc.csv")
wb = pd.read_csv(D + "water_balance.csv")
def days(a, last=True):
    a.set_xlim(0, 96); a.set_xticks(range(6, 97, 24))
    a.set_xticklabels(["14", "15", "16", "17"] if last else [])
    if last: a.set_xlabel("July 2021 (00:00 UTC)")
def legend_forcings(fig, extra=(), ncol=6, y=0.0):
    hd = [Line2D([], [], color=COL[c], ls=LS[c], lw=LW[c] + 0.6, label=LAB[c]) for c in C] + list(extra)
    fig.legend(handles=hd, loc="lower center", ncol=ncol, frameon=False, fontsize=7, bbox_to_anchor=(0.5, y))
OBS = Line2D([], [], color=FS.OBS, lw=0, marker="o", ms=3, label="Observed")

# Fig. 5: hydrographs, all forcings, LISFLOOD-FP, full wetted section; rows = infiltration setting
fig, ax = plt.subplots(2, 4, figsize=(7.2, 4.6), sharex=True)
for j, g in enumerate(G):
    o = ob[ob.gauge == g]; ymax = 0
    for i, inf in enumerate(("none", "GA")):
        a = ax[i, j]; a.plot(o.h, o.Q, color=FS.OBS, lw=0, marker="o", ms=1.8, zorder=5)
        for c in C:
            s = hy[(hy.case == c) & (hy.model == "LISFLOOD-FP") & (hy.infiltration == inf) & (hy.gauge == g) & (hy.width == "full")]
            a.plot(s.h, s.Q, color=COL[c], ls=LS[c], lw=LW[c], zorder=4 if c == "R_1km_5min" else 3); ymax = max(ymax, s.Q.max())
        days(a, i == 1); FS.panel(a, 4 * i + j)
        if i == 0: a.set_title(g)
    for i in range(2): ax[i, j].set_ylim(0, ymax * 1.05)
for i, inf in enumerate(("none", "GA")): ax[i, 0].set_ylabel(f"Discharge (m$^3$ s$^{{-1}}$)\n{INF[inf]}")
legend_forcings(fig, (OBS,), ncol=6); fig.tight_layout(rect=(0, 0.10, 1, 1))
fig.savefig(F + "fig04_hydrographs.pdf"); fig.savefig(F + "fig04_hydrographs.png", dpi=300); plt.close(fig)

# Fig. C1: same for HiPIMS
fig, ax = plt.subplots(2, 4, figsize=(7.2, 4.6), sharex=True)
for j, g in enumerate(G):
    o = ob[ob.gauge == g]; ymax = 0
    for i, inf in enumerate(("none", "GA")):
        a = ax[i, j]; a.plot(o.h, o.Q, color=FS.OBS, lw=0, marker="o", ms=1.8, zorder=5)
        for c in C:
            s = hy[(hy.case == c) & (hy.model == "HiPIMS") & (hy.infiltration == inf) & (hy.gauge == g) & (hy.width == "full")]
            a.plot(s.h, s.Q, color=COL[c], ls=LS[c], lw=LW[c]); ymax = max(ymax, s.Q.max())
        days(a, i == 1); FS.panel(a, 4 * i + j)
        if i == 0: a.set_title(g)
    for i in range(2): ax[i, j].set_ylim(0, ymax * 1.05)
for i, inf in enumerate(("none", "GA")): ax[i, 0].set_ylabel(f"Discharge (m$^3$ s$^{{-1}}$)\n{INF[inf]}")
legend_forcings(fig, (OBS,), ncol=6); fig.tight_layout(rect=(0, 0.10, 1, 1))
fig.savefig(F + "figC1_hydrographs_hipims.pdf"); plt.close(fig)

# Fig. 7: peak ratio against volume ratio; colour = forcing, marker = model
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.5), sharex=True, sharey=True)
for a, inf in zip(ax, ("none", "GA")):
    a.axhline(1, color="#999999", lw=0.8); a.axvline(1, color="#999999", lw=0.8)
    for c in C:
        s = sc[(sc.case == c) & (sc.infiltration == inf) & (sc.width == "full")]
        for m, mk in (("LISFLOOD-FP", "o"), ("HiPIMS", "s")):
            u = s[s.model == m]
            a.scatter(u.beta, u.peak_ratio, s=18, marker=mk, facecolor="white" if c.startswith("VM_") else COL[c], edgecolor=COL[c] if c.startswith("VM_") else "white", linewidth=0.9 if c.startswith("VM_") else 0.5, zorder=3)
    a.set_title(f"({'a' if inf == 'none' else 'b'}) {INF[inf][0].upper() + INF[inf][1:]}", loc="left")
    a.set_xlabel("Volume ratio $\\beta$")
ax[0].set_ylabel("Simulated / observed peak")
mk = [Line2D([], [], marker="o", ls="", color="#555555", label="LISFLOOD-FP"), Line2D([], [], marker="s", ls="", color="#555555", label="HiPIMS")]
hd = [Line2D([], [], marker="o", ls="", color=COL[c], mfc="white" if c.startswith("VM_") else COL[c], label=LAB[c]) for c in C] + mk
fig.legend(handles=hd, loc="lower center", ncol=5, frameon=False, fontsize=7); fig.tight_layout(rect=(0, 0.17, 1, 1))
fig.savefig(F + "fig06_peak_volume.pdf"); plt.close(fig)

# Fig. 3: water balance with infiltration (LISFLOOD-FP); forcing colour on the x labels' marker strip
w = wb[wb.infiltration == "GA"].set_index("case").loc[C]
fig, a = plt.subplots(figsize=(7.2, 2.9))
x = np.arange(len(C)); b0 = np.zeros(len(C))
for k, col, lab in (("infiltrated", "#525252", "Infiltrated"), ("outflow", "#9a9a9a", "Left the domain"), ("dstore", "#d9d9d9", "Added to surface storage")):
    v = (w.storage1 - w.storage0) if k == "dstore" else w[k]
    a.bar(x, v, bottom=b0, color=col, width=0.7, edgecolor="white", linewidth=1, label=lab); b0 += v.values
for i, (c, r, f) in enumerate(zip(C, w.rain, w.infiltrated / w.rain)):
    a.text(i, r + 0.6, f"{100 * f:.0f} %", ha="center", fontsize=7)
    a.add_patch(plt.Rectangle((i - 0.35, -2.2), 0.7, 1.0, color=COL[c], clip_on=False))
a.set_xticks(x); a.set_xticklabels([LAB[c].replace("WRF ", "WRF\n").replace(" VM", "\nVM") for c in C])
a.tick_params(axis="x", pad=10); a.set_ylabel("Volume (million m$^3$)"); a.set_ylim(0, 45)
a.legend(frameon=False, fontsize=7, ncol=3, loc="upper right"); fig.tight_layout(); fig.savefig(F + "fig02_water_balance.pdf"); plt.close(fig)

# Fig. 8: Walferdange water level, both models
if os.path.exists(D + "stage_walferdange.csv"):
    st = pd.read_csv(D + "stage_walferdange.csv"); so = pd.read_csv(D + "stage_obs_walferdange.csv")
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 4.4), sharex=True, sharey=True)
    for i, inf in enumerate(("none", "GA")):
        for j, m in enumerate(("LISFLOOD-FP", "HiPIMS")):
            a = ax[i, j]; a.plot(so.h, so.stage, color=FS.OBS, lw=0, marker="o", ms=1.8, zorder=5)
            for c in C:
                s = st[(st.case == c) & (st.model == m) & (st.infiltration == inf)]
                a.plot(s.h, s.stage, color=COL[c], ls=LS[c], lw=LW[c], zorder=4 if c == "R_1km_5min" else 3)
            days(a, i == 1); FS.panel(a, 2 * i + j)
            if i == 0: a.set_title(m)
        ax[i, 0].set_ylabel(f"Stage (m)\n{INF[inf]}")
    legend_forcings(fig, (OBS,), ncol=6); fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.savefig(F + "fig07_walferdange.pdf"); plt.close(fig)
print("figures written to", F)
