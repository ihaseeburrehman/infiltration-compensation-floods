"""Score figures that replace the result tables: scorecard (gauge means), gauge heatmap, sensitivity runs, surveyed
section and the radar against rain gauges. Same forcing colours as make_figures.py.
usage: python make_scorecards.py <paper_v4 dir>"""
import sys, re, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import figstyle as FS
P = sys.argv[1].rstrip("/"); D, F = P + "/data/", P + "/figures/"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
                     "axes.grid": True, "grid.color": "#e3e3e3", "grid.linewidth": 0.5})
C = ["R_1km_5min", "F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12", "VM_graphcast", "VM_aifs", "VM_fuxi", "VM_wrf12noda", "VM_wrf12"]
COL = {c: FS.col(c) for c in C}
COL["R_1km_6h"] = "#000000"
LAB = {"R_1km_5min": "Radar", "F_graphcast": "GraphCast", "F_aifs": "AIFS", "F_fuxi": "FuXi", "F_wrf12noda": "WRF before DA", "F_wrf12": "WRF after DA",
       "VM_graphcast": "GraphCast VM", "VM_aifs": "AIFS VM", "VM_fuxi": "FuXi VM", "VM_wrf12noda": "WRF before DA VM", "VM_wrf12": "WRF after DA VM", "R_1km_6h": "Radar 6 h"}
G = ["Livange", "Hesperange", "Pfaffenthal", "Steinsel"]
MK = {"LISFLOOD-FP": "o", "HiPIMS": "s"}; DY = {"LISFLOOD-FP": -0.13, "HiPIMS": 0.13}
INK2 = "#555555"
sc = pd.read_csv(D + "scores.csv"); st = pd.read_csv(D + "stage_scores.csv").rename(columns={"inf": "infiltration"})
wb = pd.read_csv(D + "water_balance.csv")
try: sar = pd.read_csv(D + "sar_scores.csv"); sar = sar[(sar.time == "sar15") & (sar.zone == "riv500") & (sar.thr == 0.1)]
except FileNotFoundError: sar = None
mean = lambda w: sc[sc.width == w].groupby(["case", "model", "infiltration"])[["peak_ratio", "nse", "beta"]].mean()

def dots(ax, rows, val, ref=None):
    """rows: list of (y, case, infiltration, filled); val(case, model, inf) -> value or nan."""
    for y, c, inf, filled in rows:
        for m in MK:
            v = val(c, m, inf)
            if np.isnan(v): continue
            k = FS.inf_col("none" if not filled else "GA"); ax.plot(v, y + DY[m], MK[m], ms=4.4, mec="white", mew=0.4, mfc=k, zorder=3 if filled else 2)
    if ref is not None: ax.axvline(ref, color=INK2, lw=0.8, zorder=1)
    ax.grid(axis="y", visible=False); ax.tick_params(axis="y", length=0)

def labels(ax, ylab):
    ax.set_yticks(list(ylab)); ax.set_yticklabels(list(ylab.values())); ax.invert_yaxis()

KEY = [Line2D([], [], marker="o", ls="", mfc=FS.INF["none"], mec="white", ms=6, label="No infiltration"),
       Line2D([], [], marker="o", ls="", mfc=FS.INF["GA"], mec="white", ms=6, label="Green–Ampt"),
       Line2D([], [], marker="o", ls="", mfc="white", mec="#555555", label="LISFLOOD-FP"),
       Line2D([], [], marker="s", ls="", mfc="white", mec="#555555", label="HiPIMS")]

def scorecard(width, fname, with_stage):
    s = mean(width)
    g = lambda col: (lambda c, m, i: s.loc[(c, m, i), col] if (c, m, i) in s.index else np.nan)
    stg = st.set_index(["case", "model", "infiltration"])
    panels = [("Peak ratio", g("peak_ratio"), 1.0), ("NSE", g("nse"), None)]
    if with_stage:
        panels.append(("Peak stage error (m)", lambda c, m, i: stg.loc[(c, m, i), "dpeak"] if (c, m, i) in stg.index else np.nan, 0.0))
        if sar is not None:
            q = sar.set_index(["case", "model", "infiltration"])
            panels.append(("Flood extent CSI", lambda c, m, i: q.loc[(c, m, i), "CSI"] if (c, m, i) in q.index else np.nan, None))
    fig, axs = plt.subplots(1, len(panels), figsize=(7.2 if len(panels) == 4 else 4.6, 3.8), sharey=True)
    ylab = {}; rows = []
    for k, c in enumerate(C):
        y = k + (0.6 if k >= 6 else 0) + (0.3 if k >= 1 else 0); ylab[y] = LAB[c]
        rows += [(y, c, "none", False), (y, c, "GA", True)]
    for k, (ax, (t, f, ref)) in enumerate(zip(axs, panels)):
        dots(ax, rows, f, ref); ax.set_title(f"({'abcd'[k]}) {t}", fontsize=8, loc="left")
    axs[1].axvline(0, color=INK2, lw=0.5, ls=":", zorder=1)
    labels(axs[0], ylab)
    fig.legend(handles=KEY, loc="lower center", ncol=4, frameon=False, fontsize=7)
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(F + fname); plt.close(fig)

scorecard("full", "fig03_scorecard.pdf", True)
scorecard("surveyed", "figD1_scorecard_mapped.pdf", False)

# Gauge heatmap: NSE and peak ratio per gauge, radar and forecasts as issued, mean of the two models (LF/HP in the cell)
s = sc[sc.width == "full"].set_index(["case", "model", "infiltration", "gauge"])
CF = C[:6]
fig, axs = plt.subplots(2, 2, figsize=(6.6, 4.5), sharex=True, sharey=True)
for j, inf in enumerate(("none", "GA")):
    for i, (col, cmap, lim) in enumerate((("peak_ratio", "PuOr_r", (0, 2)), ("nse", "Greys", (0, 1)))):
        ax = axs[i, j]; a = np.array([[s.loc[(c, "LISFLOOD-FP", inf, g), col] for g in G] for c in CF])
        b = np.array([[s.loc[(c, "HiPIMS", inf, g), col] for g in G] for c in CF]); m = (a + b) / 2
        cm = plt.get_cmap(cmap)   # cells as vector rectangles, so no raster interpolation in any viewer
        for r in range(len(CF)):
            for k in range(len(G)): ax.add_patch(plt.Rectangle((k - 0.5, r - 0.5), 1, 1, color=cm((np.clip(m[r, k], *lim) - lim[0]) / (lim[1] - lim[0])), lw=0))
        ax.set_xlim(-0.5, len(G) - 0.5); ax.set_ylim(len(CF) - 0.5, -0.5)
        for r in range(len(CF)):
            for k in range(len(G)):
                dark = (cmap == "Greys" and m[r, k] > 0.55) or (cmap == "PuOr_r" and abs(m[r, k] - 1) > 0.65)
                ax.text(k, r, f"{a[r, k]:.2f}\n{b[r, k]:.2f}".replace("-", "−"), ha="center", va="center", fontsize=7,
                        color="white" if dark else "#222222", linespacing=0.95)
        ax.set_xticks(range(4)); ax.set_xticklabels(G); ax.set_yticks(range(len(CF))); ax.set_yticklabels([LAB[c] for c in CF])
        ax.grid(False); ax.tick_params(length=0)
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.set_title(f"({'abcd'[2 * i + j]}) " + ("Peak ratio" if col == "peak_ratio" else "NSE") + (", no infiltration" if inf == "none" else ", Green–Ampt"),
                     fontsize=8, loc="left")
fig.text(0.5, 0.005, "Each cell: LISFLOOD-FP (top) and HiPIMS (bottom); colour is their mean.", ha="center", fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(F + "fig05_gauges.pdf"); plt.close(fig)

# Sensitivity runs: conductivity, 6 h radar and wet start
s = mean("full"); stg = st.set_index(["case", "model", "infiltration"])
wbi = wb.set_index(["case", "infiltration"])
SR = [("R_1km_5min", "GA_ks05", "Radar, ½ Ks"), ("R_1km_5min", "GA", "Radar, Ks"), ("R_1km_5min", "GA_ks2", "Radar, 2 Ks"),
      ("R_1km_5min", "GA_sm", "Radar, wet start"), ("R_1km_6h", "none", "Radar 6 h, none"), ("R_1km_6h", "GA", "Radar 6 h, Ks"),
      ("F_graphcast", "GA_ks05", "GraphCast, ½ Ks"), ("F_graphcast", "GA", "GraphCast, Ks"), ("F_graphcast", "GA_ks2", "GraphCast, 2 Ks"),
      ("F_graphcast", "GA_sm", "GraphCast, wet start"), ("F_aifs", "GA", "AIFS, Ks"), ("F_aifs", "GA_sm", "AIFS, wet start"),
      ("F_fuxi", "GA", "FuXi, Ks"), ("F_fuxi", "GA_sm", "FuXi, wet start")]
ylab, rows = {}, []
for k, (c, inf, lab) in enumerate(SR):
    y = k + 0.5 * sum(k >= b for b in (4, 6, 10, 12)); ylab[y] = lab; rows.append((y, c, inf, inf != "none"))
val = lambda col: (lambda c, m, i: s.loc[(c, m, i), col] if (c, m, i) in s.index else np.nan)
panels = [("Peak ratio", val("peak_ratio"), 1.0), ("NSE", val("nse"), None),
          ("Peak stage error (m)", lambda c, m, i: stg.loc[(c, m, i), "dpeak"] if (c, m, i) in stg.index else np.nan, 0.0),
          ("Infiltrated (%)", lambda c, m, i: 100 * wbi.loc[(c, i), "infiltrated"] / wbi.loc[(c, i), "rain"]
           if m == "LISFLOOD-FP" and (c, i) in wbi.index else np.nan, None)]
fig, axs = plt.subplots(1, 4, figsize=(7.2, 3.9), sharey=True)
for k, (ax, (t, f, ref)) in enumerate(zip(axs, panels)):
    dots(ax, rows, f, ref); ax.set_title(f"({'abcd'[k]}) {t}", fontsize=8, loc="left")
labels(axs[0], ylab)
fig.legend(handles=KEY, loc="lower center", ncol=4, frameon=False, fontsize=7)
fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(F + "fig08_sensitivity.pdf"); plt.close(fig)

# Radar against 18 independent rain gauges (appendix), parsed from the gauge list kept in data/
rg = pd.read_csv(D + "raingauges.csv")
fig, ax = plt.subplots(figsize=(3.2, 3.0))
ax.plot([50, 130], [50, 130], color=INK2, lw=0.8); ax.scatter(rg.gauge, rg.radar, s=14, color="#000000", zorder=3)
for _, r in rg.iterrows():
    if abs(r.radar / r.gauge - 1) > 0.15: ax.annotate(r["name"], (r.gauge, r.radar), fontsize=7, xytext=(3, 2), textcoords="offset points")
ax.set_xlim(55, 128); ax.set_ylim(55, 128); ax.set_aspect("equal")
ax.set_xlabel("Rain gauge total (mm)"); ax.set_ylabel("Radar total at gauge (mm)"); ax.spines["left"].set_visible(True)
fig.tight_layout(); fig.savefig(F + "figA1_raingauges.pdf"); plt.close(fig)
print("written", F)
