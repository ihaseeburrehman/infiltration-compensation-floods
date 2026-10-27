"""Hydrographs of the Ernz 2018 runs: Larochette (observed, CAMELS-LUX) and Mullerthal (no gauge), without and with infiltration.
Colours as in paper_v4 make_figures.py. usage: python fig_ernz_hydro.py <results dir> <out pdf> [png]"""
import sys, os, glob, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R, OUT = sys.argv[1], sys.argv[2]
import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import figstyle as FS
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "font.family": "sans-serif", "axes.linewidth": 0.6, "grid.color": "#e3e3e3", "grid.linewidth": 0.5, "lines.linewidth": 1.1})
C = ["R_1km_5min", "R_1km_1h_nocorr", "R_eurad_1h_corr", "F_graphcast", "F_aifs", "F_fuxi", "F_wrf12", "F_wrf12noda"]
COL = {c: FS.col(c) for c in C}; COL["R_eurad_1h_corr"] = "#6b6b6b"
LAB = {"R_1km_5min": "Radar (RADKLIM, gauge corrected)", "R_1km_1h_nocorr": "Radar (RADKLIM, uncorrected)", "R_eurad_1h_corr": "Radar (EURADCLIM, gauge corrected)",
       "F_graphcast": "GraphCast", "F_aifs": "AIFS", "F_fuxi": "FuXi", "F_wrf12": "WRF after DA", "F_wrf12noda": "WRF before DA"}
LS = {"R_1km_1h_nocorr": "--", "R_eurad_1h_corr": "-."}
T0 = pd.Timestamp("2018-05-31T18")
obs = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../data/camels_larochette_15min_2018.csv"), parse_dates=["Date"])
obs["h"] = (obs.Date - pd.Timedelta(hours=1) - T0).dt.total_seconds() / 3600; obs = obs[(obs.h >= 0) & (obs.h <= 48)]
fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.0), sharex=True)
for i, inf in enumerate(["noinf", "ga"]):
    for j, (site, title) in enumerate([("Larochette", "Larochette, Ernz Blanche (69 km$^2$)"), ("Mullerthal", "Mullerthal, Ernz Noire (no gauge)")]):
        a = ax[i, j]
        if site == "Larochette": a.plot(obs.h, obs.Q, color=FS.OBS, lw=0, marker="o", ms=1.8, zorder=5)
        else: pass  # the AGE design flood of the 2024 station is not cited in the paper and no longer drawn
        for c in C:
            f = f"{R}/hydro_{c}_{inf}.csv"
            if not os.path.exists(f): continue
            h = pd.read_csv(f); a.plot(h.hour, h[site], color=COL[c], ls=LS.get(c, "-"), lw=1.6 if c == "R_1km_5min" else 1.0, zorder=4 if c == "R_1km_5min" else 3)
        a.grid(True); a.set_xlim(0, 48); a.set_xticks(range(0, 49, 6))
        if i == 0: a.set_title(f"({'ab'[j]}) {title}", loc="left")
        else: FS.panel(a, 2 + j)
        if j == 0: a.set_ylabel("Discharge (m$^3$ s$^{-1}$)\n" + ("no infiltration" if inf == "noinf" else "Green–Ampt infiltration"))
        if i == 1: a.set_xlabel("Hours from 31 May 2018 18:00 UTC")
hd = [Line2D([], [], color=COL[c], ls=LS.get(c, "-"), lw=1.6, label=LAB[c]) for c in C] + [Line2D([], [], color=FS.OBS, lw=0, marker="o", ms=3, label="Observed (AGE via CAMELS-LUX)")]
fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False)
fig.tight_layout(rect=(0, 0.13, 1, 1)); fig.savefig(OUT)
if len(sys.argv) > 3: fig.savefig(sys.argv[3], dpi=120)
print("saved", OUT)
