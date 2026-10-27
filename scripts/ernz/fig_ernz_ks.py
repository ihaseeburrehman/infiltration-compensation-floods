"""Ernz 2018: Larochette and Mullerthal hydrographs of the radar run for the soil conductivity variants (tab:ernz).
usage: python fig_ernz_ks.py <results dir> <out pdf>"""
import sys, pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R, OUT = sys.argv[1:3]
import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import figstyle as FS
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False, "grid.color": "#e3e3e3", "grid.linewidth": 0.5})
C = [("noinf", "No infiltration", "#c6c6c6", "-"), ("ga", "Column $K_s$ (0.34 mm h$^{-1}$)", "#969696", "-"),
     ("ksx2_ga", "$K_s \\times 2$", "#6e6e6e", "-"), ("ksx10_ga", "$K_s \\times 10$", "#4a4a4a", "-"),
     ("ksdeep_ga", "$K_s$ of 30 cm layer (22 mm h$^{-1}$)", "#000000", "-"), ("ksdeepwet_ga", "30 cm layer, wet units pervious", "#000000", "--")]
obs = pd.read_csv(f"{R}/obs_larochette_utc.csv")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.6))
for k, lab, col, ls in C:
    h = pd.read_csv(f"{R}/hydro_R_1km_5min_{k}.csv")
    ax[0].plot(h.hour, h.Larochette, color=col, ls=ls, lw=1.2, label=lab); ax[1].plot(h.hour, h.Mullerthal, color=col, ls=ls, lw=1.2)
ax[0].plot(obs.h, obs.Q, "o", color=FS.OBS, ms=2.5, label="Observed", zorder=5)
ax[0].set_title("(a) Larochette, Ernz Blanche", loc="left"); ax[1].set_title("(b) Mullerthal, Ernz Noire", loc="left")
for a in ax: a.set_xlim(0, 36); a.set_xlabel("Hours from 31 May 2018 18:00 UTC"); a.grid(True)
ax[0].set_ylabel("Discharge (m$^3$ s$^{-1}$)"); ax[1].set_ylabel("Discharge (m$^3$ s$^{-1}$)")
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, frameon=False)
fig.tight_layout(rect=(0, 0.14, 1, 1)); fig.savefig(OUT); print("wrote", OUT)
