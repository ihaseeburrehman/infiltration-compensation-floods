"""LaTeX table of the Ernz 2018 runs: rainfall, water balance and Larochette scores.
usage: python make_ernz_table.py <results dir> <out tex>"""
import sys, re, pandas as pd, numpy as np
def m(x, fmt):
    """format a number for LaTeX with a proper minus sign and no negative zero"""
    v = float(x); v = 0.0 if abs(v) < 5e-3 else v; v = int(v) if fmt.endswith("d") else v
    return format(v, fmt).replace("-", "$-$")
R, OUT = sys.argv[1:3]
sc = pd.read_csv(f"{R}/scores_ernz.csv"); rt = pd.read_csv(f"{R}/rain_totals.csv").set_index("case"); wb = pd.read_csv(f"{R}/water_balance.csv")
ORDER = ["R_1km_5min", "R_1km_1h_nocorr", "R_eurad_1h_corr", "R_eurad_1h", "F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12", "R_1km_5min_dry", "R_1km_5min_ksx2", "R_1km_5min_ksx10", "R_1km_5min_ksdeep", "R_1km_5min_ksdeepwet", "F_wrf12noda_ksdeep", "F_aifs_ksdeep", "F_graphcast_ksdeep", "F_fuxi_ksdeep", "F_wrf12_ksdeep"]
LAB = {"R_1km_5min": "Radar, RADKLIM corrected", "R_1km_1h_nocorr": "Radar, RADKLIM uncorrected", "R_eurad_1h_corr": "Radar, EURADCLIM corrected", "R_eurad_1h": "Radar, EURADCLIM as issued",
       "F_graphcast": "GraphCast", "F_aifs": "AIFS", "F_fuxi": "FuXi", "F_wrf12noda": "WRF before DA", "F_wrf12": "WRF after DA", "R_1km_5min_dry": "Radar, RADKLIM corrected, dry start",
       "R_1km_5min_ksx2": "Radar, $K_s \\times 2$", "R_1km_5min_ksx10": "Radar, $K_s \\times 10$",
       "R_1km_5min_ksdeep": "Radar, $K_s$ of 30\\,cm layer", "R_1km_5min_ksdeepwet": "Radar, $K_s$ of 30\\,cm layer, wet units pervious",
       "F_wrf12noda_ksdeep": "WRF before DA, $K_s$ of 30\\,cm layer", "F_aifs_ksdeep": "AIFS, $K_s$ of 30\\,cm layer", "F_graphcast_ksdeep": "GraphCast, $K_s$ of 30\\,cm layer", "F_fuxi_ksdeep": "FuXi, $K_s$ of 30\\,cm layer", "F_wrf12_ksdeep": "WRF after DA, $K_s$ of 30\\,cm layer"}
L = ["\\begin{table*}[t]", "\\caption{Ernz 2018 runs: 48\\,h rainfall over the domain and the Larochette catchment, infiltrated fraction of the domain rainfall, and scores at Larochette (observed peak 11.7\\,m$^3$\\,s$^{-1}$ hourly mean). $\\beta$ is the ratio of simulated to observed flow volume. Peak timing error in hours, simulated minus observed (hourly observed peak at 03:00\\,UTC on 1 June). The dry start runs begin from a dry domain instead of the spin up state; the $K_s$ runs (RADKLIM corrected radar and the five forecasts) vary the soil conductivity of the Green--Ampt setting (Sect.~\\ref{sec:ernz}).}",
     "\\label{tab:ernz}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{llrrrrrrrr}", "\\toprule",
     "Forcing & Infiltration & Rain domain & Rain Laroch. & Infiltr. & Peak & Peak ratio & $\\Delta t$ peak & $\\beta$ & NSE \\\\",
     " & & (mm) & (mm) & (\\%) & (m$^3$\\,s$^{-1}$) & & (h) & & \\\\", "\\midrule"]
for c in ORDER:
    for inf, lab in [("noinf", "none"), ("ga", "Green--Ampt")]:
        r = sc[(sc.case == c) & (sc.infiltration == inf)]
        if r.empty: continue
        r = r.iloc[0]; w = wb[(wb.case == c) & (wb.infiltration == inf)]; fi = 100 * w.infiltrated_frac.iloc[0] if len(w) else np.nan
        key = re.sub(r"_(dry|ks\w+)$", "", c)
        L.append(f"{LAB[c]} & {lab} & {rt.loc[key, 'domain']:.1f} & {rt.loc[key, 'larochette']:.1f} & {m(fi, '.0f')} & {r.peak:.1f} & {r.peak_ratio:.2f} & {m(int(r.dt_peak), '+d')} & {r.beta:.2f} & {m(r.nse, '.2f')} \\\\")
L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
open(OUT, "w").write("\n".join(L) + "\n"); print("wrote", OUT)
