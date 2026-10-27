"""LaTeX tables for the paper from data/scores.csv and data/water_balance.csv (no hand transcription).
usage: python make_tables.py <paper_v4 dir>"""
import sys, pandas as pd
P = sys.argv[1].rstrip("/"); D = P + "/data/"; M = P + "/manuscript/"
C = ["R_1km_5min", "F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12", "VM_graphcast", "VM_aifs", "VM_fuxi", "VM_wrf12noda", "VM_wrf12"]
LAB = {"R_1km_5min": "Radar", "F_graphcast": "GraphCast", "F_aifs": "AIFS", "F_fuxi": "FuXi", "F_wrf12noda": "WRF before DA", "F_wrf12": "WRF after DA",
       "VM_fuxi": "FuXi VM", "VM_graphcast": "GraphCast VM", "VM_aifs": "AIFS VM", "VM_wrf12noda": "WRF before DA VM", "VM_wrf12": "WRF after DA VM"}
G = ["Livange", "Hesperange", "Pfaffenthal", "Steinsel"]
sc = pd.read_csv(D + "scores.csv"); wb = pd.read_csv(D + "water_balance.csv")
f = lambda x: f"{x + 0.0:.2f}".replace("-0.00", "0.00").replace("-", "$-$")

def main_table(width, label, cap):
    s = sc[sc.width == width].groupby(["case", "model", "infiltration"])[["peak_ratio", "beta", "nse", "kge", "dt_peak"]].mean()
    L = ["\\begin{table*}[t]", f"\\caption{{{cap}}}", f"\\label{{{label}}}", "\\centering", "\\begin{tabular}{llrrrrrrrrrr}", "\\toprule",
         " & & \\multicolumn{5}{c}{LISFLOOD-FP} & \\multicolumn{5}{c}{HiPIMS} \\\\",
         "\\cmidrule(lr){3-7}\\cmidrule(lr){8-12}",
         "Rainfall & Infiltration & Peak & $\\beta$ & NSE & KGE$'$ & PTE & Peak & $\\beta$ & NSE & KGE$'$ & PTE \\\\", "\\midrule"]
    for i, c in enumerate(C):
        if i in (1, 6): L.append("\\midrule")
        for inf, il in (("none", "none"), ("GA", "Green--Ampt")):
            row = [LAB[c] if inf == "none" else "", il]
            for m in ("LISFLOOD-FP", "HiPIMS"):
                v = s.loc[(c, m, inf)]
                row += [f(v.peak_ratio), f(v.beta), f(v.nse), f(v.kge), f"{v.dt_peak:+.1f}".replace("-", "$-$")]
            L.append(" & ".join(row) + " \\\\")
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
    return "\n".join(L)

def gauge_table():
    s = sc[(sc.width == "full")].set_index(["case", "model", "infiltration", "gauge"])
    L = ["\\begin{table*}[t]", "\\caption{Peak ratio (simulated over observed hourly peak) and NSE at each gauge for the radar and the forecasts as issued, station section. Each cell gives LISFLOOD-FP / HiPIMS.}",
         "\\label{tab:gauges_scores}", "\\centering", "\\begin{tabular}{llllll}", "\\toprule",
         "Rainfall & Infiltration & " + " & ".join(G) + " \\\\", "\\midrule"]
    for c in C[:6]:
        for inf, il in (("none", "none"), ("GA", "Green--Ampt")):
            row = [LAB[c] if inf == "none" else "", il]
            for g in G:
                a, b = s.loc[(c, "LISFLOOD-FP", inf, g)], s.loc[(c, "HiPIMS", inf, g)]
                row.append(f"{a.peak_ratio:.2f}/{b.peak_ratio:.2f} ({a.nse:.2f}/{b.nse:.2f})".replace("-", "$-$"))
            L.append(" & ".join(row) + " \\\\")
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
    return "\n".join(L)

def wb_table():
    L = ["\\begin{table}[t]", "\\caption{Water balance of the LISFLOOD-FP simulations over 96\\,h (million m$^3$). All runs start with 1.86 million m$^3$ of surface water. The residual of every run is below 0.01 million m$^3$.}",
         "\\label{tab:wbal}", "\\centering", "\\begin{tabular}{lrrrrrr}", "\\toprule",
         " & & & \\multicolumn{2}{c}{No infiltration} & \\multicolumn{2}{c}{Green--Ampt} \\\\", "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
         "Rainfall & mm & Rain & Outflow & Stored & Infiltrated & Outflow \\\\", "\\midrule"]
    for c in C:
        n = wb[(wb.case == c) & (wb.infiltration == "none")].iloc[0]; g = wb[(wb.case == c) & (wb.infiltration == "GA")].iloc[0]
        L.append(f"{LAB[c]} & {n.rain / 0.4677:.1f} & {n.rain:.1f} & {n.outflow:.1f} & {n.storage1:.1f} & {g.infiltrated:.1f} ({100 * g.infiltrated / g.rain:.0f}\\,\\%) & {g.outflow:.1f} \\\\")
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    return "\n".join(L)

out = [main_table("full", "tab:main", "Gauge mean scores for the 44 main simulations, station section. Peak is the ratio of simulated to observed peak, $\\beta$ the volume ratio and PTE the peak timing error in hours (negative means early). VM denotes forecasts rescaled to the radar event total."),
       main_table("surveyed", "tab:surveyed", "As Table~\\ref{tab:main}, but with discharge integrated across the surveyed channel width only (sensitivity test)."),
       main_table("valley", "tab:valley", "As Table~\\ref{tab:main}, but with discharge integrated across the whole transect at every gauge, so that at Livange the bypass flow in the ditch south east of the station is included (sensitivity test)."),
       gauge_table(), wb_table()]
open(M + "tables_generated.tex", "w").write("\n\n".join(out) + "\n")
print("written", M + "tables_generated.tex")

# Walferdange water level: peak error (m) and peak timing error (h), from data/stage_scores.csv
st = pd.read_csv(D + "stage_scores.csv").set_index(["case", "model", "inf"])
L = ["\\begin{table}[t]", "\\caption{Peak stage at Walferdange: simulated minus observed peak stage (m) and timing error (h, negative means early). Observed stage is measured above the gauge zero, simulated stage above the lowest terrain cell within 15\\,m of the gauge. The observed hourly peak stage is 3.53\\,m. LF is LISFLOOD-FP, HP is HiPIMS.}",
     "\\label{tab:stage}", "\\centering", "\\begin{tabular}{lrrrrrrrr}", "\\toprule",
     " & \\multicolumn{4}{c}{No infiltration} & \\multicolumn{4}{c}{Green--Ampt} \\\\", "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
     "Rainfall & LF m & HP m & LF h & HP h & LF m & HP m & LF h & HP h \\\\", "\\midrule"]
for c in C:
    row = [LAB[c]]
    for inf in ("none", "GA"):
        row += [f(st.loc[(c, m, inf)].dpeak) for m in ("LISFLOOD-FP", "HiPIMS")]
        row += [f"{int(st.loc[(c, m, inf)]['dt']):+d}".replace("-", "$-$") for m in ("LISFLOOD-FP", "HiPIMS")]
    L.append(" & ".join(row) + " \\\\")
L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
open(M + "table_stage.tex", "w").write("\n".join(L) + "\n"); print("written", M + "table_stage.tex")

# Sensitivity runs: Ks halved and doubled (radar, GraphCast) and radar aggregated to 6 h steps
s = sc[sc.width == "full"].groupby(["case", "model", "infiltration"])[["peak_ratio", "beta", "nse", "dt_peak"]].mean()
sp = st.reset_index().set_index(["case", "model", "inf"])
wbi = wb.set_index(["case", "infiltration"])
ROWS = [("R_1km_5min", "none", "Radar, 5\\,min", "none"), ("R_1km_5min", "GA_ks05", "", "$0.5K_s$"), ("R_1km_5min", "GA", "", "$K_s$"),
        ("R_1km_5min", "GA_ks2", "", "$2K_s$"), ("R_1km_6h", "none", "Radar, 6\\,h", "none"), ("R_1km_6h", "GA", "", "$K_s$"),
        ("F_graphcast", "none", "GraphCast", "none"), ("F_graphcast", "GA_ks05", "", "$0.5K_s$"), ("F_graphcast", "GA", "", "$K_s$"),
        ("F_graphcast", "GA_ks2", "", "$2K_s$"), ("F_graphcast", "GA_sm", "", "wet start"),
        ("R_1km_5min", "GA_sm", "Radar, 5\\,min", "wet start"), ("F_aifs", "GA_sm", "AIFS", "wet start"), ("F_fuxi", "GA_sm", "FuXi", "wet start")]
L = ["\\begin{table*}[t]", "\\caption{Sensitivity runs: radar and GraphCast with the saturated hydraulic conductivity halved and doubled, the radar aggregated to 6\\,h steps (the step of the forecasts), and a wet start in which the initial soil moisture is taken from ERA5-Land at 13 July 2021 18:00\\,UTC instead of SMAP. Discharge scores are gauge means for the station section; Inf.\\ is the infiltrated fraction of rainfall (LISFLOOD-FP); $\\Delta H$ is the simulated minus observed peak stage at Walferdange (m). Each cell gives LISFLOOD-FP / HiPIMS.}",
     "\\label{tab:sens}", "\\centering", "\\begin{tabular}{llrlllll}", "\\toprule",
     "Rainfall & Infiltration & Inf.\\ (\\%) & Peak & $\\beta$ & NSE & PTE (h) & $\\Delta H$ \\\\", "\\midrule"]
for k, (c, inf, lab, il) in enumerate(ROWS):
    if k in (4, 6, 11): L.append("\\midrule")
    w_ = wbi.loc[(c, inf)]; a, b = s.loc[(c, "LISFLOOD-FP", inf)], s.loc[(c, "HiPIMS", inf)]
    pair = lambda x, y, fmt="{:.2f}": (fmt.format(x + 0.0) + "/" + fmt.format(y + 0.0)).replace("-", "$-$")
    L.append(" & ".join([lab, il, f"{100 * w_.infiltrated / w_.rain:.0f}", pair(a.peak_ratio, b.peak_ratio), pair(a.beta, b.beta),
                         pair(a.nse, b.nse), pair(a.dt_peak, b.dt_peak, "{:+.1f}"),
                         pair(sp.loc[(c, "LISFLOOD-FP", inf)].dpeak, sp.loc[(c, "HiPIMS", inf)].dpeak)]) + " \\\\")
L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
open(M + "table_sens.tex", "w").write("\n".join(L) + "\n"); print("written", M + "table_sens.tex")

# Computing time per 96 h run, from the job timers (data/run_times.csv)
rt = pd.read_csv(D + "run_times.csv"); rt["min"] = rt.wall_s / 60
g = rt.groupby(["model", "infiltration", "gpus"])["min"].agg(["mean", "min", "max", "count"])
L = ["\\begin{table}[t]", "\\caption{Wall clock time of the 96\\,h simulations on NVIDIA A100 GPUs (minutes), double precision on the 10\\,m grid.}",
     "\\label{tab:time}", "\\centering", "\\begin{tabular}{llrrrr}", "\\toprule", "Model & Infiltration & GPUs & Runs & Mean & Range \\\\", "\\midrule"]
for (m, inf, n), r in g.iterrows():
    L.append(f"{m} & {'none' if inf == 'none' else 'Green--Ampt'} & {n} & {int(r['count'])} & {r['mean']:.1f} & {r['min']:.1f}--{r['max']:.1f} \\\\")
L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
open(M + "table_time.tex", "w").write("\n".join(L) + "\n"); print("written", M + "table_time.tex")
