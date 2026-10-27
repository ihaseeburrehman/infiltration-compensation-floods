"""Table of the transect width convergence: peak discharge as a fraction of the full section value against the total
width of the section centred on the gauge, radar forcing, LISFLOOD-FP, both infiltration settings; Alzette and Ernz gauges.
usage: python make_width_table.py <paper_v4 dir>"""
import sys, pandas as pd
P = sys.argv[1].rstrip("/")
al = pd.read_csv(P + "/data/width_convergence.csv"); er = pd.read_csv(P + "/../ernz2018/results/width_convergence.csv")
W = ["50", "100", "200", "300", "400", "600", "800", "full"]
rows = [("Alzette", "Pfaffenthal", "L84_R_1km_5min", al, "400"), ("Alzette", "Hesperange", "L84_R_1km_5min", al, "400"),
        ("Alzette", "Steinsel", "L84_R_1km_5min", al, "700"), ("Alzette", "Livange", "L84_R_1km_5min", al, "400$^{a}$"),
        ("Ernz", "Larochette", "R_1km_5min", er, "600"), ("Ernz", "Mullerthal", "R_1km_5min", er, "400")]
L = ["\\begin{table*}[t]", "\\caption{Convergence of the transect discharge with the width of the section. Peak discharge through a section of the given total width centred on the gauge, as a fraction of the peak through the whole transect, for the radar forcing in LISFLOOD-FP without and with Green--Ampt infiltration (none / GA), the width of the station section used for the scores, and the wetted width of the whole transect at the time of the peak without infiltration. $^{a}$At Livange the transect spans 1100\\,m. The step between 600 and 800\\,m is a separate ditch 350\\,m south east of the station, and the station section stops before it.}",
     "\\label{tab:width}", "\\centering", "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{llrrrrrrrrr}", "\\toprule",
     "Case & Gauge & " + " & ".join(f"{w}\\,m" for w in W[:-1]) + " & Used (m) & Wet (m) \\\\", "\\midrule"]
for case, g, run, d, full in rows:
    cells = []
    for w in W[:-1]:
        a = d[(d.run == run + "_noinf") & (d.gauge == g) & (d.width == w)].peak_frac.iloc[0]; b = d[(d.run == run + "_ga") & (d.gauge == g) & (d.width == w)].peak_frac.iloc[0]
        cells.append(f"{a:.2f} / {b:.2f}")
    wet = d[(d.run == run + "_noinf") & (d.gauge == g) & (d.width == "full")].wet_width.iloc[0]
    L.append(f"{case} & {g} & " + " & ".join(cells) + f" & {full} & {wet:.0f} \\\\")
L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
open(P + "/manuscript/table_width.tex", "w").write("\n".join(L) + "\n"); print("\n".join(L))
