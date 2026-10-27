"""Water balance of Ernz runs from the LISFLOOD-FP mass file (columns Time Tstep MinTstep NumTsteps Area Vol Qin Hds Qout Qerror Verror Rain-(Inf+Evap)).
Rain volume = domain mean of the forcing (rain_totals.csv) x domain area; infiltration = -(last column) at the end; outflow = integral of Qout.
usage: python wb_ernz.py <runs dir> <rain_totals.csv> <out csv>   (run on the HPC login node)"""
import sys, re, glob, os, numpy as np, csv
RUNS, RT, OUT = sys.argv[1:4]; A = 203.1e6
tot = {r["case"]: float(r["domain"]) for r in csv.DictReader(open(RT))}
rows = []
for d in sorted(glob.glob(RUNS + "/E18_*")):
    m = glob.glob(d + "/results/*.mass")
    if not m or not os.path.exists(d + "/time.txt") or "Elapsed" not in open(d + "/time.txt").read(): continue
    a = np.loadtxt(m[0], skiprows=1); t = a[:, 0]; vol = a[:, 5]; qout = a[:, 8]; col = a[:, 11]
    mm = re.match(r"E18_(.+)_(noinf|ga)(_ks\w+)?$", os.path.basename(d)); case = mm.group(1) + (mm.group(3) or ""); infl = mm.group(2)
    rain = tot.get(re.sub(r"_(dry|ks\w+)$", "", case), np.nan) * 1e-3 * A
    # last mass column: -Inf in the binary before 25 Oct 2026 07:42 (first 20 event runs), Rain-Inf after (the _ks runs and the regenerated radar GA run; identical hydrographs)
    inf = (rain - col) if (col[-1] > 0 or mm.group(3)) else -col
    out = (np.trapezoid if hasattr(np, "trapezoid") else np.trapz)(qout, t); dv = vol[-1] - vol[0]
    rows.append(dict(case=case, infiltration=infl, rain_Mm3=rain / 1e6, infiltrated_Mm3=inf[-1] / 1e6, infiltrated_frac=inf[-1] / rain if rain else np.nan,
                     outflow_Mm3=out / 1e6, storage_change_Mm3=dv / 1e6, start_storage_Mm3=vol[0] / 1e6, end_storage_Mm3=vol[-1] / 1e6, peak_outflow=qout.max(),
                     closure_Mm3=(rain - inf[-1] - out - dv) / 1e6, sim_h=t[-1] / 3600))
with open(OUT, "w") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); [w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}) for r in rows]
for r in rows: print("%-26s %-5s rain %6.2f inf %6.2f (%4.0f %%) out %6.2f dS %6.2f closure %7.3f Mm3" % (r["case"], r["infiltration"], r["rain_Mm3"], r["infiltrated_Mm3"], 100 * r["infiltrated_frac"], r["outflow_Mm3"], r["storage_change_Mm3"], r["closure_Mm3"]))
