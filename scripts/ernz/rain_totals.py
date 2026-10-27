"""Rainfall totals of every Ernz forcing over the domain, the two valleys and the Larochette catchment (masks in domain/masks.npz).
usage: python rain_totals.py <out csv> <case ...>"""
import sys, numpy as np, netCDF4 as nc, pandas as pd
M = dict(np.load("domain/masks.npz")); act = M["domain"]; NR, NC = act.shape
rows = []
for c in sys.argv[2:]:
    d = nc.Dataset(f"forcing/rain_{c}.nc"); v = d["rainfall_depth"]; tot = np.zeros((NR, NC)); ser = []
    for i in range(v.shape[0]):
        fr = np.array(v[i]); fr = fr.T if fr.shape == (NC, NR) else fr; tot += fr; ser.append(float(fr[act].mean()))
    dt = float(d["time"][1] - d["time"][0]); ser = np.array(ser)
    if dt < 1: step = int(round(1 / dt)); n = (len(ser) - 1) // step * step; hourly = np.add.reduceat(ser[:n], np.arange(0, n, step))
    else: hourly = ser
    rows.append({"case": c, **{k: round(float(tot[m].mean()), 1) for k, m in M.items()}, "max_interval_mm": round(float(hourly.max()), 1),
                 "interval_h": max(dt, 1.0), "hour_of_max": int(hourly.argmax()) * max(dt, 1.0), "intervals_gt_1mm": int((hourly > 1).sum())})
df = pd.DataFrame(rows); df.to_csv(sys.argv[1], index=False); print(df.to_string(index=False))
