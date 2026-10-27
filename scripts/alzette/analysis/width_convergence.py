"""Width convergence of the transect discharge: peak discharge (and NSE at the Alzette gauges) as a function of the
total transect width centred on the gauge, from the per segment flux files (segflux_fv.npz of xs_fullvalley.py, segflux_seg.npz
of xs_ernz_seg.py). usage: python width_convergence.py <npz dir> <out csv> [obs csv with gauge,h,Q]"""
import sys, glob, os, numpy as np, pandas as pd
D, OUT = sys.argv[1], sys.argv[2]; OBS = pd.read_csv(sys.argv[3]) if len(sys.argv) > 3 else None
W = [20, 50, 100, 150, 200, 300, 400, 600, 800, 1200, None]
def qwidth(z, g, w):
    t, L, qn = z[g + "_t"], z[g + "_L"], z[g + "_qn"]; a, b = (-w / 2, w / 2) if w else (-1e9, 1e9)
    ov = np.clip(np.minimum(t + L / 2, b) - np.maximum(t - L / 2, a), 0, None); return (qn * ov).sum(1)
rows = []
for f in sorted(glob.glob(D + "/*.npz")):
    z = np.load(f); run = os.path.basename(f)[:-4]; gauges = sorted({k[:-2] for k in z.files if k.endswith("_t")})
    for g in gauges:
        qf = qwidth(z, g, None); ipk = int(qf.argmax())
        for w in W:
            q = qwidth(z, g, w); r = dict(run=run, gauge=g, width=w if w else "full", peak=round(float(q.max()), 2), peak_frac=round(float(q.max() / qf.max()), 4),
                                         q_at_full_peak_frac=round(float(q[ipk] / qf[ipk]), 4), wet_width=None)
            wd = z[g + "_wd"][ipk]; L = z[g + "_L"]; t = z[g + "_t"]
            m = (np.abs(t) <= (w / 2 if w else 1e9)) & (wd > 0.05); r["wet_width"] = round(float(L[m].sum()), 1)
            if OBS is not None and (OBS.gauge == g).any():
                o = OBS[OBS.gauge == g].set_index("h").Q; h = [k for k in z["steps"].tolist() if k in o.index and not np.isnan(o[k]) and 1 <= k <= 96]
                s = np.array([q[list(z["steps"]).index(k)] for k in h]); ob = o.loc[h].values
                r["nse"] = round(float(1 - ((s - ob) ** 2).sum() / ((ob - ob.mean()) ** 2).sum()), 3); r["peak_ratio"] = round(float(s.max() / ob.max()), 3)
            rows.append(r)
df = pd.DataFrame(rows); df.to_csv(OUT, index=False); print(df.to_string(index=False, max_rows=400))
