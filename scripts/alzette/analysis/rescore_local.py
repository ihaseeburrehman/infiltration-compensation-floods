"""Rescore every run from data/hydrographs.csv with the observations in data/obs_utc.csv (same scores as score_all.py,
which needs the HPC flux files). Used after the half hour correction of the AGE hourly files (26 Oct 2026).
usage: python rescore_local.py <hydrographs.csv> <obs_utc.csv> <scores_out.csv>"""
import sys, csv, numpy as np, pandas as pd
HY, OB, OUT = sys.argv[1:4]
obs = {g: dict(zip(d.h, d.Q)) for g, d in pd.read_csv(OB).groupby("gauge")}
def score(qs, qo):                                   # identical to score_all.score
    h = [k for k in range(1, 97) if k in qs and not np.isnan(qo.get(k, np.nan))]
    s = np.array([qs[k] for k in h]); o = np.array([qo[k] for k in h])
    nse = 1 - ((s - o) ** 2).sum() / ((o - o.mean()) ** 2).sum()
    r = np.corrcoef(s, o)[0, 1]; beta = s.mean() / o.mean(); gam = (s.std() / s.mean()) / (o.std() / o.mean())
    kge = 1 - np.sqrt((r - 1) ** 2 + (beta - 1) ** 2 + (gam - 1) ** 2)
    ho = h[int(o.argmax())]
    return dict(peak=s.max(), peak_obs=o.max(), peak_ratio=s.max() / o.max(), t_peak=h[int(s.argmax())], dt_peak=h[int(s.argmax())] - ho,
                beta=beta, r=r, gamma=gam, nse=nse, kge=kge)
hy = pd.read_csv(HY); rows = []
for (c, m, inf, g, w), d in hy.groupby(["case", "model", "infiltration", "gauge", "width"], sort=False):
    rows.append(dict(case=c, model=m, infiltration=inf, gauge=g, width=w, **{k: round(float(v), 4) for k, v in score(dict(zip(d.h, d.Q)), obs[g]).items()}))
with open(OUT, "w") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
print(len(rows), "rows")
