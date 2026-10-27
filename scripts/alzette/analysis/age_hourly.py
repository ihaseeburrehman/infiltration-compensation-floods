"""AGE native records -> hourly discharge on model hours h=0..96 (h0 = 13 Jul 2021 18:00 UTC).
Two clock assumptions: AGE clock = UTC ('utc') or AGE clock = UTC+1 ('utc1', as the manuscript states).
The AGE clock has no summer time (year file of Pfaffenthal 2021: 96 values on 28 Mar and 31 Oct) and matches CAMELS-LUX (UTC+1).
15 min records are averaged over the hour centred on h. The hourly files (Q60) give the mean of the hour ending at the time
stamp (they match CAMELS-LUX one hour later), so the hour centred on h is the mean of the two hourly values ending at h and h+1."""
import sys, csv, datetime as dt, numpy as np
D = sys.argv[1]; OUT = sys.argv[2]
F = {"Pfaffenthal": "Pfaffenthal.Q15.VO.07.2021.csv", "Steinsel": "Steinsel.Q15.VO.07.2021.csv",
     "Hesperange": "Hesperange.Q60_07.2021.csv", "Livange": "Livange.Q60_07.2021.csv"}
H0 = dt.datetime(2021, 7, 13, 18)
rows = []
for g, f in F.items():
    t, q = [], []
    for line in open(f"{D}/{f}", encoding="latin-1"):
        p = line.strip().split(";")
        if len(p) == 3 and p[0][:1].isdigit():
            try: v = float(p[2].replace(",", "."))
            except ValueError: continue
            t.append(dt.datetime.strptime(p[0] + p[1], "%d.%m.%y%H:%M:%S")); q.append(v)
    t = np.array([(x - H0).total_seconds() / 3600 for x in t]); q = np.array(q)
    for lab, off in (("utc", 0), ("utc1", 1), ("utc2", 2)):
        tt = t - off
        for h in range(0, 97):
            if "Q60" in f:                                 # hourly mean ending at tt: half of the hours ending at h and h+1
                m = (np.abs(tt - h) < 1e-6) | (np.abs(tt - h - 1) < 1e-6)
                rows.append((g, lab, h, q[m].mean() if m.sum() == 2 else np.nan))
            else:
                m = (tt > h - 0.5) & (tt <= h + 0.5)      # mean of records within the hour centred on h
                rows.append((g, lab, h, q[m].mean() if m.any() else np.nan))
    i = q[(t >= 0) & (t <= 96)].argmax(); print(g, "native peak", q[(t >= 0) & (t <= 96)][i], "at clock h", t[(t >= 0) & (t <= 96)][i])
with open(OUT, "w") as fh:
    w = csv.writer(fh); w.writerow(["gauge", "clock", "h", "Q"]); w.writerows(rows)
