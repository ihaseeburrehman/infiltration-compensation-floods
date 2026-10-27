"""AGE water level at Walferdange (15 min, cm above the gauge datum, clock UTC+1 all year) -> hourly mean level
as stage in m above the gauge zero (and in m above sea level) on model hours h = 0..96 (h0 = 13 July 2021 18:00 UTC). Gauge datum 225.32 m (file header MPH).
usage: stage_obs.py <AGE dir> <out csv>"""
import sys, csv, datetime as dt, numpy as np
D, OUT = sys.argv[1], sys.argv[2]
H0 = dt.datetime(2021, 7, 13, 18); t, z = [], []
for line in open(f"{D}/Walferdange.W15.07.2021.csv", encoding="latin-1"):
    p = line.strip().split(";")
    if len(p) == 3 and p[0][:1].isdigit():
        try: v = float(p[2].replace(",", "."))
        except ValueError: continue
        t.append((dt.datetime.strptime(p[0] + p[1], "%d.%m.%y%H:%M:%S") - H0).total_seconds() / 3600 - 1); z.append(225.32 + v / 100)
t, z = np.array(t), np.array(z)
with open(OUT, "w") as f:
    w = csv.writer(f); w.writerow(["h", "wse", "stage"])
    for h in range(97):
        m = (t > h - 0.5) & (t <= h + 0.5)
        w.writerow([h, round(z[m].mean(), 3), round(z[m].mean() - 225.32, 3)] if m.any() else [h, "", ""])
i = np.argmax(np.where((t >= 0) & (t <= 96), z, -1)); print("peak %.2f m at h %.2f (UTC)" % (z[i], t[i]))
