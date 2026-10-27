"""Channel geometry at the five Alzette gauges from the 0.5 m LiDAR terrain model of 2024.

usage: python lidar_sections.py <MNT_Lidar2024.tif> <Alzette 10 m model DEM .asc> <transects400.json> <OUTDIR>
Profiles follow the gauge transects of the discharge extraction (AGE station coordinates, normal to the flow),
+-200 m, sampled every 0.25 m (bilinear) on the LiDAR and by nearest cell on the 10 m model terrain.
LiDAR does not see below water: the lowest flat stretch is the water surface on the flight day.
Reports the flight-day water width, the crest of the lower bank, the width at that crest (bankfull),
the model channel width and the model bed, and at Walferdange the width at the observed peak level."""
import sys, json, numpy as np, rasterio
from rasterio.windows import from_bounds

TIF, ASC, TJ, OUT = sys.argv[1:5]
MPH = {"Livange": 265.86, "Hesperange": 255.74, "Pfaffenthal": 235.25, "Walferdange": 225.32, "Steinsel": 222.26}
OBS_PEAK_WSE = {"Walferdange": 225.32 + 3.528}
T = json.load(open(TJ))
s = np.arange(-200, 200.01, 0.25)
hdr = {}
with open(ASC) as f:
    for _ in range(6):
        k, v = f.readline().split(); hdr[k.lower()] = float(v)
X0, CS = hdr["xllcorner"], hdr["cellsize"]; YTOP = hdr["yllcorner"] + hdr["nrows"] * CS

def model_rows(rows):
    out = {}
    with open(ASC) as f:
        for _ in range(6): f.readline()
        for i, line in enumerate(f):
            if i in rows: out[i] = np.array(line.split(), float)
            if i > max(rows): break
    return out

def stretch(z, i, lim):
    a = b = i
    while a > 0 and z[a - 1] <= lim: a -= 1
    while b < len(z) - 1 and z[b + 1] <= lim: b += 1
    return a, b

res = {}
with rasterio.open(TIF) as src:
    for g, v in T.items():
        x0, y0, ux, uy = v["x"], v["y"], v["ux"], v["uy"]; px, py = -uy, ux
        xs, ys = x0 + s * px, y0 + s * py
        win = from_bounds(x0 - 210, y0 - 210, x0 + 210, y0 + 210, src.transform)
        a = src.read(1, window=win).astype(float); a[a == src.nodata] = np.nan
        tr = src.window_transform(win)
        fc = (xs - tr.c) / 0.5 - 0.5; fr = (tr.f - ys) / 0.5 - 0.5
        c = np.floor(fc).astype(int); r = np.floor(fr).astype(int); dc, dr = fc - c, fr - r
        z = (a[r, c] * (1 - dc) * (1 - dr) + a[r, c + 1] * dc * (1 - dr) + a[r + 1, c] * (1 - dc) * dr + a[r + 1, c + 1] * dc * dr)
        rows = np.floor((YTOP - ys) / CS).astype(int); cols = np.floor((xs - X0) / CS).astype(int)
        M = model_rows(set(rows.tolist())); zm = np.array([M[rr][cc] for rr, cc in zip(rows, cols)])
        near = np.abs(s) <= 25
        i = np.flatnonzero(near)[np.argmin(z[near])]; wl = z[i]
        a0, b0 = stretch(z, i, wl + 0.15)                       # flight-day water
        L = z[max(a0 - 120, 0):a0 + 1]; R = z[b0:b0 + 121]      # 30 m beyond each water edge
        crest_l, crest_r = float(np.nanmax(L)), float(np.nanmax(R)); crest = min(crest_l, crest_r)
        a1, b1 = stretch(z, i, crest)
        im = np.flatnonzero(near)[np.argmin(zm[near])]; am, bm = stretch(zm, im, zm[im] + 0.3)
        d = {"gauge_zero": MPH[g], "water_surface": round(float(wl), 2),
             "water_width": (b0 - a0) * 0.25, "bank_crest_left": round(crest_l, 2), "bank_crest_right": round(crest_r, 2),
             "bankfull_depth_above_water": round(crest - float(wl), 2), "bankfull_width": (b1 - a1) * 0.25,
             "model_bed": round(float(zm[im]), 2), "model_channel_width": round((bm - am) * 0.25, 1),
             "model_bank_crest": round(float(min(np.nanmax(zm[max(am - 120, 0):am + 1]), np.nanmax(zm[bm:bm + 121]))), 2)}
        if g in OBS_PEAK_WSE:
            lv = OBS_PEAK_WSE[g]; a2, b2 = stretch(z, i, lv)
            d.update(obs_peak_wse=round(lv, 2), width_at_obs_peak=(b2 - a2) * 0.25,
                     obs_peak_above_bank_crest=round(lv - crest, 2), wet_from_to=[float(s[a2]), float(s[b2])])
        res[g] = d
        np.savez_compressed(f"{OUT}/profile_{g}.npz", s=s, z_lidar=z, z_model=zm)
json.dump(res, open(f"{OUT}/lidar_sections.json", "w"), indent=1)
for g, d in res.items(): print(g, d)
