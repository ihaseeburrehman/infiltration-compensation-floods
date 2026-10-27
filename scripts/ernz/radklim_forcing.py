"""RADKLIM RW (hourly, 1 km, gauge adjusted with German gauges only) -> LISFLOOD-FP dynamicrainfile on the Ernz Noire 10 m grid.
Conventions as in the Alzette paper: mm per interval, labelled at the interval START, closing zero frame, block mapping of each
radar pixel onto the 10 m cells (no interpolation). RADKLIM stamps are interval ENDS (hh:50 = [hh-1:50, hh:50)).
Gap fill: RADKLIM has no data west of about 6.25 E (41 % of the domain, the upstream half). Missing pixels take the inverse
distance weighted mean (power 2) of the valid pixels within 15 km in the same hour.
Gauge correction (Mathias 2019): RADKLIM underestimates the convective core by about 2 in the peak hour; the 6 h total at
Waldbillig (100.6 mm gauge vs 60.2 mm radar) is explained by that hour alone, so only the peak hour [23:50, 00:50) is
corrected: multiplicative factor field = IDW (power 2) of gauge/radar ratios at the stations on valid pixels, clipped to
[0.3, 3]; in that hour the gap is filled from the corrected pixels plus the gauges inside the gap (Godbrange) as pseudo pixels.
5 min option (--yw <dir of raa01-yw2017.002_10000-YYMMDDhhmm-dwd---bin files>): the hourly value of each pixel (corrected, gap
filled) is distributed over its twelve 5 min YW intervals in proportion to the YW depths (as make_forcing.py of paper_v4);
hours without YW files or without YW echo are spread evenly. RADOLAN binary rows are stored south to north, the RW NetCDF rows
in the same order, so no flip. The YW file stamped hh:mm covers [hh:mm-5, hh:mm).
Gap constraint (--catchref <CAMELS-LUX 15 min csv> <catchment shapefile> <gauge_id>): CAMELS-LUX gives for the Larochette
catchment (69 km2, almost entirely inside the RADKLIM gap) the 15 min catchment mean rain RR_rad (RADKLIM where available,
interpolated AGE/ASTA/LIST stations elsewhere; stamps UTC+1, interval end). Each hour, the gap pixels are scaled by one
factor so that the forcing mean over that catchment equals the CAMELS hourly value (factor clipped to [0.1, 5]; zero fill
with CAMELS rain > 0 becomes uniform CAMELS rain). Without it the IDW fill carries the eastern core over the Ernz Blanche
(82.7 mm instead of 38.5 mm for 31 May 18 UTC + 48 h).
usage: python radklim_forcing.py <domain dir> <RW nc dir> <start YYYYmmddHH> <hours> <case> <out dir> [--nocorr] [--yw <dir>]
       [--catchref <csv> <shp> <gauge_id>]"""
import sys, json, datetime as dt, numpy as np, netCDF4 as nc, rasterio
from pyproj import Transformer
D, RW, T0, H, CASE, O = sys.argv[1:7]; CORR = "--nocorr" not in sys.argv
YW = sys.argv[sys.argv.index("--yw") + 1] if "--yw" in sys.argv else None
CREF = sys.argv[sys.argv.index("--catchref") + 1:sys.argv.index("--catchref") + 4] if "--catchref" in sys.argv else None
import re, os
def read_radolan(path):
    b = open(path, "rb").read(); k = b.index(b"\x03"); h = b[:k].decode("ascii", "ignore"); a = np.frombuffer(b[k + 1:], dtype="<i2")
    m = re.search(r"GP\s*(\d+)x\s*(\d+)", h); ny, nx = int(m.group(1)), int(m.group(2)); prec = 10.0 ** float(re.search(r"PR\s*E(-?\d+)", h).group(1))
    a = a.reshape(ny, nx); v = (a & 0x0FFF).astype("float64") * prec; v[(a & 0x2000) > 0] = np.nan; return v
T0 = dt.datetime.strptime(T0, "%Y%m%d%H"); H = int(H)
with rasterio.open(f"{D}/dem10_mean.tif") as s: dem = s.read(1); tf = s.transform
NR, NC = dem.shape; act = dem > -9000
x = tf.c + (np.arange(NC) + .5) * tf.a; y = tf.f - (np.arange(NR) + .5) * tf.a
def radar_index():
    """row*nx+col of the RADOLAN pixel (polar stereographic, 1 km) containing each 10 m cell centre"""
    d = nc.Dataset(sorted(__import__("glob").glob(f"{RW}/*/RW_2017.002_*.nc"))[0]); rx, ry = d["x"][:], d["y"][:]
    X, Y = np.meshgrid(x, y); px, py = Transformer.from_crs("EPSG:2169", d["crs"].proj4.replace("+units=km", "+units=m"), always_xy=True).transform(X.ravel(), Y.ravel())
    ci = np.round((px - rx[0]) / (rx[1] - rx[0])).astype(int); ri = np.round((py - ry[0]) / (ry[1] - ry[0])).astype(int)
    return (ri * len(rx) + ci).reshape(NR, NC)
try:
    flat = np.load(f"{D}/radklim_index.npy"); assert flat.shape == (NR, NC)
except (FileNotFoundError, AssertionError):
    flat = radar_index(); np.save(f"{D}/radklim_index.npy", flat); print("radar index rebuilt")
used = np.unique(flat[act])
# peak hour gauge values (mm, [23:50, 00:50) 31 May/1 Jun 2018 UTC), Mathias (2019) Fig. 6 and text
GAUGES = {"Waldbillig": (6.2773, 49.79806, 69.4), "Christnach": (6.27, 49.80, 73.2), "Echternach": (6.44337, 49.8031, 25.5),
          "Godbrange": (6.23569, 49.73637, 7.8), "Bettendorf": (6.2095, 49.8741, 2.2)}
PEAK = dt.datetime(2018, 6, 1, 0, 50)   # the hour assigned to [00:00, 01:00) 1 Jun
files = {}
def rr_at(tend):
    f = f"{RW}/{tend:%Y}/RW_2017.002_{tend:%Y%m}.nc"
    if f not in files: files[f] = nc.Dataset(f)
    d = files[f]; t = nc.num2date(d["time"][:], d["time"].units)
    i = [k for k, v in enumerate(t) if (v.year, v.month, v.day, v.hour, v.minute) == (tend.year, tend.month, tend.day, tend.hour, tend.minute)]
    assert len(i) == 1, (tend, len(i)); return np.ma.filled(d["RR"][i[0]].astype("float64"), np.nan), d
def idw(px, py, vx, vy, vv, p=2, rmax=15000.0):
    out = np.empty(len(px))
    for k in range(len(px)):
        dd = np.hypot(vx - px[k], vy - py[k]); m = dd < rmax
        if not m.any(): m = dd <= dd.min()
        w = 1 / np.maximum(dd[m], 500.0) ** p; out[k] = (w * vv[m]).sum() / w.sum()
    return out
frames, series, filled, yw_missing, cfac = [], [], [], [], []
if CREF:
    import pandas as pd, geopandas as gpd
    from rasterio.features import rasterize
    cam = pd.read_csv(CREF[0], parse_dates=["Date"]); cam["utc_end"] = cam.Date - pd.Timedelta(hours=1)
    cpoly = gpd.read_file(CREF[1]); cpoly = cpoly[cpoly.gauge_id == int(CREF[2])].to_crs(2169)
    cmask = rasterize(((g_, 1) for g_ in cpoly.geometry), out_shape=(NR, NC), transform=tf, fill=0, dtype="uint8").astype(bool) & act
    def cam_hour(t_start):
        w = cam[(cam.utc_end > t_start) & (cam.utc_end <= t_start + dt.timedelta(hours=1))]
        return float(w.RR_rad.sum()) if len(w) == 4 and w.RR_rad.notna().all() else np.nan
tr = Transformer.from_crs("EPSG:4326", "EPSG:2169", always_xy=True)
for h in range(H):
    tend = T0 + dt.timedelta(hours=h + 1) - dt.timedelta(minutes=10)   # RADKLIM hour hh:50 is assigned to [hh:00, hh+1:00)
    a, d = rr_at(tend)
    if h == 0:
        lon, lat = d["lon"][:].ravel(), d["lat"][:].ravel(); PX, PY = tr.transform(lon[used], lat[used])
    v = a.ravel()[used]; ok = np.isfinite(v); filled.append(float((~ok).mean()))
    gx, gy = tr.transform([g[0] for g in GAUGES.values()], [g[1] for g in GAUGES.values()]); gv = np.array([g[2] for g in GAUGES.values()])
    if CORR and tend == PEAK:
        # 1. correction of the valid radar pixels: IDW of gauge/radar ratios at the gauges that sit on valid pixels
        kk = [int(np.argmin(np.hypot(PX - xg, PY - yg))) for xg, yg in zip(gx, gy)]; onrad = np.array([ok[k] for k in kk])
        ratios = np.array([gv[i] / v[kk[i]] if onrad[i] and v[kk[i]] > 1 else np.nan for i in range(len(kk))])
        for n, i in zip(GAUGES, range(len(kk))): print(f"  {n:11s} gauge {gv[i]:5.1f} radar {'gap' if not onrad[i] else f'{v[kk[i]]:5.1f}'} ratio {ratios[i]:.2f}")
        m = np.isfinite(ratios); F = np.clip(idw(PX, PY, np.array(gx)[m], np.array(gy)[m], ratios[m], rmax=30000.0), 0.3, 3.0)
        v[ok] = v[ok] * F[ok]; print(f"  peak hour factor on valid pixels: mean {F[ok].mean():.2f}, max {F[ok].max():.2f}")
        # 2. gap fill with the corrected pixels and the gauges in the gap as pseudo pixels
        sx = np.concatenate([PX[ok], np.array(gx)[~onrad]]); sy = np.concatenate([PY[ok], np.array(gy)[~onrad]]); sv = np.concatenate([v[ok], gv[~onrad]])
        if (~ok).any(): v[~ok] = idw(PX[~ok], PY[~ok], sx, sy, sv)
    elif (~ok).any(): v[~ok] = idw(PX[~ok], PY[~ok], PX[ok], PY[ok], v[ok])
    if CREF:
        t_start = tend + dt.timedelta(minutes=10) - dt.timedelta(hours=1); ref = cam_hour(t_start)
        lut0 = np.zeros(flat.max() + 1); lut0[used] = v; cur = float(lut0[flat][cmask].mean())
        gapm = ~ok
        if np.isfinite(ref):
            if cur > 0.01: fct = float(np.clip(ref / cur, 0.1, 5.0)); v[gapm] = v[gapm] * fct
            elif ref > 0: v[gapm] = ref; fct = np.nan
            else: fct = 1.0
            cfac.append(fct)
    def to10(vals):
        lut = np.zeros(flat.max() + 1, "float32"); lut[used] = vals; g10 = lut[flat]; g10[~act] = 0; return g10
    if YW is None:
        g10 = to10(v); frames.append(g10); series.append(float(g10[act].mean()))
    else:
        t_start = tend + dt.timedelta(minutes=10) - dt.timedelta(hours=1)          # hour [hh:00, hh+1:00)
        w = []
        for k5 in range(12):
            te5 = t_start + dt.timedelta(minutes=5 * (k5 + 1)); f5 = f"{YW}/raa01-yw2017.002_10000-{te5:%y%m%d%H%M}-dwd---bin"
            w.append(np.nan_to_num(read_radolan(f5).ravel()[used]) if os.path.exists(f5) else None)
        if any(x is None for x in w): w = [np.ones(len(used)) for _ in range(12)]; yw_missing.append(tend)
        w = np.array(w); sw = w.sum(0); w = np.where(sw > 0, w / np.where(sw > 0, sw, 1), 1.0 / 12)
        for k5 in range(12):
            g10 = to10(v * w[k5]); frames.append(g10); series.append(float(g10[act].mean()))
out = f"{O}/rain_{CASE}.nc"; ds = nc.Dataset(out, "w", format="NETCDF4")
DT = 1.0 / 12 if YW else 1.0
ds.createDimension("time", len(frames) + 1); ds.createDimension("x", NC); ds.createDimension("y", NR)
vt = ds.createVariable("time", "f8", ("time",)); vt.units = "hour"; vt.long_name = "time"; vt.axis = "T"
vx = ds.createVariable("x", "f8", ("x",)); vx.units = "m"; vx.long_name = "x"; vx.axis = "X"; vx.spatial_ref = "EPSG:2169"
vy = ds.createVariable("y", "f8", ("y",)); vy.units = "m"; vy.long_name = "y"; vy.axis = "Y"; vy.spatial_ref = "EPSG:2169"
vr = ds.createVariable("rainfall_depth", "f4", ("time", "y", "x"), zlib=True, complevel=1, chunksizes=(1, NR, NC)); vr.units = "mm"; vr.long_name = "rainfall_depth"
vx[:] = x; vy[:] = y
for i, f in enumerate(frames + [np.zeros((NR, NC), "float32")]): vr[i] = f; vt[i] = i * DT
ds.setncatts({"crs": "EPSG:2169", "source": "DWD RADKLIM RW 2017.002 hourly", "time_stamp_convention": "interval_start",
              "window_start_utc": T0.strftime("%Y-%m-%dT%H:%M:%SZ"), "window_end_utc": (T0 + dt.timedelta(hours=H)).strftime("%Y-%m-%dT%H:%M:%SZ"),
              "gap_fill": "IDW p2 of valid pixels within 15 km", "gauge_correction": "peak hour [00:00,01:00) 1 Jun 2018 only (RADKLIM stamp 00:50), IDW of gauge/radar ratios clipped 0.3-3, gap filled with gauges" if CORR else "none",
              "domain_total_mm": float(sum(series)), "time_step_h": DT, "yw_5min": "none" if YW is None else f"{YW}; hours spread evenly (no YW files): {len(yw_missing)}",
              "gap_constraint": "none" if not CREF else f"CAMELS-LUX gauge {CREF[2]} hourly RR_rad; gap factors: {np.round(cfac, 2).tolist()}"}); ds.close()
json.dump({"case": CASE, "start": T0.isoformat(), "hours": H, "domain_total_mm": sum(series), "series_mm": series, "gap_fraction_pixels": filled}, open(f"{O}/rain_{CASE}.json", "w"), indent=1)
print(CASE, "dt h", round(DT, 4), "domain total %.1f mm" % sum(series), "max interval %.1f mm" % max(series), "mean gap fraction %.2f" % np.mean(filled), "hours without YW:", [t.strftime("%d %H") for t in yw_missing])
if CREF: print(" gap factors (hours with rain):", [round(f_, 2) for f_, s_ in zip(cfac, (series if DT == 1 else np.add.reduceat(series, np.arange(0, len(series), 12)))) if s_ > 0.05])
print(" series:", [round(s, 1) for s in (series if DT == 1 else np.add.reduceat(series, np.arange(0, len(series), 12)))])
