"""EURADCLIM v3.0 (KNMI, OPERA composite, gauge adjusted, 2 km hourly) as a LISFLOOD-FP dynamic rain file on the 10 m domain grid.
usage: python eurad_forcing.py <domain dir> <eurad h5 root> <start YYYY-mm-ddTHH UTC> <hours> <case> <out dir> [--corr]
--corr: gauge correction of the two hours [23, 01) UTC that hold the gauge peak hour [23:50, 00:50): factor field = IDW (p2, 30 km) of
  gauge peak hour / EURADCLIM 2 h sum at the gauge pixel (a lower bound of the true ratio), clipped [0.3, 3], as radklim_forcing.py does.
File stamped HH:00 covers [HH-1:05, HH:05) UTC and is assigned to [HH-1, HH); nearest 2 km pixel; nodata -> 0 (none inside the domain).
Same NetCDF convention as radklim_forcing.py (interval start stamps, closing zero frame)."""
import sys, glob, json, datetime as dt, numpy as np, h5py, rasterio, netCDF4 as nc
from pyproj import Transformer, CRS
DOM, ROOT, T0, H, CASE, O = sys.argv[1], sys.argv[2], dt.datetime.strptime(sys.argv[3], "%Y-%m-%dT%H"), int(sys.argv[4]), sys.argv[5], sys.argv[6]; CORR = "--corr" in sys.argv
GAUGES = {"Waldbillig": (6.2773, 49.79806, 69.4), "Christnach": (6.27, 49.80, 73.2), "Echternach": (6.44337, 49.8031, 25.5), "Godbrange": (6.23569, 49.73637, 7.8), "Bettendorf": (6.2095, 49.8741, 2.2)}
CH = [dt.datetime(2018, 5, 31, 23), dt.datetime(2018, 6, 1, 0)]   # hours corrected
with rasterio.open(f"{DOM}/dem10_mean.tif") as r: tf, NR, NC = r.transform, r.height, r.width
act = np.load(f"{DOM}/domain_mask.npy")
files = {f[-15:-3]: f for f in glob.glob(f"{ROOT}/**/*.h5", recursive=True)}
with h5py.File(next(iter(files.values()))) as h: W = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in h["where"].attrs.items()}
proj = CRS.from_proj4(W["projdef"]); ll2p = Transformer.from_crs("EPSG:4326", proj, always_xy=True); lu2p = Transformer.from_crs("EPSG:2169", proj, always_xy=True)
xul, yul = ll2p.transform(W["UL_lon"], W["UL_lat"]); dx, dy = W["xscale"], W["yscale"]
cols, rows = np.meshgrid(np.arange(NC), np.arange(NR)); xc, yc = tf * (cols + 0.5, rows + 0.5); px, py = lu2p.transform(xc, yc)
J = np.floor((px - xul) / dx).astype(int); I = np.floor((yul - py) / dy).astype(int)
frames, series, nmiss = [], [], 0
F = np.ones((NR, NC), "float32"); corr_note = "none"
if CORR:
    raw = {}
    for t in CH:
        with h5py.File(files[f"{t + dt.timedelta(hours=1):%Y%m%d%H%M}"]) as f: d = f["dataset1/data1/data"][:].astype("float32"); d[d < -8000000] = 0; raw[t] = d
    two = raw[CH[0]] + raw[CH[1]]
    gx, gy, ratios = [], [], []
    for n, (lo, la, gv) in GAUGES.items():
        x_, y_ = ll2p.transform(lo, la); i_, j_ = int((yul - y_) // dy), int((x_ - xul) // dx); r2 = float(two[i_, j_])
        print(f"  {n:11s} gauge 1 h {gv:5.1f}  EURADCLIM 2 h {r2:5.1f}  ratio {gv / r2:.2f}")
        xx, yy = Transformer.from_crs("EPSG:4326", "EPSG:2169", always_xy=True).transform(lo, la); gx.append(xx); gy.append(yy); ratios.append(gv / r2)
    gx, gy, ratios = map(np.array, (gx, gy, ratios))
    dd = np.hypot(xc[..., None] - gx, yc[..., None] - gy); w = 1 / np.maximum(dd, 500.0) ** 2; w[dd > 30000] = 0
    F = np.clip((w * ratios).sum(-1) / np.maximum(w.sum(-1), 1e-12), 0.3, 3.0).astype("float32"); F[w.sum(-1) == 0] = 1.0
    print("  factor on domain: mean %.2f, min %.2f, max %.2f" % (F[act].mean(), F[act].min(), F[act].max()))
    corr_note = "hours [23,01) UTC x IDW(gauge 1 h / EURADCLIM 2 h) clipped [0.3,3]; ratios " + ", ".join(f"{n} {r:.2f}" for n, r in zip(GAUGES, ratios))
for h in range(H):
    stamp = f"{T0 + dt.timedelta(hours=h + 1):%Y%m%d%H%M}"
    if stamp not in files: nmiss += 1; g = np.zeros((NR, NC), "float32")
    else:
        with h5py.File(files[stamp]) as f: d = f["dataset1/data1/data"][:].astype("float32")
        d[d < -8000000] = 0.0; g = d[I, J]; g[~act] = 0.0
        if CORR and T0 + dt.timedelta(hours=h) in CH: g = g * F
    frames.append(g); series.append(float(g[act].mean()))
x = tf.c + (np.arange(NC) + 0.5) * tf.a; y = tf.f + (np.arange(NR) + 0.5) * tf.e
out = f"{O}/rain_{CASE}.nc"; ds = nc.Dataset(out, "w", format="NETCDF4")
ds.createDimension("time", H + 1); ds.createDimension("x", NC); ds.createDimension("y", NR)
vt = ds.createVariable("time", "f8", ("time",)); vt.units = "hour"; vt.long_name = "time"; vt.axis = "T"
vx = ds.createVariable("x", "f8", ("x",)); vx.units = "m"; vx.long_name = "x"; vx.axis = "X"; vx.spatial_ref = "EPSG:2169"
vy = ds.createVariable("y", "f8", ("y",)); vy.units = "m"; vy.long_name = "y"; vy.axis = "Y"; vy.spatial_ref = "EPSG:2169"
vr = ds.createVariable("rainfall_depth", "f4", ("time", "y", "x"), zlib=True, complevel=1, chunksizes=(1, NR, NC)); vr.units = "mm"; vr.long_name = "rainfall_depth"
vx[:] = x; vy[:] = y
for i, f in enumerate(frames + [np.zeros((NR, NC), "float32")]): vr[i] = f; vt[i] = float(i)
ds.setncatts({"crs": "EPSG:2169", "source": "EURADCLIM v3.0 KNMI (OPERA composite, gauge adjusted, 2 km, CC BY 4.0), nearest pixel", "time_stamp_convention": "interval_start",
              "window_start_utc": T0.strftime("%Y-%m-%dT%H:%M:%SZ"), "window_end_utc": (T0 + dt.timedelta(hours=H)).strftime("%Y-%m-%dT%H:%M:%SZ"),
              "eurad_hour_offset": "EURADCLIM hour [HH-1:05, HH:05) assigned to [HH-1, HH)", "domain_total_mm": float(sum(series)), "time_step_h": 1.0, "missing_hours": nmiss, "gauge_correction": corr_note}); ds.close()
json.dump({"case": CASE, "start": T0.isoformat(), "hours": H, "domain_total_mm": sum(series), "series_mm": series, "missing_hours": nmiss}, open(f"{O}/rain_{CASE}.json", "w"), indent=1)
print(CASE, "domain total %.1f mm, max hour %.1f mm, missing hours %d" % (sum(series), max(series), nmiss))
