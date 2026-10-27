"""WRF 12 km rapid cycle (6 h forecasts, with 3D-VAR cycling = After_DA, or free cycling = Before_DA) -> LISFLOOD-FP
dynamicrainfile on the Ernz 10 m grid. Each wrfout_d01_<valid> file is the +6 h state of the cycle started 6 h earlier,
so RAINNC+RAINC in it is the 6 h accumulation over [valid-6h, valid). Mapped like the data driven forecasts: linear interpolation of the 12 km cell
centres to the 10 m grid (same as ai_to_lisflood_rain.py). mm per interval, labelled at the interval start, closing zero.
usage: python wrf12_forcing.py <wrfout dir> <dem asc> <start YYYY-mm-ddTHH> <hours> <out nc>"""
import sys, glob, os, datetime as dt, numpy as np, netCDF4 as nc
from pyproj import Transformer
from scipy.interpolate import griddata
WD, DEM, T0, H, OUT = sys.argv[1:6]; T0 = dt.datetime.strptime(T0, "%Y-%m-%dT%H"); H = int(H)
hdr = dict(l.split() for l in open(DEM).readlines()[:6]); NC_, NR = int(hdr["ncols"]), int(hdr["nrows"]); x0, y0, c = float(hdr["xllcorner"]), float(hdr["yllcorner"]), float(hdr["cellsize"])
x = x0 + (np.arange(NC_) + .5) * c; y = (y0 + (np.arange(NR) + .5) * c)[::-1]; X, Y = np.meshgrid(x, y)
tr = Transformer.from_crs("EPSG:4326", "EPSG:2169", always_xy=True)
frames, series = [], []
for k in range(H // 6):
    valid = T0 + dt.timedelta(hours=6 * (k + 1))
    f = f"{WD}/wrfout_d01_{valid:%Y-%m-%d_%H_%M_%S}"; d = nc.Dataset(f)
    # After_DA files carry START_DATE = valid - 6 h (the cycle start). Before_DA files (the background of the next
    # analysis) carry START_DATE = valid, but hold the same [valid-6h, valid) accumulation: in 2018 the Before_DA field
    # stamped T correlates with the After_DA field stamped T (r 0.5 to 0.7) and not with the one stamped T+6h (r 0.2),
    # and the 2021 Alzette before-DA forcing was built with this convention.
    assert d.START_DATE in ((valid - dt.timedelta(hours=6)).strftime("%Y-%m-%d_%H:%M:%S"), valid.strftime("%Y-%m-%d_%H:%M:%S")), (f, d.START_DATE)
    r = np.asarray(d["RAINNC"][0] + d["RAINC"][0], float); lon, lat = np.asarray(d["XLONG"][0]), np.asarray(d["XLAT"][0])
    m = (lat > 48.5) & (lat < 51) & (lon > 4.5) & (lon < 8); px, py = tr.transform(lon[m], lat[m])
    g = np.clip(griddata((px, py), r[m], (X, Y), method="linear", fill_value=0.0), 0, None).astype("float32")
    frames.append(g); series.append(float(g.mean())); print(f"  {valid:%Y-%m-%d %H} UTC  6 h domain mean {series[-1]:.2f} mm (src max {r[m].max():.1f})", flush=True)
ds = nc.Dataset(OUT, "w", format="NETCDF4"); ds.createDimension("time", len(frames) + 1); ds.createDimension("x", NC_); ds.createDimension("y", NR)
vt = ds.createVariable("time", "f8", ("time",)); vt.units = "hour"; vt.long_name = "time"; vt.axis = "T"
vx = ds.createVariable("x", "f8", ("x",)); vx.units = "m"; vx.long_name = "x"; vx.axis = "X"; vx.spatial_ref = "EPSG:2169"
vy = ds.createVariable("y", "f8", ("y",)); vy.units = "m"; vy.long_name = "y"; vy.axis = "Y"; vy.spatial_ref = "EPSG:2169"
vr = ds.createVariable("rainfall_depth", "f4", ("time", "y", "x"), zlib=True, complevel=4); vr.units = "mm"; vr.long_name = "rainfall_depth"
vx[:] = x; vy[:] = y
for i, g in enumerate(frames + [np.zeros((NR, NC_), "float32")]): vr[i] = g; vt[i] = 6 * i
ds.setncatts({"crs": "EPSG:2169", "source": f"WRF 12 km rapid cycle RAINNC+RAINC, {WD}", "time_stamp_convention": "interval_start",
              "window_start_utc": T0.strftime("%Y-%m-%dT%H:%M:%SZ"), "domain_total_mm": float(sum(series))}); ds.close()
print(os.path.basename(OUT), "domain total %.1f mm" % sum(series))
