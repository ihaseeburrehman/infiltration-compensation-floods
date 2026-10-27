"""soil_ernz.py with KS_MODE env: v3b (default), x2, x10, deep (Ks from sl4, the 30 cm layer, only), deepwet (deep and wet units pervious)."""
"""Green-Ampt v3b parameters on the Ernz Noire 10 m grid, same recipe as the Alzette 2021 case (soil_params.py,
soil_wetmask.py, soil_v3a.py, soil_v3b.py, soil_sm.py), in one script:
  classes   INSPIRE soil map 1:100,000 (Geoportail, CC0): impervious (urban, lakes) and wet units (alluvium, colluvium
            gleyed/seepage: 'alluvions', 'suintement', 'fortement') get Ks = 0 and theta_i = porosity.
  Ks, THS   EU-SoilHydroGrids 250 m (Toth et al. 2017), sl1-sl4 (0, 5, 15, 30 cm), trapezoid weights; Ks harmonic mean;
            THS reduced by SoilGrids coarse fragments (mean of 0-30 cm).
  psi       Saxton and Rawls (2006) wetting front suction from SoilGrids sand, clay, SOC, bulk density, cfvo (0-30 cm).
  theta_i   ERA5-Land 0-28 cm relative saturation at the start time (H-TESSEL saturation by soil type) x porosity.
usage: python soil_ernz.py <dem asc or tif> <soil dir> <start YYYY-mm-ddTHH> <outdir>
  soil dir holds so.soil-100000.gml, esh_{KS,THS}_M_sl{1..4}_ernz.tif, sg/sg_{var}_{depth}.tif, era5land_ernz_*.nc, slt_ernz.nc"""
import sys, glob, os, datetime as dt, numpy as np, geopandas as gpd, rasterio, netCDF4
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling
from affine import Affine
from pyproj import Transformer
DEM, S, T0, O = sys.argv[1:5]; MODE = os.environ.get("KS_MODE", "v3b"); T0 = dt.datetime.strptime(T0, "%Y-%m-%dT%H"); os.makedirs(O, exist_ok=True)
if DEM.endswith(".tif"):
    with rasterio.open(DEM) as r: dem = r.read(1); tf = r.transform
    NR, NC = dem.shape; X0, Y0, C = tf.c, tf.f - NR * tf.a, tf.a
else:
    hdr = dict(l.split() for l in open(DEM).readlines()[:6]); NC, NR = int(hdr["ncols"]), int(hdr["nrows"]); X0, Y0, C = float(hdr["xllcorner"]), float(hdr["yllcorner"]), float(hdr["cellsize"])
    tf = Affine(C, 0, X0, 0, -C, Y0 + NR * C); dem = np.loadtxt(DEM, skiprows=6, dtype="float32")
act = dem > -9000
def to10(path, scale, nod, res=Resampling.nearest):
    with rasterio.open(path) as r:
        a = r.read(1).astype("float32"); a[(a == nod) | (a < 0)] = np.nan; a *= scale
        o = np.full((NR, NC), np.nan, "float32")
        reproject(a, o, src_transform=r.transform, src_crs=r.crs, dst_transform=tf, dst_crs="EPSG:2169", resampling=res, src_nodata=np.nan, dst_nodata=np.nan); return o
# 1. soil map classes
g = gpd.read_file(S + "/so.soil-100000.gml").to_crs(2169); l = g.soilBodyLabel.str.lower()
imp = l.str.contains("bâtie") | l.str.contains("urbaine") | l.str.startswith("lac")
wet = l.str.contains("alluvions") | l.str.contains("suintement") | l.str.contains("fortement")
cls = rasterize(((geom, 1) for geom in g.geometry[imp]), out_shape=(NR, NC), transform=tf, fill=0, dtype="uint8") == 1
wetm = rasterize(((geom, 1) for geom in g.geometry[wet]), out_shape=(NR, NC), transform=tf, fill=0, dtype="uint8") == 1
mapped = rasterize(((geom, 1) for geom in g.geometry), out_shape=(NR, NC), transform=tf, fill=0, dtype="uint8") == 1
print("soil map: mapped %.1f %% of domain, impervious %.1f %%, wet units %.1f %%" % tuple(100 * m[act].mean() for m in (mapped, cls, wetm)))
# 2. Ks and porosity (v3b)
W = np.array([2.5, 7.5, 12.5, 7.5])
KS = [to10(S + f"/esh_KS_M_sl{i}_ernz.tif", 0.01 * 10 / 24, -32767) for i in (1, 2, 3, 4)]      # cm/d x100 -> mm/h
TH = [to10(S + f"/esh_THS_M_sl{i}_ernz.tif", 0.01, 255) for i in (1, 2, 3, 4)]
CF = {d: to10(S + f"/sg/sg_cfvo_{d}.tif", 0.001, -32768) for d in ("0-5cm", "5-15cm", "15-30cm")}
cf = np.nanmean(list(CF.values()), axis=0)
Ks = W.sum() / sum(w / k for w, k in zip(W, KS))
if MODE == "x2": Ks = 2 * Ks
elif MODE == "x10": Ks = 10 * Ks
elif MODE in ("deep", "deepwet"): Ks = KS[3]                                               # 30 cm layer only (first layer without the class floor)
if MODE == "deepwet": wetm = np.zeros_like(wetm)                                        # wet units treated as pervious
print("KS_MODE", MODE); por = sum(w * t for w, t in zip(W, TH)) / W.sum() * (1 - np.nan_to_num(cf))
# 3. suction head (v3a, Saxton and Rawls 2006)
DZ = {"0-5cm": 5, "5-15cm": 10, "15-30cm": 15}
def saxton(S_, C_, OM, BD, Rv):
    t1500 = -0.024*S_ + 0.487*C_ + 0.006*OM + 0.005*S_*OM - 0.013*C_*OM + 0.068*S_*C_ + 0.031; t1500 = t1500 + (0.14*t1500 - 0.02)
    t33 = -0.251*S_ + 0.195*C_ + 0.011*OM + 0.006*S_*OM - 0.027*C_*OM + 0.452*S_*C_ + 0.299; t33 = t33 + (1.283*t33**2 - 0.374*t33 - 0.015)
    ts33 = 0.278*S_ + 0.034*C_ + 0.022*OM - 0.018*S_*OM - 0.027*C_*OM - 0.584*S_*C_ + 0.078; ts33 = ts33 + (0.636*ts33 - 0.107)
    ths = t33 + ts33 - 0.097*S_ + 0.043; thsDF = 1 - BD/2.65; t33DF = t33 - 0.2*(ths - thsDF); ts33DF = thsDF - t33DF
    t1500 = np.minimum(t1500, t33DF - 0.01); B = (np.log(1500) - np.log(33)) / (np.log(t33DF) - np.log(t1500)); lam = 1/B
    pe = -21.67*S_ - 27.93*C_ - 81.97*ts33DF + 71.12*S_*ts33DF + 8.29*C_*ts33DF + 14.05*S_*C_ + 27.16
    pe = np.clip(pe + (0.02*pe**2 - 0.113*pe - 0.70), 0.5, None)
    return (2 + 3*lam)/(1 + 3*lam) * pe/2 * 10.2                      # cm
psi_l = []
for d, z in DZ.items():
    sa = to10(S + f"/sg/sg_sand_{d}.tif", 0.001, -32768, Resampling.bilinear); cl = to10(S + f"/sg/sg_clay_{d}.tif", 0.001, -32768, Resampling.bilinear)
    soc = to10(S + f"/sg/sg_soc_{d}.tif", 0.01 * 1.724, -32768, Resampling.bilinear); bd = to10(S + f"/sg/sg_bdod_{d}.tif", 0.01, -32768, Resampling.bilinear)
    psi_l.append(z * saxton(sa, cl, soc, bd, np.nan_to_num(CF[d])))
psi = sum(psi_l) / sum(DZ.values()) / 100                             # m
# 4. initial moisture from ERA5-Land
SAT = {1: 0.403, 2: 0.439, 3: 0.430, 4: 0.520, 5: 0.614, 6: 0.766, 7: 0.472}
f = sorted(glob.glob(S + "/era5land_ernz_*.nc"))
if f:
    a = netCDF4.Dataset(f[0]); tv = [v for v in ("valid_time", "time") if v in a.variables][0]; t = netCDF4.num2date(a[tv][:], a[tv].units)
    k = [i for i, d_ in enumerate(t) if (d_.year, d_.month, d_.day, d_.hour) == (T0.year, T0.month, T0.day, T0.hour)][0]
    lat, lon = np.asarray(a["latitude"][:]), np.asarray(a["longitude"][:]); sm = np.asarray((7 * a["swvl1"][k] + 21 * a["swvl2"][k]) / 28)
    st = netCDF4.Dataset(S + "/slt_ernz.nc"); sl = np.asarray(st["slt"][0]).round().astype(int)
    Srel = np.clip(sm / np.vectorize(lambda c: SAT.get(c, np.nan))(sl), 0, 1)
    x = X0 + (np.arange(NC) + .5) * C; y = Y0 + NR * C - (np.arange(NR) + .5) * C; X, Y = np.meshgrid(x, y)
    lo, la = Transformer.from_crs("EPSG:2169", "EPSG:4326", always_xy=True).transform(X, Y)
    Sg = Srel[np.abs(la[..., None] - lat).argmin(-1), np.abs(lo[..., None] - lon).argmin(-1)]
    print("ERA5-Land %s: relative saturation domain mean %.3f (%.3f-%.3f)" % (t[k], Sg[act].mean(), Sg[act].min(), Sg[act].max()))
else:
    Sg = np.full((NR, NC), np.nan, "float32"); print("ERA5-Land file missing: theta_i not written")
# 5. assemble
ok = act & np.isfinite(Ks) & np.isfinite(por) & np.isfinite(psi) & ~cls & ~wetm
Ks_o = np.where(ok, Ks, 0.0); por_o = np.where(act & np.isfinite(por), por, 0.45); psi_o = np.where(ok, psi, 0.0)
th = np.where(ok, np.minimum(por_o, Sg * por_o), por_o)
print("pervious active %.1f %% of domain; Ks p10/50/90 %s mm/h; porosity %s; psi %s m; mean dtheta %.3f" % (100 * ok[act].mean(),
      np.percentile(Ks[ok], [10, 50, 90]).round(2), np.percentile(por[ok], [10, 50, 90]).round(3), np.percentile(psi[ok], [10, 50, 90]).round(3), np.nanmean(por_o[ok] - th[ok])))
np.savez_compressed(f"{O}/green_ampt_10m_v3b.npz", Ks=Ks_o / 3.6e6, porosity=por_o, psi=psi_o, theta_i=th, impervious=cls, wetmask=wetm)
hd = f"ncols {NC}\nnrows {NR}\nxllcorner {X0}\nyllcorner {Y0}\ncellsize {C}\nNODATA_value -9999"
for name, arr, fmt in [("ga_ks_mmh", Ks_o, "%.5g"), ("ga_psi_mm", psi_o * 1000, "%.5g"), ("ga_porosity", por_o, "%.4f")] + ([("ga_theta_init", th, "%.4f")] if f else []):
    np.savetxt(f"{O}/{name}.asc", np.where(act, arr, -9999), fmt=fmt, header=hd, comments="")
print("written", O)
