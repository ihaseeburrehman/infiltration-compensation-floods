"""Green-Ampt v3a: Ks, porosity and suction head from SoilGrids 2.0 (0-30 cm) with Saxton and Rawls (2006)
pedotransfer functions, including bulk density and coarse fragments. Impervious cells, initially wet cells and the
v2 wet-area mask are kept from v2; initial moisture is the same satellite product (capped at porosity)."""
import numpy as np, os, rasterio
from rasterio.warp import reproject, Resampling
from affine import Affine
S = "soil/"; O = "/ec/res4/scratch/lux0804/lf84_cases/inputs_v3a"; os.makedirs(O, exist_ok=True)
NR, NC, X0, Y0, C = 4362, 2440, 59799.282445142424, 49254.0, 10.0
tf = Affine(C, 0, X0, 0, -C, Y0 + NR * C)
DZ = {"0-5cm": 5, "5-15cm": 10, "15-30cm": 15}
def rd(v, d):
    with rasterio.open(S + f"sg_{v}_{d}.tif") as r:
        a = r.read(1).astype("float64"); a[a < 0] = np.nan; return a, r.transform, r.crs
def saxton(S_, C_, OM, BD, Rv):
    """Saxton and Rawls (2006), eqs 1-3, 5, 6, 16-19, 22: S_, C_ fractions, OM %, BD g/cm3, Rv gravel volume fraction."""
    t1500 = -0.024*S_ + 0.487*C_ + 0.006*OM + 0.005*S_*OM - 0.013*C_*OM + 0.068*S_*C_ + 0.031
    t1500 = t1500 + (0.14*t1500 - 0.02)
    t33 = -0.251*S_ + 0.195*C_ + 0.011*OM + 0.006*S_*OM - 0.027*C_*OM + 0.452*S_*C_ + 0.299
    t33 = t33 + (1.283*t33**2 - 0.374*t33 - 0.015)
    ts33 = 0.278*S_ + 0.034*C_ + 0.022*OM - 0.018*S_*OM - 0.027*C_*OM - 0.584*S_*C_ + 0.078
    ts33 = ts33 + (0.636*ts33 - 0.107)
    ths = t33 + ts33 - 0.097*S_ + 0.043
    # density adjustment with the SoilGrids bulk density (eqs 6-9)
    thsDF = 1 - BD/2.65; t33DF = t33 - 0.2*(ths - thsDF); ts33DF = thsDF - t33DF
    t1500 = np.minimum(t1500, t33DF - 0.01)
    B = (np.log(1500) - np.log(33)) / (np.log(t33DF) - np.log(t1500)); lam = 1/B
    Ks = 1930 * np.clip(thsDF - t33DF, 1e-4, None)**(3 - lam)          # mm/h, matrix soil
    pe = -21.67*S_ - 27.93*C_ - 81.97*ts33DF + 71.12*S_*ts33DF + 8.29*C_*ts33DF + 14.05*S_*C_ + 27.16
    pe = np.clip(pe + (0.02*pe**2 - 0.113*pe - 0.70), 0.5, None)        # air entry, kPa
    psi_f = (2 + 3*lam)/(1 + 3*lam) * pe/2 * 10.2                        # wetting-front suction, cm (Rawls et al. 1983)
    a = BD/2.65; Rw = Rv/(a + Rv*(1 - a))
    Kb = Ks*(1 - Rw)/(1 - Rw*(1 - 1.5*a))                                # gravel correction (eq 22)
    return Kb, thsDF*(1 - Rv), psi_f
L = {}
for d, z in DZ.items():
    sa, T, crs = rd("sand", d); cl, _, _ = rd("clay", d); soc, _, _ = rd("soc", d); bd, _, _ = rd("bdod", d); cf, _, _ = rd("cfvo", d)
    L[d] = saxton(sa/1000, cl/1000, soc/100*1.724, bd/100, np.nan_to_num(cf)/1000)
w = np.array(list(DZ.values()), float)
Ks = w.sum() / sum(wi/L[d][0] for wi, d in zip(w, DZ))                  # harmonic (vertical) mean over 0-30 cm
por = sum(wi*L[d][1] for wi, d in zip(w, DZ)) / w.sum()
psi = sum(wi*L[d][2] for wi, d in zip(w, DZ)) / w.sum()
print("SoilGrids 250 m, basin box: Ks p10/50/90 %s mm/h, porosity %s, psi %s cm" % tuple(np.round(np.nanpercentile(x, [10, 50, 90]), 3) for x in (Ks, por, psi)))
def to10(a):
    o = np.full((NR, NC), np.nan, "float32")
    reproject(a.astype("float32"), o, src_transform=T, src_crs=crs, dst_transform=tf, dst_crs="EPSG:2169",
              resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan); return o
K10, P10, F10 = to10(Ks), to10(por), to10(psi)
ga = dict(np.load(S + "green_ampt_10m_v2.npz"))
b = np.load("/ec/res4/scratch/lux0804/radar_scale_2026/basin_mask.npy").astype(bool)
act = (ga["Ks"] > 0) & np.isfinite(K10) & np.isfinite(P10)
print("active (pervious, not wet) %.1f%% of basin; v2 Ks there: mean %.2f, v3a mean %.2f mm/h" % (100*act[b].mean(), ga["Ks"][act & b].mean()*3.6e6, K10[act & b].mean()))
for c in np.unique(ga["cls"][act & b]):
    m = act & b & (ga["cls"] == c); print("   class %2d  %4.1f%%  Ks v2 %6.2f  v3a %6.2f mm/h" % (c, 100*m.sum()/b.sum(), ga["Ks"][m].mean()*3.6e6, K10[m].mean()))
ga["Ks"][act] = K10[act]/3.6e6; ga["porosity"][act] = P10[act]; ga["psi"][act] = F10[act]/100
ga["theta_i"][act] = np.minimum(ga["theta_i"][act], ga["porosity"][act])
ga["dtheta"] = np.where(ga["Ks"] > 0, np.clip(ga["porosity"] - ga["theta_i"], 0, None), 0).astype("float32")
print("basin-mean Ks v3a %.2f mm/h; mean dtheta %.3f" % (ga["Ks"][b].mean()*3.6e6, ga["dtheta"][b].mean()))
np.savez_compressed(S + "green_ampt_10m_v3a.npz", **ga)
with open("/ec/res4/scratch/lux0804/radar_scale_2026/inputs/Alzette_sub_basin_10m_bridge_burn.asc") as f: hdr = "".join(f.readline() for _ in range(6))
for name, arr, fmt in [("ga_ks_mmh", ga["Ks"]*3.6e6, "%.5g"), ("ga_psi_mm", ga["psi"]*1000, "%.5g"),
                       ("ga_porosity", ga["porosity"], "%.4f"), ("ga_theta_init", ga["theta_i"], "%.4f")]:
    np.savetxt(O + f"/{name}.asc", arr, fmt=fmt, header=hdr.rstrip("\n"), comments="")
print("written", O)
