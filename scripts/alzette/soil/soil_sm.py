"""Green-Ampt v3b with ERA5-Land initial soil moisture (wet start): relative saturation of the 0-28 cm layer
(swvl1 0-7 cm, swvl2 7-28 cm) at 13 July 2021 18 UTC, nearest ERA5-Land cell, H-TESSEL saturation by soil type.
theta_i = min(porosity, S * porosity) on active cells; Ks, porosity and suction head as v3b."""
import numpy as np, os, netCDF4
from pyproj import Transformer
S = "soil/"; E = S + "era5land/"; O = "/ec/res4/scratch/lux0804/lf84_cases/inputs_v3b_sm"; os.makedirs(O, exist_ok=True)
NR, NC, X0, Y0, C = 4362, 2440, 59799.282445142424, 49254.0, 10.0
SAT = {1: 0.403, 2: 0.439, 3: 0.430, 4: 0.520, 5: 0.614, 6: 0.766, 7: 0.472}   # H-TESSEL saturation by soil type
a = netCDF4.Dataset(E + "era5land_lux_20210712_13.nc"); t = netCDF4.num2date(a["valid_time"][:], a["valid_time"].units)
k = [i for i, d in enumerate(t) if (d.day, d.hour) == (13, 18)][0]
lat, lon = a["latitude"][:], a["longitude"][:]
sm = (7 * a["swvl1"][k] + 21 * a["swvl2"][k]) / 28
st = netCDF4.Dataset(E + "slt.nc"); sl = np.asarray(st["slt"][0]).round().astype(int)          # soil type grid differs: nearest cell
slt = sl[np.abs(np.asarray(lat)[:, None] - st["latitude"][:]).argmin(1)][:, np.abs(np.asarray(lon)[:, None] - st["longitude"][:]).argmin(1)]
Srel = np.clip(sm / np.vectorize(lambda c: SAT.get(c, np.nan))(slt), 0, 1)
ga = dict(np.load(S + "green_ampt_10m_v3b.npz"))
b = np.load("/ec/res4/scratch/lux0804/radar_scale_2026/basin_mask.npy").astype(bool)
x = X0 + (np.arange(NC) + 0.5) * C; y = Y0 + NR * C - (np.arange(NR) + 0.5) * C
X, Y = np.meshgrid(x, y); lo, la = Transformer.from_crs("EPSG:2169", "EPSG:4326", always_xy=True).transform(X, Y)
ii = np.abs(la[..., None] - lat).argmin(-1); jj = np.abs(lo[..., None] - lon).argmin(-1)
Sg = np.asarray(Srel)[ii, jj]; act = (ga["Ks"] > 0) & np.isfinite(Sg)
print("ERA5-Land %s; S basin mean %.3f (range %.3f-%.3f); theta_i v3b mean %.3f" % (t[k], Sg[b].mean(), Sg[b].min(), Sg[b].max(), ga["theta_i"][b & act].mean()))
ga["theta_i"][act] = np.minimum(ga["porosity"][act], Sg[act] * ga["porosity"][act])
ga["dtheta"] = np.where(ga["Ks"] > 0, np.clip(ga["porosity"] - ga["theta_i"], 0, None), 0).astype("float32")
print("theta_i sm mean %.3f; mean dtheta %.3f" % (ga["theta_i"][b & act].mean(), ga["dtheta"][b].mean()))
np.savez_compressed(S + "green_ampt_10m_v3b_sm.npz", **ga)
with open("/ec/res4/scratch/lux0804/radar_scale_2026/inputs/Alzette_sub_basin_10m_bridge_burn.asc") as f: hdr = "".join(f.readline() for _ in range(6))
for name in ("ga_ks_mmh", "ga_psi_mm", "ga_porosity"): os.system(f"cp /ec/res4/scratch/lux0804/lf84_cases/inputs_v3b/{name}.asc {O}/")
np.savetxt(O + "/ga_theta_init.asc", ga["theta_i"], fmt="%.4f", header=hdr.rstrip("\n"), comments="")
print("written", O)
