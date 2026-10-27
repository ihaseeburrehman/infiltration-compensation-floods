"""Green-Ampt v2: independent wet-area mask from the soil map (no calibration).
Ks = 0 and theta_init = porosity on alluvium/colluvium, strongly gleyed units and seepage zones,
i.e. areas with a shallow water table where infiltrated water is not a loss (Fenicia et al. 2005)."""
import numpy as np, geopandas as gpd, os, pandas as pd
from rasterio.features import rasterize
from affine import Affine
S = "soil/"; O = "/ec/res4/scratch/lux0804/lf84_cases/inputs_v2"; os.makedirs(O, exist_ok=True)
NR, NC, X0, Y0, C = 4362, 2440, 59799.282445142424, 49254.0, 10.0
tf = Affine(C, 0, X0, 0, -C, Y0 + NR * C)
g = gpd.read_file(S + "so.soil-100000.gml").to_crs(2169)
l = g.soilBodyLabel.str.lower()
wet = l.str.contains("alluvions") | l.str.contains("suintement") | l.str.contains("fortement")
print("wet units:"); [print("  ", s) for s in sorted(g.soilBodyLabel[wet].unique())]
m = rasterize(((geom, 1) for geom in g.geometry[wet]), out_shape=(NR, NC), transform=tf, fill=0, dtype="uint8").astype(bool)
ga = dict(np.load(S + "green_ampt_10m.npz"))
b = np.load("/ec/res4/scratch/lux0804/radar_scale_2026/basin_mask.npy").astype(bool)
K0 = ga["Ks"][b].mean()
ga["Ks"][m] = 0.0; ga["theta_i"][m] = ga["porosity"][m]; ga["dtheta"][m] = 0.0
print("wet mask %.1f%% of basin; basin-mean Ks %.2f -> %.2f mm/h" % (100 * m[b].mean(), K0 * 3.6e6, ga["Ks"][b].mean() * 3.6e6))
np.savez_compressed(S + "green_ampt_10m_v2.npz", wetmask=m, **ga)
with open("/ec/res4/scratch/lux0804/radar_scale_2026/inputs/Alzette_sub_basin_10m_bridge_burn.asc") as f: hdr = "".join(f.readline() for _ in range(6))
for name, arr, fmt in [("ga_ks_mmh", ga["Ks"] * 3.6e6, "%.5g"), ("ga_psi_mm", ga["psi"] * 1000, "%.5g"),
                       ("ga_porosity", ga["porosity"], "%.4f"), ("ga_theta_init", ga["theta_i"], "%.4f")]:
    np.savetxt(O + f"/{name}.asc", arr, fmt=fmt, header=hdr.rstrip("\n"), comments="")
print("written", O)
