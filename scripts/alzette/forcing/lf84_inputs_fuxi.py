"""FuXi (corrected latitude order) and FuXi rescaled to the radar total: 50 m block-mean rain for LISFLOOD-FP 8.4.
"""
import numpy as np, netCDF4, os
W = "/ec/res4/scratch/lux0804/radar_scale_2026"; H = "/ec/res4/scratch/lux0804/hipims"; O = "/ec/res4/scratch/lux0804/lf84_cases/inputs"
os.makedirs(O, exist_ok=True)
NR, NC = 4362, 2440
act = np.loadtxt  # placeholder to keep namespace simple
B = 5; bi, bj = np.arange(NR) // B, np.arange(NC) // B
blk = bi[:, None] * (NC // B + 1) + bj[None, :]
import pandas as pd
dem = pd.read_csv(W + "/inputs/Alzette_sub_basin_10m_bridge_burn.asc", skiprows=6, sep=r"\s+", header=None, dtype=np.float32, engine="c").values
act = dem > -9000
u, inv = np.unique(blk[act], return_inverse=True); cnt = np.bincount(inv)
src = {"F_fuxi": "inputs_18utc/rain_fuxi_noflip_96h.nc", "VM_fuxi": "forcing/rain_VM_fuxi2.nc"}   # noflip = corrected FuXi latitude order
for case, f in src.items():
    a = netCDF4.Dataset(W + "/" + f); out = netCDF4.Dataset(O + f"/rain_{case}_b50.nc", "w", format="NETCDF4")
    for d in a.dimensions: out.createDimension(d, len(a.dimensions[d]))
    for vn in ("time", "x", "y"):
        v = out.createVariable(vn, a[vn].dtype, a[vn].dimensions); v.setncatts({k: a[vn].getncattr(k) for k in a[vn].ncattrs()}); v[:] = a[vn][:]
    r = out.createVariable("rainfall_depth", "f4", ("time", "y", "x"), zlib=True, complevel=1, chunksizes=(1, NR, NC))
    r.setncatts({k: a["rainfall_depth"].getncattr(k) for k in a["rainfall_depth"].ncattrs() if k != "_FillValue"})
    tot = 0.0
    for i in range(a["rainfall_depth"].shape[0]):
        x = np.asarray(a["rainfall_depth"][i], dtype=np.float64)
        m = np.bincount(inv, x[act]) / cnt; g = np.zeros((NR, NC), np.float32); g[act] = m[inv]
        r[i] = g; tot += m[inv].mean()
    out.close(); print(case, "basin mean %.2f mm" % tot, flush=True)
