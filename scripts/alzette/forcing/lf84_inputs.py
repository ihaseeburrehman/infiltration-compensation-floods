"""LISFLOOD-FP 8.4 inputs identical to the HiPIMS cases: 50 m block-mean forecast rain on the 10 m grid,
and Green-Ampt grids (Ks mm/h, psi mm, porosity, theta_init) from the same soil file."""
import numpy as np, netCDF4, os
W = "/ec/res4/scratch/lux0804/radar_scale_2026"; H = "/ec/res4/scratch/lux0804/hipims"; O = "/ec/res4/scratch/lux0804/lf84_cases/inputs"
os.makedirs(O, exist_ok=True)
hdr = "".join(open(W + "/inputs/Alzette_sub_basin_10m_bridge_burn.asc").readline() for _ in range(1))
with open(W + "/inputs/Alzette_sub_basin_10m_bridge_burn.asc") as f: hdr = "".join(f.readline() for _ in range(6))
NR, NC = 4362, 2440
dem_nd = None
ga = np.load(H + "/soil/green_ampt_10m.npz")
for name, arr, fmt in [("ga_ks_mmh", ga["Ks"] * 3.6e6, "%.5g"), ("ga_psi_mm", ga["psi"] * 1000, "%.5g"),
                       ("ga_porosity", ga["porosity"], "%.4f"), ("ga_theta_init", ga["theta_i"], "%.4f")]:
    np.savetxt(O + f"/{name}.asc", arr, fmt=fmt, header=hdr.rstrip("\n"), comments="")
    print(name, "written, mean %.4g" % arr.mean(), flush=True)
act = np.loadtxt  # placeholder to keep namespace simple
B = 5; bi, bj = np.arange(NR) // B, np.arange(NC) // B
blk = bi[:, None] * (NC // B + 1) + bj[None, :]
import pandas as pd
dem = pd.read_csv(W + "/inputs/Alzette_sub_basin_10m_bridge_burn.asc", skiprows=6, sep=r"\s+", header=None, dtype=np.float32, engine="c").values
act = dem > -9000
u, inv = np.unique(blk[act], return_inverse=True); cnt = np.bincount(inv)
src = {"F_graphcast": "inputs_18utc/rain_graphcast_v2_96h.nc", "F_aifs": "inputs_18utc/rain_aifs_96h.nc",
       "VM_graphcast": "forcing/rain_VM_graphcast.nc", "VM_aifs": "forcing/rain_VM_aifs.nc"}
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
