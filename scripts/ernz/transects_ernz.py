"""Discharge transects for the Ernz grid, same format as paper_v4 transects_fullvalley.json (x, y, W0, ux, uy, segs).
Main stems from the WFD surface water body lines (geoportail collection 2091, Name_2021: 'Ernz noire', 'Ernz blanche',
'Halerbaach'). Centre = conditioned channel cell nearest to the point on the stem; flow direction = local tangent of the stem,
oriented downhill; transect perpendicular, +-EXT m, cut into cell segments every 5 cm. Outlets sit 300 m above the downstream
end of the stem. Mullerthal = AGE gauge assumed at the Heringer Millen bridge (49.7896 N, 6.3055 E).
usage: python transects_ernz.py <domain dir> <col_2091 geojson> <out json>"""
import sys, json, numpy as np, geopandas as gpd, rasterio
from shapely.geometry import Point
from shapely.ops import linemerge, unary_union
from pyproj import Transformer
D, LINES, OUT = sys.argv[1:4]
with rasterio.open(f"{D}/dem10_cond.tif") as s: z = s.read(1); tf = s.transform
with rasterio.open(f"{D}/dem10_mean.tif") as s: zm = s.read(1)
NR, NC = z.shape; CS = tf.a; X0, YTOP = tf.c, tf.f
ch = (z > -9000) & (np.abs(zm - z) > 0.001); rr, cc = np.where(ch); cx = X0 + (cc + .5) * CS; cy = YTOP - (rr + .5) * CS
tr = Transformer.from_crs("EPSG:4326", "EPSG:2169", always_xy=True)
g = gpd.read_file(LINES).to_crs(2169)
def stem(name, near=None):
    """merged stem; if it stays multipart, the part nearest to `near` (a Point) or, for outlets, the part with the lowest end"""
    M = unary_union(g[g.Name_2021 == name].geometry); M = linemerge(M) if M.geom_type == "MultiLineString" else M
    if M.geom_type != "MultiLineString": return M
    if near is not None: return min(M.geoms, key=lambda l: l.distance(near))
    return min(M.geoms, key=lambda l: min(np.nan_to_num(zat(l.interpolate(300)), nan=1e9), np.nan_to_num(zat(l.interpolate(l.length - 300)), nan=1e9)))
def zat(p):
    r_, c_ = int((YTOP - p.y) / CS), int((p.x - X0) / CS)
    return zm[r_, c_] if 0 <= r_ < NR and 0 <= c_ < NC and zm[r_, c_] > -9000 else np.nan
# name: (lon, lat) or (x, y) of the gauge, half width, stem, 'outlet'
G = {"Larochette": ((83689, 95240), 300, "Ernz blanche", None), "Beaufort": ((91030, 98490), 150, "Halerbaach", None),
     "Mullerthal": (tr.transform(6.3055, 49.7896), 200, "Ernz noire", None), "Grundhof": (None, 200, "Ernz noire", "outlet"),
     "Reisdorf": (None, 250, "Ernz blanche", "outlet")}
T = {}
for n, (pt, ext, river, mode) in G.items():
    L = stem(river, None if mode == "outlet" else Point(*pt)); e0, e1 = Point(L.coords[0]), Point(L.coords[-1])
    z0, z1 = zat(L.interpolate(300)), zat(L.interpolate(L.length - 300))         # 300 m inside each end (ends may be outside the domain)
    downstream_end = e1 if np.nan_to_num(z1, nan=1e9) < np.nan_to_num(z0, nan=1e9) else e0
    if mode == "outlet":
        s_ = L.project(downstream_end); s_ = s_ - 300 if s_ > L.length / 2 else s_ + 300; pl = L.interpolate(s_); x, y = pl.x, pl.y
    else:
        x, y = pt; pl = L.interpolate(L.project(Point(x, y)))
    k = np.argmin(np.hypot(cx - pl.x, cy - pl.y)); x0, y0 = cx[k], cy[k]
    s_ = L.project(pl); a, b = L.interpolate(max(s_ - 60, 0)), L.interpolate(min(s_ + 60, L.length))
    ux, uy = b.x - a.x, b.y - a.y; nrm = np.hypot(ux, uy); ux, uy = ux / nrm, uy / nrm
    if L.project(downstream_end) < s_: ux, uy = -ux, -uy                      # point towards the downstream end
    ds = 0.05; sv = np.arange(-ext, ext, ds) + ds / 2; px, py = x0 - uy * sv, y0 + ux * sv
    r = np.floor((YTOP - py) / CS).astype(int); c = np.floor((px - X0) / CS).astype(int)
    cut = np.flatnonzero((np.diff(r) != 0) | (np.diff(c) != 0)) + 1
    segs = [[int(r[i]), int(c[i]), float(sv[i] - ds / 2), float(sv[j - 1] + ds / 2)] for i, j in zip(np.r_[0, cut], np.r_[cut, len(sv)])]
    segs = [sg for sg in segs if 0 <= sg[0] < NR and 0 <= sg[1] < NC and z[sg[0], sg[1]] > -9000]
    T[n] = dict(x=float(x0), y=float(y0), W0=10.0, ux=float(ux), uy=float(uy), gauge_x=float(x), gauge_y=float(y), stem=river, segs=segs)
    print(f"{n:11s} point ({x:.0f},{y:.0f}) -> channel cell ({x0:.0f},{y0:.0f}) {np.hypot(x0 - x, y0 - y):.0f} m away, {river}, dir ({ux:.2f},{uy:.2f}), bed {z[rr[k], cc[k]]:.1f} m, {len(segs)} segments")
json.dump(T, open(OUT, "w")); print("written", OUT)
