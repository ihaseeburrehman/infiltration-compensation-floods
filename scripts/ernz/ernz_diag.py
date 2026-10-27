import numpy as np, rasterio, geopandas as gpd, pandas as pd, netCDF4, json, csv
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, ListedColormap
from scipy import ndimage
E='/Users/haseeb.rehman/Documents/Phd_thesis/Research_papers/WRF_vs_AI_LISFLOOD_v1/ernz2018/'
W='ernz_wd/'
with rasterio.open(E+'domain/dem10_cond.tif') as r: dem=r.read(1).astype('float64'); tf=r.transform
with rasterio.open(E+'domain/dem10_mean.tif') as r: dmean=r.read(1).astype('float64')
NR,NC=dem.shape; X0,YT,C=tf.c,tf.f,tf.a; act=dem>-9000
burnt=act&(dem<dmean-1e-3)
M=np.load(E+'domain/masks.npz'); lar=M['larochette'].astype(bool)
def rd(f): return pd.read_csv(W+f,sep=r'\s+',header=None,skiprows=6,dtype='float32').values
G9,G12,N12,GM=rd('ga_0009.wd'),rd('ga_0012.wd'),rd('noinf_0012.wd'),rd('ga.max')
for a in (G9,G12,N12,GM): a[~act]=np.nan
print('loaded', G12.shape)
T=json.load(open(E+'scripts/transects_ernz.json'))['Larochette']; gx,gy=T['gauge_x'],T['gauge_y']
rows=np.where(lar.any(1))[0]; cols=np.where(lar.any(0))[0]
r0,r1,c0,c1=rows.min()-20,rows.max()+20,cols.min()-20,cols.max()+20
ext=[X0+c0*C,X0+c1*C,YT-r1*C,YT-r0*C]
lines=gpd.read_file(E+'domain/lines/col_749_23.geojson').to_crs(2169)
stem=gpd.read_file(E+'domain/lines/col_2091.geojson').to_crs(2169)
eb=stem[stem.Name_2021.str.contains('blanche',case=False)]
ls=LightSource(315,45); sub=dem[r0:r1,c0:c1].copy(); sub[sub<-9000]=np.nan
hs=ls.hillshade(np.nan_to_num(sub,nan=np.nanmin(sub)),vert_exag=2,dx=C,dy=C)
def panel(ax,wd,title):
    ax.imshow(hs,cmap='gray',extent=ext,vmin=0,vmax=1.2,alpha=0.9)
    d=wd[r0:r1,c0:c1].copy(); d[d<0.05]=np.nan
    im=ax.imshow(d,extent=ext,cmap='Blues',vmin=0.05,vmax=2,interpolation='nearest')
    lines.plot(ax=ax,color='k',lw=0.2,alpha=0.5)
    ax.contour(np.linspace(ext[0]+C/2,ext[1]-C/2,c1-c0),np.linspace(ext[3]-C/2,ext[2]+C/2,r1-r0),lar[r0:r1,c0:c1].astype(float),[0.5],colors='r',linewidths=0.8)
    ax.plot(gx,gy,'*',color='orange',ms=12,mec='k'); ax.set_xlim(ext[0],ext[1]); ax.set_ylim(ext[2],ext[3]); ax.set_title(title,fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    v=np.nansum(np.where(lar,wd,0))*C*C/1e6; ax.text(0.02,0.02,f'water in catchment {v:.2f} Mm3',transform=ax.transAxes,fontsize=8,bbox=dict(fc='w',alpha=.8))
    return im
fig,axs=plt.subplots(2,2,figsize=(12,13))
panel(axs[0,0],G9,'Green-Ampt, hour 9 (03 UTC, observed peak)')
panel(axs[0,1],G12,'Green-Ampt, hour 12 (06 UTC, simulated peak)')
panel(axs[1,0],N12,'no infiltration, hour 12')
im=panel(axs[1,1],GM,'Green-Ampt, maximum depth 48 h')
fig.colorbar(im,ax=axs,shrink=0.5,label='depth (m), cells > 0.05 m'); fig.suptitle('Larochette catchment (red), radar 5 min forcing; star = gauge',fontsize=12)
fig.savefig('ernz_depth_catchment.png',dpi=130,bbox_inches='tight'); plt.close(fig)
# zoom around gauge 2.5 km
fig,axs=plt.subplots(1,2,figsize=(14,7))
for ax,wd,t in zip(axs,(G12,GM),('Green-Ampt hour 12','Green-Ampt max depth')):
    R=125; rr0=int((YT-gy)/C)-R; cc0=int((gx-X0)/C)-R; e2=[X0+cc0*C,X0+(cc0+2*R)*C,YT-(rr0+2*R)*C,YT-rr0*C]
    s=dem[rr0:rr0+2*R,cc0:cc0+2*R].copy(); s[s<-9000]=np.nan
    ax.imshow(ls.hillshade(np.nan_to_num(s,nan=np.nanmin(s)),vert_exag=2,dx=C,dy=C),cmap='gray',extent=e2,vmin=0,vmax=1.2)
    cs=ax.contour(np.linspace(e2[0]+5,e2[1]-5,2*R),np.linspace(e2[3]-5,e2[2]+5,2*R),s,levels=np.arange(200,400,10),colors='k',linewidths=0.3,alpha=0.6); ax.clabel(cs,fontsize=6,fmt='%d')
    d=wd[rr0:rr0+2*R,cc0:cc0+2*R].copy(); d[d<0.05]=np.nan
    im=ax.imshow(d,extent=e2,cmap='Blues',vmin=0.05,vmax=3,interpolation='nearest')
    lines.plot(ax=ax,color='cyan',lw=0.5); eb.plot(ax=ax,color='magenta',lw=0.8,ls='--')
    # section line
    s_=np.linspace(-300,300,2); ax.plot(T['x']-T['uy']*s_,T['y']+T['ux']*s_,'orange',lw=1.5)
    ax.plot(gx,gy,'*',color='orange',ms=14,mec='k'); ax.set_xlim(e2[0],e2[1]); ax.set_ylim(e2[2],e2[3]); ax.set_title(t+' - Larochette ±1.25 km'); ax.set_aspect('equal')
fig.colorbar(im,ax=axs,shrink=0.6,label='depth (m)'); fig.savefig('ernz_depth_larochette_zoom.png',dpi=130,bbox_inches='tight'); plt.close(fig)
# profile along Ernz blanche stem from upstream end to gauge, then beyond
geom=eb.geometry.unary_union
from shapely.ops import linemerge; from shapely.geometry import LineString, Point
if geom.geom_type!='LineString': geom=linemerge(geom)
if geom.geom_type!='LineString': geom=max(geom.geoms,key=lambda g:g.length)
# orient: upstream end = higher DEM
def samp(x,y,arr,win=1):
    i=int((YT-y)/C); j=int((x-X0)/C); w=arr[max(i-win,0):i+win+1,max(j-win,0):j+win+1]; return np.nanmin(w) if np.isfinite(w).any() else np.nan
p0,p1=geom.coords[0],geom.coords[-1]
if samp(*p0,dem)<samp(*p1,dem): geom=LineString(geom.coords[::-1])
ds=np.arange(0,geom.length,20); pts=[geom.interpolate(d) for d in ds]
bed=np.array([samp(p.x,p.y,dem) for p in pts])
# channel cell: argmin dem in window; take wd at that cell
def chan(p,wd):
    i=int((YT-p.y)/C); j=int((p.x-X0)/C); w=dem[i-1:i+2,j-1:j+2]; w=np.where(w<-9000,np.nan,w)
    if not np.isfinite(w).any(): return np.nan
    k=np.unravel_index(np.nanargmin(w),w.shape); return wd[i-1+k[0],j-1+k[1]]
wse12=bed+np.array([chan(p,G12) for p in pts]); wse9=bed+np.array([chan(p,G9) for p in pts]); dmax=np.array([chan(p,GM) for p in pts])
dg=np.argmin([(p.x-gx)**2+(p.y-gy)**2 for p in pts]); print('gauge at chainage',ds[dg],'m of',geom.length)
fig,ax=plt.subplots(2,1,figsize=(14,8),sharex=True)
ax[0].plot(ds/1000,bed,'k',lw=1,label='conditioned bed (10 m)'); ax[0].plot(ds/1000,wse9,'c',lw=0.8,label='water surface GA h9'); ax[0].plot(ds/1000,wse12,'b',lw=0.8,label='water surface GA h12')
ax[0].axvline(ds[dg]/1000,color='orange',label='Larochette gauge'); ax[0].legend(); ax[0].set_ylabel('m a.s.l.'); ax[0].set_title('Ernz Blanche main stem profile (WFD line col_2091)')
ok=np.isfinite(bed); rise=np.diff(bed); rise[~np.isfinite(rise)]=0; bad=np.where(rise>0.5)[0]; print('bed rises >0.5 m downstream (20 m step):',len(bad)); 
for b in bad[:30]: print(f'  chainage {ds[b]/1000:.2f} km rise {rise[b]:.2f} m  xy {pts[b].x:.0f},{pts[b].y:.0f}')
ax[1].plot(ds/1000,dmax,'r',lw=0.8,label='max depth in channel cell'); ax[1].plot(ds/1000,wse12-bed,'b',lw=0.8,label='depth h12'); ax[1].axvline(ds[dg]/1000,color='orange'); ax[1].legend(); ax[1].set_ylabel('m'); ax[1].set_xlabel('chainage from upstream end (km)')
fig.savefig('ernz_blanche_profile.png',dpi=130,bbox_inches='tight'); plt.close(fig)
# ponding check: connected water bodies > 0.3 m not touching burnt cells within catchment
for name,wd in (('GA h12',G12),('noinf h12',N12),('GA max',GM)):
    deep=(np.nan_to_num(wd)>0.3)&lar
    near=ndimage.binary_dilation(burnt,iterations=5)
    lab,n=ndimage.label(deep); touch=np.unique(lab[deep&near]); 
    ids=np.setdiff1d(np.arange(1,n+1),touch); areas=ndimage.sum(deep,lab,ids)*C*C/1e4; vols=ndimage.sum(np.nan_to_num(wd),lab,ids)*C*C
    o=np.argsort(areas)[::-1][:8]; print(name,'off-channel bodies >0.3 m:',len(ids),'total ha %.1f, volume %.3f Mm3'%(areas.sum(),vols.sum()/1e6), ' total catchment water %.3f Mm3'%(np.nansum(np.where(lar,wd,0))*C*C/1e6))
    for k in o:
        cy,cx=ndimage.center_of_mass(deep,lab,ids[k]); print(f'   {areas[k]:.1f} ha, {vols[k]/1e3:.0f} dam3, maxd {ndimage.maximum(np.nan_to_num(wd),lab,ids[k]):.2f} m at {X0+cx*C:.0f},{YT-cy*C:.0f}')
# rain over Larochette catchment per hour vs hydrographs
d=netCDF4.Dataset(E+'forcing/rain_R_1km_5min.nc'); rain=d['rainfall_depth']; n=rain.shape[0]; hr=np.zeros(48)
for t in range(min(n,576)):
    a=rain[t][lar]; hr[t//12]+=np.nanmean(a)
print('Larochette catchment hourly rain mm:', np.round(hr[:20],1), 'total', hr.sum().round(1))
obs=pd.read_csv(E+'results/obs_larochette_utc.csv'); sn=pd.read_csv(E+'results/hydro_R_1km_5min_noinf.csv'); sg=pd.read_csv(E+'results/hydro_R_1km_5min_ga.csv')
fig,ax=plt.subplots(figsize=(10,5)); ax2=ax.twinx(); ax2.bar(np.arange(48)+0.5,hr,color='grey',alpha=.4,label='catchment rain (mm/h)'); ax2.invert_yaxis(); ax2.set_ylim(40,0)
ax.plot(obs.h,obs.Q,'k',lw=2,label='observed'); ax.plot(sn.hour,sn.Larochette,'r',label='no infiltration'); ax.plot(sg.hour,sg.Larochette,'b',label='Green-Ampt'); ax.set_xlim(0,30); ax.set_xlabel('hour from 31 May 18 UTC'); ax.set_ylabel('m3/s'); ax.legend(loc='upper right'); ax2.set_ylabel('mm/h')
fig.savefig('ernz_larochette_timing.png',dpi=130,bbox_inches='tight')
np.savez_compressed('ernz_diag_series.npz',rain=hr,bed=bed,ds=ds,wse12=wse12)
