"""Reference values from the Python implementation (see compare.py).

Usage: python ref_python.py OUT.csv
The same test fields are built in ref_r.R and ref_octave.m (1-based indices there).
"""
import sys, numpy as np, geoareaweight as g
out=[]
def put(k,a):
    for i,v in enumerate(np.atleast_1d(np.asarray(a,float)).ravel()): out.append(f"{k}[{i}],{v:.15e}")
lat=np.arange(-90,90.1,2.5); lon=np.arange(0,360,2.5); nt=5
I=np.arange(1,lon.size+1)[None,None,:]; J=np.arange(1,lat.size+1)[None,:,None]; T=np.arange(1,nt+1)[:,None,None]
x=np.sin(0.3*I+0.17*J+1.1*T)+0.01*J
x[(I*7+J*3+T)%11==0]=np.nan
put("w_cos",g.area_weights(lat,method="cos")); put("w_band",g.area_weights(lat)); put("w_ell",g.area_weights(lat,method="ellipsoid"))
put("area_sum_sphere",g.cell_area(lat,lon).sum()); put("area_sum_ell",g.cell_area(lat,lon,ellipsoid=True).sum())
for m in ("cos","band","ellipsoid","none"): put("mean_"+m,g.area_mean(x,lat,lon,method=m))
mk=np.zeros((lat.size,lon.size),bool); mk[20:50,10:90]=True
put("mean_mask",g.area_mean(x,lat,lon,mask=mk)); put("mean_skipna_false",g.area_mean(x,lat,lon,skipna=False))
put("integral",g.area_integral(np.nan_to_num(x[0]),lat,lon)/1e12)
lat2=np.array([85,60,33,10,-12,-41,-70,-88.]); lon2=np.array([0,7,30,31,100,200,300.])
J2=np.arange(1,lat2.size+1)[:,None]; I2=np.arange(1,lon2.size+1)[None,:]
x2=np.cos(0.5*I2+0.9*J2)+0.1*J2
put("w2d_irreg",g.weights_2d(lat2,lon2)); put("mean_irreg_band",g.area_mean(x2,lat2,lon2)); put("mean_irreg_cos",g.area_mean(x2,lat2,lon2,method="cos"))
put("mean_irreg_ell",g.area_mean(x2,lat2,lon2,method="ellipsoid"))
cond=(x[0]>0.5); put("frac",g.area_fraction(np.where(np.isnan(x[0]),np.nan,cond.astype(float)),lat,lon))
open(sys.argv[1],"w").write("\n".join(out)+"\n")
