import sys,numpy as np
def load(f):
    d={}
    for l in open(f):
        k,v=l.strip().rsplit(",",1); d[k]=float("nan" if v=="NA" else v)
    return d
ref=load("ref_python.csv"); worst={}
for other in ("ref_r.csv","ref_octave.csv"):
    o=load(other); assert o.keys()==ref.keys(),(other,set(o)^set(ref))
    mx=0; wk=None
    for k in ref:
        a,b=ref[k],o[k]
        if np.isnan(a) and np.isnan(b): continue
        r=abs(a-b)/max(abs(a),1e-300) if a!=0 else abs(b)
        if r>mx: mx,wk=r,k
    print(f"python vs {other[4:-4]}: {len(ref)} numbers, max rel diff {mx:.2e} at {wk}")
