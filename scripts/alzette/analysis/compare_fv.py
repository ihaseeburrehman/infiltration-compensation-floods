import numpy as np, os
S="/ec/res4/scratch/lux0804/"
G=["Pfaffenthal","Steinsel","Hesperange","Livange","Walferdange"]; OBS=dict(Walferdange=float("nan"),Pfaffenthal=134.5,Steinsel=131.5,Hesperange=122.6,Livange=98.75)
C=["R_1km_5min","F_graphcast","F_aifs","VM_graphcast","VM_aifs"]
SET=[("LF noinf",S+"radar_scale_2026/runs/L84_{c}_noinf_dp_1g/segflux_fv.npz"),("LF inf",S+"radar_scale_2026/runs/L84_{c}_ga3b_dp/segflux_fv.npz"),
     ("HP noinf",S+"hipims/runs/{c}_mg4/segflux_fv.npz"),("HP inf",S+"hipims/runs/{c}_inf_mg4/segflux_fv.npz")]
for g in G:
    print(f"\n{g} (obs {OBS[g]})"); print(f"{"case":13s}"+"".join(f"{s:>22s}" for s,_ in SET))
    for c in C:
        row=f"{c:13s}"
        for s,p in SET:
            p=p.format(c=c)
            if not os.path.exists(p): row+=f"{"-":>22s}"; continue
            z=np.load(p); q=(z[g+"_qn"]*z[g+"_L"]).sum(1); i=q.argmax()
            row+=f"{q[i]:7.1f} ({100*q[i]/OBS[g]:3.0f}%) h{z["steps"][i]:2d}".rjust(22)
        print(row)
