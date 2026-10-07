def T(x): return (3*x+1)//2 if x & 1 else x//2
def traj(x, cap=100000):
    seen={}
    s=0
    while x not in seen and x!=1:
        seen[x]=s; x=T(x); s+=1
    seen[1]=s
    return seen
def tst(x):
    s=0
    while x!=1:
        x = 3*x+1 if x&1 else x//2; s+=1
    return s
for c in (1,2):
    print(f"=== c={c}: coalescence of 3^A+c with 3^(A+1)+c ===")
    for A in range(4,32):
        x,y=3**A+c,3**(A+1)+c
        sx=traj(x); sy=traj(y)
        common=[v for v in sx if v in sy]
        v=min(common, key=lambda v: (sx[v]+sy[v]))
        print(f"A={A:3d}: meet at {v} (steps: {sx[v]} from 3^A+c, {sy[v]} from 3^(A+1)+c) | tst {tst(x)},{tst(y)}")
