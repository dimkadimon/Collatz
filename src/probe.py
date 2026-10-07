import sys
from functools import lru_cache

def T(x):  # accelerated step
    return (3*x+1)//2 if x & 1 else x//2

def acc_info(j):
    """returns (h, o, traj) where h = #T-steps to reach 1, o = #odd values among T^0..T^{h-1}"""
    x=j; h=0; o=0; traj=[]
    while x!=1:
        traj.append(x)
        if x&1: o+=1
        x=T(x); h+=1
    return h,o

def tst(x):
    s=0
    while x!=1:
        x = 3*x+1 if x&1 else x//2
        s+=1
    return s

# 1) verify formula tst(2^n+j) = n + o_n(j) + tst(3^o_n(j) + T^n(j))
print("== formula check ==")
bad=0
for n in range(1,26):
    for j in range(1,40):
        x0=(1<<n)+j
        # real trajectory for n steps, then tst of remainder
        x=x0; steps=0; o=0
        for t in range(n):
            if x&1: o+=1
            x = (3*x+1)//2 if x&1 else x//2
        lhs = tst(x0)
        rhs = n + o + tst(x)
        if lhs!=rhs: bad+=1; print("MISMATCH",n,j,lhs,rhs)
print("formula check mismatches:", bad)

# 2) run structure: run at 2^n+1, and signature data
print("== run structure ==")
for n in list(range(6,26))+[30,40]:
    base=1<<n
    t0=tst(base+1)
    k=1
    heights={}
    while True:
        t=tst(base+k+1)
        if t!=t0: break
        k+=1
    print(f"n={n:3d} a(n)={k:8d}  tst(2^n+1)={t0}")
