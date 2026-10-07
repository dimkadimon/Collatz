#!/usr/bin/env python3
"""an_fast.py -- compute A277109(n) = a(n) for a range of n, exactly.

Theory (class reduction):
  For j with h(j) := accelerated steps of j to 1, o(j) := #odd values before 1,
  Q := 2o-h, p := h mod 2, r := (n-p) mod 2, and q := T^n(j) (accelerated):
      tst(2^n + j) = n + o_n(j)[+extra] ... which for h(j) <= n collapses to
      tst(2^n + j) = (n + h - r)/2 + phi(A, 1+r),   A = (n + Q + r)/2,
  where phi(A,c) = A + tst(3^A + c)  (verified separately).
  Hence tst(2^n+j) depends only on Q (since p == Q mod 2), and
      a(n) = min{ minj(Q) : Q bad } - 1,
  bad Q  <=>  phi((n+Q+r)/2, 1+r) != phi((n+1)//2, 1+(n mod 2)).
  minj(Q) = exact smallest j in class Q, from the reverse-tree DP (dp_classes.py).
"""
import numpy as np, struct, sys, argparse

def load_phi(path):
    d=open(path,'rb').read()
    Amax,=struct.unpack('<I',d[:4])
    V=np.frombuffer(d,dtype='<u8',offset=4,count=2*(Amax+1))
    return int(Amax), V

def load_classmin(path, Qlim):
    dp={}
    for line in open(path):
        if line[0]=='#': continue
        f=line.split()
        Q=int(f[0]); M=int(f[2]); h=int(f[3])
        if abs(Q)<=Qlim:
            dp[Q]=[M,h]
    return dp

def T(x): return (3*x+1)//2 if x&1 else x//2
def tst_direct(x):
    s=0
    while x!=1:
        if x&1: x=3*x+1; s+=1
        x>>=1; s+=1
    return s
def tst_2n_plus_j(n,j):
    x=j;o=0
    for _ in range(n):
        if x&1: o+=1
        x=T(x)
    return n+o+tst_direct(3**o+x)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('phi'); ap.add_argument('classmin')
    ap.add_argument('nlo',type=int); ap.add_argument('nhi',type=int); ap.add_argument('out')
    ap.add_argument('--qlim',type=int,default=900)
    a=ap.parse_args()
    Amax,V=load_phi(a.phi)
    dp=load_classmin(a.classmin,a.qlim)
    Qs=np.arange(-a.qlim,a.qlim+1,dtype=np.int64)
    dpQ=np.array([dp.get(int(Q),[-1,10**9])[0] for Q in Qs],dtype=object)
    dph=np.array([dp.get(int(Q),[-1,10**9])[1] for Q in Qs],dtype=np.int64)
    # coverage threshold: smallest j with h(j) == Hlimit  (from the class-min DP extremes)
    Hmax=max(h for _,h in dp.values())
    tau=min([M for (M,h) in dp.values() if h==Hmax] or [10**30])
    f=open(a.out,'w'); f.write("n,a_n,first_bad_j,Q_bad,kind,exact,note\n")
    stats={}
    for n in range(a.nlo,a.nhi+1):
        A1=(n+1)//2; c1=1+(n&1); L=V[2*A1+c1-1]
        r=(n-Qs)&1
        A=(n+Qs+r)//2
        idx=2*A+r
        ok=(A>=0)&(A<=Amax)
        vals=np.where(ok, V[np.where(ok,idx,0)], 0)
        bad=(vals!=L)|(~ok)
        hok=(dph<=n)          # class minimum must be a "landed" j (h<=n)
        usable=bad&hok&(dpQ!=-1)
        note=''
        if usable.any():
            cand=[(dpQ[i],int(Qs[i])) for i in np.nonzero(usable)[0]]
            first_bad,Qbad=min(cand)
            exact = first_bad < tau
            kind='class'
        else:
            first_bad=None; Qbad=None; exact=False; kind='none'
        if not exact: note='coverage-limited'
        if exact and dph[int(np.nonzero(usable)[0][0])]>2000: note='h-margin-low'
        f.write(f"{n},{'' if first_bad is None else first_bad-1},{'' if first_bad is None else first_bad},"
                f"{'' if Qbad is None else Qbad},{kind},{'yes' if exact else 'no'},{note}\n")
        stats[kind]=stats.get(kind,0)+1
    f.close()
    print("tau(DP coverage) =",tau, "stats:",stats, file=sys.stderr)

if __name__=='__main__': main()
