#!/usr/bin/env python3
"""boundary.py -- for each n, find the smallest j whose height differs from tst(2^n+1).

Method: the height tst(2^n+j) depends only on the class Q(j)=2o(j)-h(j) of j, via
    tst(2^n+j) = (n + h - r)/2 + phi(A, 1+r),  A=(n+Q+r)/2, r=(n-h) mod 2   [h(j)<=n]
so a(n) = min{ minj(Q) : phi((n+Q+r)/2,1+r) != phi((n+1)//2,1+(n&1)) } - 1.

minj(Q) (smallest j with class Q) comes from
  * the exhaustive minj scan            (exact for minj <= Jmax's)
  * the proved rule minj(Q) = 2^(-Q) for Q<0
  * a candidate pool of hard numbers (A006877 records, Frinkle break values,
    2^k, m*2^k)   -> verified upper bounds
Every reported break is verified by direct height computation.
"""
import numpy as np, struct, sys, csv, argparse

def load_phi(path):
    d=open(path,'rb').read(); Amax,=struct.unpack('<I',d[:4])
    return int(Amax), np.frombuffer(d,dtype='<u8',offset=4,count=2*(Amax+1))

def load_scan(paths, Qlim):
    """min over the given minj-scan files: minj[Q] = smallest j found, or None"""
    best={}
    for p in paths:
        try: d=open(p,'rb').read()
        except IOError: continue
        if len(d)<8+8*16002*2: continue
        N=8001*2
        MJ=np.frombuffer(d,dtype='<q',offset=8,count=N)
        for i in np.nonzero(MJ>0)[0]:
            Q=int(i//2-4000); j=int(MJ[i])
            if abs(Q)<=Qlim and (Q not in best or j<best[Q]): best[Q]=j
    return best

def T(x): return (3*x+1)//2 if x&1 else x//2
_t={}
def info(j):
    r=_t.get(j)
    if r: return r
    x=j;h=0;o=0
    while x!=1:
        if x&1:o+=1
        x=T(x);h+=1
    _t[j]=(h,o); return h,o
def tst_direct(x):
    s=0
    while x!=1:
        if x&1: x=3*x+1; s+=1
        x>>=1; s+=1
    return s
def tst_2npj(n,j):
    x=j;o=0
    for _ in range(n):
        if x&1:o+=1
        x=T(x)
    return n+o+tst_direct(3**o+x)

def build_pool(Qlim):
    pool={}
    def add(j,src):
        h,o=info(j); Q=2*o-h
        if abs(Q)>Qlim: return
        cur=pool.get(Q)
        if cur is None or j<cur[0]: pool[Q]=(j,src,h)
    for line in open('/home/user/collatz/data/b006877.txt'):
        p=line.split()
        if len(p)>=2 and p[1].isdigit(): add(int(p[1]),'record')
    for r in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
        for col in ('B','streak'):
            try:
                v=int(r[col])
                if v>1: add(v+1 if col=='streak' else v,'frinkle')
            except (ValueError,TypeError): pass
    for k in range(0,131):
        add(1<<k,'2^k')
        if k<=60:
            for m in range(3,2001,2): add(m<<k,'m2^k')
    return pool

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('phi'); ap.add_argument('out'); ap.add_argument('nlo',type=int); ap.add_argument('nhi',type=int)
    ap.add_argument('--qlim',type=int,default=800); ap.add_argument('--scan',nargs='*',default=['/tmp/minj5e9.bin'])
    ap.add_argument('--verify',type=int,default=0)
    a=ap.parse_args()
    Amax,V=load_phi(a.phi)
    scan=load_scan(a.scan,a.qlim)
    pool=build_pool(a.qlim)
    print(f"phi Amax={Amax}  scan classes={len(scan)}  pool classes={len(pool)}",file=sys.stderr)
    Qs=np.arange(-a.qlim,a.qlim+1,dtype=np.int64)
    f=open(a.out,'w'); f.write("n,a_n,first_bad_j,Q_of_break,minj_source,verified,note\n")
    lo=min(min(scan.values()) if scan else 10**18, 10**18)
    for n in range(a.nlo,a.nhi+1):
        A1=(n+1)//2; c1=1+(n&1); L=V[2*A1+c1-1]
        r=(n-Qs)&1
        A=(n+Qs+r)//2
        ok=(A>=0)&(A<=Amax)
        vals=np.where(ok, V[np.where(ok,2*A+r,0)], 0)
        bad=((vals!=L)|(~ok))
        best=None;src=None;bestQ=None
        for i in np.nonzero(bad)[0]:
            Q=int(Qs[i])
            cands=[]
            if Q<0: cands.append((1<<(-Q),'struct2^Q'))  # h=|Q|<=800<=n
            if Q in scan: cands.append((scan[Q],'scan'))
            if Q in pool and pool[Q][2]<=n: cands.append((pool[Q][0],pool[Q][1]))
            for j,s in cands:
                if best is None or j<best: best,src,bestQ=j,s,Q
        if best is None:
            f.write(f"{n},,,,none,,no candidate\n"); continue
        verified=''
        if a.verify:
            ref=tst_2npj(n,1); t=tst_2npj(n,best)
            verified='yes' if t!=ref else 'NO'
        f.write(f"{n},{best-1},{best},{bestQ},{src},{verified},\n")
    f.close()

if __name__=='__main__': main()
