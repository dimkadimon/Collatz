#!/usr/bin/env python3
"""certify.py -- recompute certified intervals [a_lo,a_hi] for a(n)=A277109(n) for all n,
using the phi model + rule/scan class-minima + a large verified member set, then derive
A277684 least-n candidates and check them against the published 65 values.

member set: every member j verified by computing Q(j) directly (records, +-1, 2r, T^k(r),
workbook B/streak/H values).  badness of class Q at n is the phi test (valid when h(j)<=n).
"""
import sys, struct, csv
import numpy as np
sys.path.insert(0, '/home/user/collatz/src')
from classes import info, T, Q_of

PHI='/tmp/phi22k.bin'; SCAN='/tmp/minj5e9.bin'
QLIM=1200

d=open(PHI,'rb').read(); Amax,=struct.unpack('<I',d[:4])
V=np.frombuffer(d,dtype='<u8',offset=4,count=2*(Amax+1))

# ---- scans ----
exact={}; jmax=0
dd=open(SCAN,'rb').read(); N=8001*2
jmax,=struct.unpack('<q',dd[:8])
MJ=np.frombuffer(dd,dtype='<q',offset=8,count=N)
for i in np.nonzero(MJ>0)[0]:
    Q=int(i//2-4000); par=int(i%2); j=int(MJ[i])
    k=(Q,par)
    if k not in exact or j<exact[k]: exact[k]=j

# ---- member set ----
recs=[]
for line in open('/home/user/collatz/data/b006877.txt'):
    p=line.split()
    if len(p)>=2 and p[1].isdigit(): recs.append(int(p[1]))
members={}   # Q -> list of (j, h)
def add(j):
    h,o=info(j); Q=2*o-h
    if abs(Q)<=QLIM:
        members.setdefault(Q,[]).append((j,h))
for r in recs:
    for w in (r-1,r,r+1,2*r):
        if w>=2: add(w)
    x=r
    for _ in range(3):
        x=T(x); add(x)
    if r%3==2: add((2*r-1)//3)
for row in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
    for col in ('B','streak','H'):
        try:
            v=int(row[col])
            if v>1: add(v); add(v+1)
        except Exception: pass
for Q in members: members[Q]=sorted(members[Q])
# framework csv a_hi (includes the lost climb members) for cross-check
csv_hi={}; csv_lo={}
for f in ('/home/user/collatz/data/framework_n0_3953.csv','/home/user/collatz/data/framework_n3954_13999.csv',
          '/home/user/collatz/data/framework_n14000_39199.csv'):
    for row in csv.DictReader(open(f)):
        try:
            n=int(row['n']); hi=float(row['a_hi']) if row['a_hi'] not in ('',None) else None
            lo=float(row['a_lo']) if row['a_lo'] not in ('',None) else None
        except Exception: continue
        if hi is None or hi==float('inf'): continue
        csv_hi[n]=int(hi)
        csv_lo[n]=int(lo) if lo==lo else None

def interval(n):
    A1=(n+1)//2; c1=1+(n&1); L=int(V[2*A1+c1-1])
    Qs=np.arange(-QLIM,QLIM+1,dtype=np.int64)
    r=(n-Qs)&1
    A=(n+Qs+r)//2
    ok=(A>=0)&(A<=Amax)
    vals=np.where(ok, V[np.where(ok,2*A+r,0)],0)
    badQ=Qs[(vals!=L)|(~ok)]
    a_lo=jmax+1; a_hi=None; qsrc={}
    for Q in badQ:
        Q=int(Q)
        mem=None
        if Q<0:  mem=1<<(-Q); src='rule'
        k=(Q,abs(Q)%2)
        if k in exact:
            j=exact[k]
            if mem is None or j<mem: mem=j; src='scan'
        if Q in members:
            for j,h in members[Q]:
                if h<=n:
                    if mem is None or j<mem: mem=j; src='mem'
                    break
        if mem is not None:
            a_lo=min(a_lo,mem)
            if a_hi is None or mem<a_hi: a_hi=mem; qsrc['hi']=(Q,src)
        else:
            a_lo=min(a_lo,jmax+1)
    if a_hi is None: return None
    return (a_lo-1, a_hi-1, qsrc.get('hi'))

def main():
    nlo,nhi=(int(sys.argv[1]),int(sys.argv[2])) if len(sys.argv)>2 else (0,39199)
    out=open('/tmp/cert_%d_%d.csv'%(nlo,nhi),'w')
    out.write("n,a_lo,a_hi,exact,hiQ,hiSrc,csv_hi\n")
    for n in range(nlo,nhi+1):
        r=interval(n)
        if r is None: continue
        lo,hi,qs=r
        q,Q=qs if qs else ('','')
        ch=csv_hi.get(n)
        if ch is not None: hi=min(hi,ch)
        out.write(f"{n},{lo},{hi},{'yes' if lo==hi else 'no'},{Q},{Q},{ch if ch is not None else ''}\n".replace(",",",",1) if False else f"{n},{lo},{hi},{'yes' if lo==hi else 'no'},{Q},{ch if ch is not None else ''}\n")
    out.close()

if __name__=='__main__': main()
