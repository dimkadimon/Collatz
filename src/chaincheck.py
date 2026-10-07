#!/usr/bin/env python3
"""chaincheck.py -- decide first-break j for tst(2^n+j)=tst(2^n+1) (class method).

class(j) = (Q,p), Q=2o-h, p=h mod 2 (accelerated steps).  tst(2^n+j) is constant on
classes: tst(2^n+j) = n+o+tst(3^A+c), A=(n+Q+r)/2, c=1+r, r=(n-p)&1.
Class is good iff phi(A,c) == phi(A1,c1).  Smallest bad j = min over bad classes of
their minimum member: <=5e9 exact from the minj scan; larger ones from a candidate pool
(A006877 records, Frinkle breaks/streaks, 2^k, m*2^k), each verified by direct height
computation.  kind=scan -> EXACT answer; kind=pool -> upper bound (verified breaker).
"""
import struct, sys, csv
from collections import defaultdict

T=lambda x:(3*x+1)//2 if x&1 else x//2
def tst(x):
    s=0
    while x!=1:
        x=3*x+1 if x&1 else x//2; s+=1
    return s
_ic={}
def info(j):
    r=_ic.get(j)
    if r is not None: return r
    x=j;h=0;o=0
    while x!=1:
        if x&1:o+=1
        x=T(x);h+=1
    _ic[j]=(h,o); return h,o

def load_phi(path):
    phi={}
    with open(path,'rb') as f:
        while True:
            b=f.read(13)
            if len(b)<13: break
            A,c,p=struct.unpack('<ibq',b); phi[(A,c)]=p
    return phi
def load_minj(path):
    with open(path,'rb') as f:
        JMAX=struct.unpack('<Q',f.read(8))[0]
        N=8001*2
        MJ=struct.unpack('<%dq'%N,f.read(8*N))
        MJO=struct.unpack('<%dq'%N,f.read(8*N))
        maxQ,minQ=struct.unpack('<qq',f.read(16))
    return JMAX,MJ,MJO,maxQ,minQ

def build_pool():
    pool=set()
    for line in open('/home/user/collatz/data/b006877.txt'):
        p=line.split()
        if len(p)>=2 and p[1].isdigit(): pool.add(int(p[1]))
    for r in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
        for col in ('streak','B'):
            try:
                v=int(r[col])
                if v>0: pool.add(v)
                if col=='streak' and v>0: pool.add(v+1)
            except (ValueError,TypeError): pass
    for k in range(0,131):
        pool.add(1<<k)
        if k<=120:
            for m in range(3,64,2): pool.add(m<<k)
    return pool

def main():
    phi=load_phi(sys.argv[1]); JMAX,MJ,MJO,maxQ,minQ=load_minj(sys.argv[2])
    nlo,nhi=int(sys.argv[3]),int(sys.argv[4])
    out=open(sys.argv[5],'w') if len(sys.argv)>5 else sys.stdout
    pool=build_pool()
    sig={}
    for j in pool:
        h,o=info(j); sig[j]=(2*o-h,h&1,h,o)
    byclass=defaultdict(list)
    for j,(Q,p,h,o) in sig.items(): byclass[(Q,p)].append(j)
    for k in byclass: byclass[k].sort()
    scanned=[]                      # (Q,p,j) with exact class minimum j <= JMAX
    for i in range(len(MJ)):
        j=MJ[i]
        if j>0:
            Q=i//2-4000; p=i%2
            scanned.append((Q,p,j))
    print(f"pool classes={len(byclass)} scanned classes={len(scanned)}",file=sys.stderr)
    out.write("n,P,aP,first_bad_j,a_n,kind,verified,note\n")
    for n in range(nlo,nhi+1):
        A1=(n+1)//2; c1=1+(n&1); L=phi.get((A1,c1))
        if L is None:
            out.write(f"{n},,,NA,NA,phi-missing,,\n"); continue
        gcache={}
        def isgood(Q,p):
            v=gcache.get((Q,p))
            if v is None:
                r=(n-p)&1; A=(n+Q+r)//2
                if A<0: v=False
                else:
                    t=phi.get((A,1+r))
                    v = None if t is None else (t==L)
                gcache[(Q,p)]=v
            return v
        best=None; kind=None
        for (Q,p,j) in scanned:
            if isgood(Q,p) is False and (best is None or j<best):
                best=j; kind='scan'
        if best is None or kind=='pool':
            for (Q,p),lst in byclass.items():
                cand=[v for v in lst if sig[v][2]<=n]
                if not cand: continue
                if isgood(Q,p) is not False: continue
                j=cand[0]
                if best is None or j<best: best=j; kind='pool'
        if best is None:
            out.write(f"{n},,,-,none,,\n"); continue
        h,o=info(best)
        # direct verification: tst(2^n+best) must differ from tst(2^n+1)
        def Nplus(j):
            x=j;o2=0
            for _ in range(n):
                if x&1:o2+=1
                x=T(x)
            return n+o2+tst(3**o2+x)
        t1=Nplus(1); t2=Nplus(best)
        ok=(t2!=t1)
        out.write(f"{n},{2*o-h},{h&1},{best},{best-1},{kind},{'yes' if ok else 'NO'},"
                  f"{'' if ok else 'VERIFY-FAILED'}\n")
        out.flush()
    out.close()

if __name__=='__main__': main()
