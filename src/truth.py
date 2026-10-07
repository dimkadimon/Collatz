#!/usr/bin/env python3
"""truth.py -- determine the true first break j for tst(2^n + j) = tst(2^n + 1),
using  phi-table  (from phi.c) + a candidate list of "hard" numbers + the
minj(Q,p) table (from minj.c) + powers of two + m*2^k families.

Prints, for each n: the structure of good/bad Q classes, the smallest bad j found
(exact if it comes from the scanned table, else an upper bound verified directly).
"""
import struct, sys, json

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
        n=8001*2
        MJ=struct.unpack('<%dq'%n, f.read(8*n))
        MJO=struct.unpack('<%dq'%n, f.read(8*n))
    return JMAX,MJ,MJO

def T(x): return (3*x+1)//2 if x&1 else x//2
def tst(x):
    s=0
    while x!=1:
        x=3*x+1 if x&1 else x//2; s+=1
    return s
_cache={}
def info(j):
    r=_cache.get(j)
    if r: return r
    x=j;h=0;o=0
    while x!=1:
        if x&1:o+=1
        x=T(x);h+=1
    _cache[j]=(h,o)
    return h,o

def tst_2n_plus_j(n,j):
    """exact tst(2^n+j) by the identity (needs j < 2^n)"""
    x=j;o=0
    for t in range(n):
        if x&1: o+=1
        x=T(x)
    return n + o + tst(3**o + x)

def window(phi,n,qmax=1200):
    """returns dict p -> set of good Q, and (A1,c1,L)"""
    A1=(n+1)//2; c1=1+(n&1); L=phi[(A1,c1)]
    good={0:set(),1:set()}
    for p in (0,1):
        r=(n-p)&1; cc=1+r
        for Q in range(-qmax,qmax+1):
            if (Q&1)!=p: continue
            A=(n+Q+r)//2
            if A<0: continue
            v=phi.get((A,cc))
            if v is not None and v==L: good[p].add(Q)
    return good,(A1,c1,L)

def class_of(j):
    h,o=info(j); return 2*o-h, h&1, h, o

def main():
    phi=load_phi(sys.argv[1] if len(sys.argv)>1 else '/tmp/phi8k.bin')
    JMAX,MJ,MJO=load_minj(sys.argv[2] if len(sys.argv)>2 else '/tmp/minj5e9.bin')
    # candidate hard numbers: records (A006877) + Frinkle boundaries + extra known
    cands=[]
    try:
        for line in open('/home/user/collatz//home/user/collatz/data/b006877.txt'):
            p=line.split()
            if len(p)>=2 and p[1].isdigit(): cands.append(int(p[1]))
    except IOError: pass
    try:
        import csv
        for r in csv.DictReader(open('/home/user/collatz//home/user/collatz/data/collatz_streaks.csv')):
            try: s=int(r['streak']); cands.append(s+1)
            except: pass
    except IOError: pass
    cands=sorted(set(cands))
    # precompute signatures
    sig={}
    for j in cands:
        if j < (1<<62): sig[j]=class_of(j)
    nlo=int(sys.argv[3]); nhi=int(sys.argv[4])
    out=[]
    for n in range(nlo,nhi+1):
        good,(A1,c1,L)=window(phi,n)
        # 1) smallest bad class rep from scan table
        best=None; src=None
        for p in (0,1):
            r=(n-p)&1
            for Q in range(-1200,1201):
                if (Q&1)!=p: continue
                if Q in good[p]: continue
                A=(n+Q+r)//2
                if A<0: continue
                j=0
                if -4000<=Q<4000:
                    j=MJ[(Q+4000)*2+p]
                if j>0 and (best is None or j<best): best=j; src='scan'
        # 2) candidates from hard numbers
        for j,(Q,p,h,o) in sig.items():
            if h>n: continue
            if Q not in good[p]:
                if best is None or j<best: best=j; src='hard'
        # 3) powers of two
        for t in range(1,121):
            j=1<<t
            Q=-t; p=t&1
            if Q<-600: 
                if best is None or j<best: best=j; src='2^k'
                continue
            if Q not in good[p]:
                if best is None or j<best: best=j; src='2^k'
        out.append((n,best,src))
    print("n, first_bad_j_ub, a(n)_ub, source, verified_by_direct_height")
    for n,best,src in out:
        if best is None:
            print(f"{n}, -, -, -, no breaker found"); continue
        # direct verification of the breaker
        ref=tst_2n_plus_j(n,1); t=tst_2n_plus_j(n,best)
        ok = (t!=ref)
        print(f"{n}, {best}, {best-1}, {src}, {'verified' if ok else 'FAILED'}")

if __name__=='__main__':
    main()
