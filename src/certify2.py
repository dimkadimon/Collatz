#!/usr/bin/env python3
"""certify2.py -- definitive certified intervals for A277109(n) = a(n), using the
validated convention (classes.py: L = V[2*ceil(n/2) + (n&1)]; class idx = 2*((n+Q+r)//2)+r,
r=(n-Q)&1; bad <=> idx value != L).

Bounds per bad class Q:
  lb(Q): 2^|Q| for Q<0 (minj(-j)=2^j, proved); scan minj when j<=5e9 covered; else 5e9+1
  ub(Q): same as lb when exact; else least verified member (rule/scan/records/climb)
a(n) = min over bad Q of minj(Q) - 1  (requires h<=n per member)
"""
import numpy as np, struct, sys, csv

PHI='/tmp/phi22k.bin'; SCAN='/tmp/minj5e9.bin'
QLIM=1300

def load():
    d=open(PHI,'rb').read(); Amax,=struct.unpack('<I',d[:4])
    V=np.frombuffer(d,dtype='<u8',offset=4,count=2*(Amax+1))
    dd=open(SCAN,'rb').read(); jmax,=struct.unpack('<q',dd[:8])
    MJ=np.frombuffer(dd,dtype='<q',offset=8,count=8001*2)
    scan={}
    for i in np.nonzero(MJ>0)[0]:
        scan[(int(i//2-4000),int(i%2))]=int(MJ[i])
    return Amax,V,jmax,scan

Amax,V,jmax,scan=load()

def info(j):
    x=j;h=0;o=0
    while x!=1:
        if x&1: o+=1; x=(3*x+1)//2
        else: x>>=1
        h+=1
    return h,o

def build_members(extra_paths=()):
    mem={}  # Q -> sorted [(j,h)]
    def add(j):
        if j<2: return
        h,o=info(j); Q=2*o-h
        if -QLIM<=Q<=QLIM: mem.setdefault(Q,[]).append((j,h))
    recs=[int(l.split()[1]) for l in open('/home/user/collatz/data/b006877.txt') if len(l.split())>=2 and l.split()[1].isdigit()]
    for r in recs:
        for w in (r-1,r,r+1,2*r):
            add(w)
        x=r
        for _ in range(3):
            x=(3*x+1)//2 if x&1 else x>>1
            add(x)
    for row in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
        for col in ('B','streak','H'):
            try:
                v=int(row[col])
                if v>1: add(v); add(v+1)
            except Exception: pass
    for p in extra_paths:
        try:
            for line in open(p):
                if line.startswith('#'): continue
                a=line.split()
                if len(a)>=2:
                    Q=int(a[0]); j=int(a[1])
                    add(j)
        except IOError: pass
    for Q in mem: mem[Q]=sorted(mem[Q])
    return mem

def Lval(n):
    return int(V[2*((n+1)//2) + (n&1)])

def bad_interval(n, mem):
    """returns (a_lo,a_hi,exact,detail list of (Q,lb,ub,src))"""
    L=Lval(n)
    lo=None; hi=None; det=[]
    for Q in range(-QLIM,QLIM+1):
        r=(n-Q)&1; A=(n+Q+r)//2
        if A<0 or A>Amax: continue
        if int(V[2*A+r])==L: continue     # good class
        lb=None; ub=None; src='gap'
        if Q<0 and (-Q)<=n:
            lb=ub=1<<(-Q); src='rule'
        k=(Q,abs(Q)%2)
        if k in scan:
            j=scan[k]
            if info(j)[0]<=n:
                lb=j if lb is None else min(lb,j); ub=j if ub is None else min(ub,j)
                if src=='gap': src='scan'
        if Q in mem:
            for j,h in mem[Q]:
                if h<=n:
                    if ub is None or j<ub: ub=j; src='mem'
                    break
        if lb is None: lb=jmax+1
        if ub is None: ub=None
        lo=lb if lo is None else min(lo,lb)
        if ub is not None: hi=ub if hi is None else min(hi,ub)
        det.append((Q,lb,ub,src))
    if hi is None: return None
    return (lo-1, hi-1, lo==hi, det)

if __name__=='__main__':
    cmd=sys.argv[1] if len(sys.argv)>1 else 'validate'
    mem=build_members([p for p in sys.argv[2:] if p.endswith('.txt')])
    if cmd=='validate':
        pub={}
        for l in open('/home/user/collatz/data/b277109.txt'):
            p=l.split()
            if len(p)>=2: pub[int(p[0])]=int(p[1])
        bad_lo=bad_hi=exact=match_hi=0
        exs=[]
        for n in sorted(pub):
            r=bad_interval(n,mem)
            if r is None: continue
            lo,hi,ex,_=r
            if lo>pub[n] or hi<pub[n]:
                bad_lo+=1 if lo>pub[n] else 0
                bad_hi+=1 if hi<pub[n] else 0
                if len(exs)<10: exs.append((n,lo,hi,pub[n]))
            if ex: exact+=1
            if hi==pub[n]: match_hi+=1
        print(f"n=0..3953: contradictions lo>pub {bad_lo}, hi<pub {bad_hi}; exact(a_lo==a_hi) rows {exact}; a_hi==published {match_hi}/{len(pub)}")
        for e in exs: print("  CONTRADICTION",e)
    elif cmd=='row':
        for a in sys.argv[2:]:
            n=int(a); r=bad_interval(n,mem)
            lo,hi,ex,det=r
            print(f"n={n}: a in [{lo},{hi}] exact={ex}")
            det.sort(key=lambda t:(t[2] if t[2] is not None else 1<<200))
            for Q,lb,ub,src in det[:12]:
                print(f"   Q={Q:5d} lb={lb} ub={ub} ({src})")
    elif cmd=='scan':
        nlo,nhi=int(sys.argv[2]),int(sys.argv[3])
        out=open(sys.argv[4],'w') if len(sys.argv)>4 else sys.stdout
        if out is not sys.stdout: out.write("n,a_lo,a_hi,exact\n")
        for n in range(nlo,nhi+1):
            r=bad_interval(n,mem)
            if r is None: continue
            lo,hi,ex,_=r
            print(f"{n},{lo},{hi},{'yes' if ex else 'no'}",file=out)
        if out is not sys.stdout: out.close()
