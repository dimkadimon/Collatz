#!/usr/bin/env python3
"""boundary2.py -- certified interval analysis for a(n) = A277109(n).

For each n, the bad classes Q = 2o(j)-h(j) (with p = h mod 2) are those where the
landing value phi((n+Q+r)/2, 1+r) differs from phi((n+1)//2, 1+(n&1)); then
      a(n) = min{ minj(Q) : Q bad, h<=n } - 1.
Every class gets a LOWER bound (from exhaustive scans: none <= Jmax) and an UPPER bound
(scan hit, the proved-for-|Q|<=32 / heuristic-for-all-Q<0 rule minj=2^|Q|, climbed
witnesses, or a candidate pool).  We report the interval [lo,hi]; lo==hi means exact.

Output: n,a_lo,a_hi,best_Q,ub_src,exact,all_ub,note
"""
import numpy as np, struct, sys, csv, argparse, os

def load_phi(path):
    d=open(path,'rb').read(); Amax,=struct.unpack('<I',d[:4])
    return int(Amax), np.frombuffer(d,dtype='<u8',offset=4,count=2*(Amax+1))

def load_scans(paths):
    """returns exact[(Q,p)]=j, jmax = max Jmax among scans"""
    exact={}; jmax=0
    for p in paths:
        try: d=open(p,'rb').read()
        except IOError: continue
        if len(d)<8+16*16002: continue
        N=8001*2
        JM,=struct.unpack('<q',d[:8]); jmax=max(jmax,JM)
        MJ=np.frombuffer(d,dtype='<q',offset=8,count=N)
        for i in np.nonzero(MJ>0)[0]:
            Q=int(i//2-4000); par=int(i%2); j=int(MJ[i])
            if par!=abs(Q)%2: continue
            k=(Q,par)
            if k not in exact or j<exact[k]: exact[k]=j
    return exact,jmax

def T(x): return (3*x+1)//2 if x&1 else x//2
def info(j):
    x=j;h=0;o=0
    while x!=1:
        if x&1: o+=1
        x=T(x); h+=1
    return h,o

def load_climb(paths):
    ub={}
    for p in paths:
        try: f=open(p)
        except IOError: continue
        for line in f:
            if line.startswith('#'): continue
            t=line.split()
            if len(t)!=2: continue
            Q,j=int(t[0]),int(t[1])
            h,o=info(j)
            if 2*o-h!=Q: print(f"WARN {p}: {j} has Q={2*o-h} != {Q}",file=sys.stderr); continue
            k=(Q,h&1)
            if k not in ub or j<ub[k][0]: ub[k]=(j,h,'climb')
    return ub

def build_pool(qlim):
    pool={}
    def add(j,src,hmax):
        h,o=info(j); Q=2*o-h
        if abs(Q)>qlim: return
        k=(Q,h&1)
        cur=pool.get(k)
        if cur is None or j<cur[0]: pool[k]=(j,h,src)
    for line in open('/home/user/collatz/data/b006877.txt'):
        p=line.split()
        if len(p)>=2 and p[1].isdigit(): add(int(p[1]),'record',0)
    for r in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
        for col in ('B','streak'):
            try:
                v=int(r[col])
                if v>1: add(v+1 if col=='streak' else v,'frinkle',0)
            except Exception: pass
    return pool

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('phi'); ap.add_argument('out'); ap.add_argument('nlo',type=int); ap.add_argument('nhi',type=int)
    ap.add_argument('--qlim',type=int,default=800)
    ap.add_argument('--scan',nargs='*',default=[])
    ap.add_argument('--climb',nargs='*',default=[])
    ap.add_argument('--verify',type=int,default=0)
    ap.add_argument('--hmax',type=int,default=1)   # enforce h(candidate) <= n
    ap.add_argument('--refute',default='')   # old csv to compare against
    ap.add_argument('--oldplus1',type=int,default=0)
    a=ap.parse_args()
    Amax,V=load_phi(a.phi)
    exact,jmax=load_scans(a.scan)
    climb=load_climb(a.climb)
    pool=build_pool(a.qlim)
    old={}
    if a.refute:
        for r in csv.DictReader(open(a.refute)): old[int(r['n'])]=(int(r['a_n']) if r['a_n'] else None)
    print(f"phi Amax={Amax} scan-exact={len(exact)} jmax={jmax} climb={len(climb)} pool={len(pool)}",file=sys.stderr)
    Qs=np.arange(-a.qlim,a.qlim+1,dtype=np.int64)
    f=open(a.out,'w'); f.write("n,a_lo,a_hi,best_Q,ub_src,exact,all_ub,old_claim,verdict,note\n")
    stats=dict(exact=0,refuted=0,consistent=0,none=0,partial=0)
    for n in range(a.nlo,a.nhi+1):
        A1=(n+1)//2; c1=1+(n&1); L=V[2*A1+c1-1]
        r=(n-Qs)&1
        A=(n+Qs+r)//2
        ok=(A>=0)&(A<=Amax)
        vals=np.where(ok, V[np.where(ok,2*A+r,0)], 0)
        bad=((vals!=L)|(~ok))
        lo=None;hi=None;bestQ=None;bsrc=None;allub=True
        for i in np.nonzero(bad)[0]:
            Q=int(Qs[i]); par=abs(Q)%2 if Q<0 else Q%2
            h_rule=abs(Q)
            lb=jmax+1; ub=None; src=None
            key=(Q,par)
            if key in exact:
                hj=exact[key]
                ok_h=True
                if a.hmax: 
                    hh=info(hj)[0]
                    ok_h = hh<=n
                if ok_h: lb=ub=hj; src='scan'
            if Q<0:
                v=1<<(-Q)
                if h_rule<=n:
                    lb=min(lb,v); ub=v if ub is None else min(ub,v); src='rule2^Q'
            c=climb.get(key)
            if c and c[1]<=n:
                if ub is None or c[0]<ub: ub=c[0]; src=(src+'+climb') if (src and src!='climb') else 'climb'
            pp=pool.get(key)
            if pp and pp[1]<=n:
                if ub is None or pp[0]<ub: ub=pp[0]; src=(src+'+pool') if src else 'pool'
            if lo is None or lb<lo: lo=lb
            if ub is None: allub=False
            else:
                if hi is None or ub<hi: hi=ub; bestQ=Q; bsrc=src
        lo=lo-1 if lo is not None else 0
        hi=hi-1 if hi is not None else jmax-1
        if lo>hi: lo=hi
        oldclaim=old.get(n)
        verdict=''
        if oldclaim is not None:
            if hi<oldclaim: verdict='OLD-REFUTED'
            elif lo==hi==oldclaim: verdict='OLD-EXACT'
            elif lo<=oldclaim<=hi: verdict='OLD-CONSISTENT'
            else: verdict='OLD-BELOW-LO'
        exactflag = 'yes' if lo==hi else 'no'
        if lo==hi: stats['exact']+=1
        if verdict=='OLD-REFUTED': stats['refuted']+=1
        elif verdict in ('OLD-EXACT','OLD-CONSISTENT'): stats['consistent']+=1
        if not allub: stats['partial']+=1
        f.write(f"{n},{lo},{hi},{bestQ if bestQ is not None else ''},{bsrc or ''},{exactflag},{allub},{oldclaim if oldclaim is not None else ''},{verdict},\n")
    f.close()
    print("stats:",stats,file=sys.stderr)

if __name__=='__main__': main()
