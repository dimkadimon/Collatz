#!/usr/bin/env python3
"""boundary2.py -- certified intervals for a(n) = A277109(n).

a(n) = min{ m>=1 : tst(2^n+m) != tst(2^n+1) } - 1   (b-file scale; verified vs b-file).

Class model: tst(2^n+j) depends on the class (h(j),o(j)) via
    Q = 2o-h ,  p = h mod 2 ,  r = (n-h)&1 ,  A = (n+Q+r)//2 ,  c = 1+r
    tst(2^n+j) = (n+h-r)/2 + phi(A,c)          [VALID ONLY IF h(j) <= n]
So a class is "bad" iff phi(A,c) != phi(A1,c1) (with (A1,c1) the class of j=1).
For class-min candidates with h(j) > n the phi comparison is NOT valid: badness is
decided by direct computation of tst(2^n+j) (tst_2npj).

a(n) = min over bad classes of (class-min member j) - 1.
Candidate class-mins come from: exhaustive minj scans (exact <= Jmax), the proved
rule minj(Q)=2^|Q| for Q<0, climbed witnesses, and a candidate pool (upper bounds).

Output: n,a_lo,a_hi,best_Q,ub_src,exact,all_ub,old_claim,verdict,note
"""
import numpy as np, struct, sys, csv, argparse

def load_phi(path):
    d=open(path,'rb').read(); Amax,=struct.unpack('<I',d[:4])
    return int(Amax), np.frombuffer(d,dtype='<u8',offset=4,count=2*(Amax+1))

def load_scans(paths):
    """exact[(Q,p)] = smallest j found; jmax = largest scan limit"""
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
            k=(Q,par)
            if k not in exact or j<exact[k]: exact[k]=j
    return exact,jmax

def T(x): return (3*x+1)//2 if x&1 else x//2
_ic={}
def info(j):
    r=_ic.get(j)
    if r is not None: return r
    x=j;h=0;o=0
    while x!=1:
        if x&1: o+=1
        x=T(x); h+=1
    _ic[j]=(h,o); return h,o
def tst(x):
    s=0
    while x!=1:
        x=3*x+1 if x&1 else x>>1
        s+=1
    return s
def tst_2npj(n,j):
    """exact standard-step count of 2^n+j via n accelerated steps on j"""
    x=j;o=0
    for _ in range(n):
        if x&1: o+=1
        x=T(x)
    return n+o+tst(3**o+x)

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
    def add(j,src):
        h,o=info(j); Q=2*o-h
        if abs(Q)>qlim: return
        k=(Q,h&1)
        cur=pool.get(k)
        if cur is None or j<cur[0]: pool[k]=(j,h,src)
    for line in open('/home/user/collatz/data/b006877.txt'):
        p=line.split()
        if len(p)>=2 and p[1].isdigit(): add(int(p[1]),'record')
    for r in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
        for col in ('B','streak'):
            try:
                v=int(r[col])
                if v>1: add(v+1 if col=='streak' else v,'frinkle')
            except Exception: pass
    return pool

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('phi'); ap.add_argument('out'); ap.add_argument('nlo',type=int); ap.add_argument('nhi',type=int)
    ap.add_argument('--qlim',type=int,default=800)
    ap.add_argument('--scan',nargs='*',default=[])
    ap.add_argument('--climb',nargs='*',default=[])
    ap.add_argument('--refute',default='')
    a=ap.parse_args()
    Amax,V=load_phi(a.phi)
    scan,jmax=load_scans(a.scan)
    climb=load_climb(a.climb)
    pool=build_pool(a.qlim)
    old={}
    if a.refute:
        for r in csv.DictReader(open(a.refute)):
            try: old[int(r['n'])]=(int(r['a_n']) if r.get('a_n') not in (None,'') else None)
            except Exception: pass
    print(f"phi Amax={Amax} scan={len(scan)} jmax={jmax} climb={len(climb)} pool={len(pool)} old={len(old)}",file=sys.stderr)
    Qs=np.arange(-a.qlim,a.qlim+1,dtype=np.int64)
    NQ=len(Qs)
    # ---- per-class candidates (n-independent): rule(2^|Q| for Q<0), scan, climb, pool
    jmin=[None]*NQ; hmin=[0]*NQ; lb=[0]*NQ; src=['']*NQ
    for i,Q in enumerate(Qs):
        Q=int(Q); cands=[]
        if Q<0: cands.append((1<<(-Q),abs(Q),'rule2^Q'))
        k=(Q,abs(Q)%2)
        if k in scan: j=scan[k]; cands.append((j,info(j)[0],'scan'))
        if k in climb: cands.append(climb[k][:2]+('climb',))
        if k in pool:  cands.append(pool[k][:2]+('pool',))
        if cands:
            j,h,src_=min(cands)
            jmin[i]=j; hmin[i]=h; src[i]=src_
            lb[i]=min([c[0] for c in cands if c[2] in ('rule2^Q','scan')]) if any(c[2] in ('rule2^Q','scan') for c in cands) else (jmax+1 if jmax>0 else None)
        else:
            lb[i]=(jmax+1 if jmax>0 else None)
    order=[i for i in sorted(range(NQ), key=lambda i:(jmin[i] is None, jmin[i] if jmin[i] is not None else 0))]
    f=open(a.out,'w'); f.write("n,a_lo,a_hi,best_Q,ub_src,exact,all_ub,old_claim,verdict,note\n")
    stats=dict(exact=0,refuted=0,consistent=0,none=0,partial=0,dircheck=0)
    for n in range(a.nlo,a.nhi+1):
        A1=(n+1)//2; c1=1+(n&1); L=V[2*A1+c1-1]
        r=(n-Qs)&1
        A=(n+Qs+r)//2
        ok=(A>=0)&(A<=Amax)
        vals=np.where(ok, V[np.where(ok,2*A+r,0)], 0)
        badphi=((vals!=L)|(~ok))
        # a_hi : ascending class-min scan, first verified-bad class
        best=None; bestQ=None; bsrc=''; note=[]; ndir=0; nbad=0
        excl=set(); addbad=set()
        for i in order:
            if jmin[i] is None: break
            if best is not None and jmin[i]>=best: break
            Q=int(Qs[i]); h=hmin[i]
            if h<=n:
                if badphi[i]:
                    nbad+=1
                    if best is None: best=jmin[i]; bestQ=Q; bsrc=src[i]
                    else: break
            else:
                j=jmin[i]
                t1=tst_2npj(n,j); t0=tst_2npj(n,1); ndir+=1
                if t1!=t0:
                    nbad+=1; addbad.add(i)
                    if best is None: best=jmin[i]; bestQ=Q; bsrc=src[i]+'+dir'
                    else: break
                else:
                    excl.add(i)
        # a_lo : min lower bound over phi-bad classes, minus exonerated, plus direct-bad
        bidx=np.nonzero(badphi)[0]
        lo=None
        for i in list(bidx)+list(addbad):
            if i in excl: continue
            v=lb[i]
            if v is not None and (lo is None or v<lo): lo=v
        if ndir: note.append('dir')
        nocand=np.array([jmin[i] is None for i in range(NQ)])
        allub = not (badphi & ~nocand).any() or not (badphi & ~nocand).any() if False else not (badphi & nocand & np.array([i not in excl for i in range(NQ)])).any()
        if best is None:
            note.append('nobad'); stats['none']+=1
            hi=''; bestQ=''; bsrc=''
            lo=''
        else:
            stats['dircheck']+=ndir
            a_lo = (lo-1) if lo is not None else 0
            a_hi = best-1
            if a_lo>a_hi: a_lo=a_hi
            hi=a_hi; lo=a_lo
        oldclaim=old.get(n)
        verdict=''
        if oldclaim is not None and best is not None:
            if hi<oldclaim: verdict='OLD-REFUTED'
            elif lo==hi==oldclaim: verdict='OLD-EXACT'
            elif lo<=oldclaim<=hi: verdict='OLD-CONSISTENT'
            else: verdict='OLD-BELOW-LO'
        ex='yes' if (best is not None and lo==hi) else 'no'
        if ex=='yes': stats['exact']+=1
        if verdict=='OLD-REFUTED': stats['refuted']+=1
        elif verdict in ('OLD-EXACT','OLD-CONSISTENT'): stats['consistent']+=1
        if not allub: stats['partial']+=1
        if not allub: note.append('nocand')
        f.write(f"{n},{lo},{hi},{bestQ},{bsrc},{ex},{allub},{oldclaim if oldclaim is not None else ''},{verdict},{'+'.join(note)}\n")
    f.close()
    print("stats:",stats,file=sys.stderr)

if __name__=='__main__': main()
