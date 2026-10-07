#!/usr/bin/env python3
"""walkup.py -- raise Q by walking backward in the accelerated reverse tree from known
high-Q seeds.  Reverse step: y -> (2^a*y-1)/3, valid for a odd iff y%3==2, a even iff y%3==1;
new class Q = old Q + (2-a) (a=1 gives +1, a=2 gives 0, a=3 gives -1, ...), new height h += a.
Goal: reach a target class (Q,p) with p = Q mod 2 and h <= nmax, keeping j as small as we can.
"""
import random, sys, csv, heapq

def info(j):
    x=j;h=0;o=0
    while x!=1:
        if x&1: x=(3*x+1)//2; o+=1
        else: x//=2
        h+=1
    return h,o

def Q_of(j):
    h,o=info(j); return 2*o-h, h&1, h

def seeds():
    S=[]
    for line in open('/home/user/collatz/data/b006877.txt'):
        p=line.split()
        if len(p)>=2 and p[1].isdigit(): S.append(int(p[1]))
    try:
        for r in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
            for col in ('B','streak'):
                try:
                    v=int(r[col])
                    if v>1: S.append(v); S.append(v+1)
                except Exception: pass
    except IOError: pass
    return sorted(set(S))

def walk(target_Q, seed, rng, hmax=40000, jmax=None, nodes=4000):
    """randomized best-first: from seed, try to hit class (target_Q). Return j or None."""
    h0,o0 = info(seed); Q0 = 2*o0-h0
    stack=[(seed,Q0,h0)]
    best=None
    visited=set()
    n=0
    while stack and n<nodes:
        n+=1
        y,Q,h = stack.pop()
        if Q==target_Q:
            if best is None or y<best: best=y
            continue
        if (y,Q) in visited: continue
        visited.add((y,Q))
        if Q>target_Q: continue
        r=y%3
        if r==0: continue
        cands=[]
        for a in (1,2,3,4,5,6,7,8):
            if (r==2) != (a%2==1): continue
            y2=((1<<a)*y-1)//3
            if y2<3: continue
            Q2=Q+2-a
            if Q2>target_Q: continue
            h2=h+a
            if h2>hmax: continue
            if jmax and y2>jmax: continue
            cands.append((y2,Q2,h2,a))
        cands.sort(key=lambda c:(c[1],-c[0]))   # prefer Q growth, then small y
        for c in cands[::-1]: stack.append(c[:3])
    return best,n

if __name__=='__main__':
    T=int(sys.argv[1]); tries=int(sys.argv[2]) if len(sys.argv)>2 else 3
    S=seeds()
    print(f"# {len(S)} seeds; target Q={T}",file=sys.stderr)
    rng=random.Random(12345)
    # order seeds by Q desc (closest to target)
    S2=sorted(S,key=lambda j:-Q_of(j)[0])[:400]
    best=None
    for s in S2:
        j,n=walk(T,s,rng)
        if j is not None and (best is None or j<best):
            best=j
            print(f"# seed {s} -> j={j} bits={j.bit_length()}",file=sys.stderr)
    if best:
        h,o=info(best)
        print(f"Q={T}: witness j={best} (bits {best.bit_length()}) h={h} o={o} Qcheck={2*o-h}")
    else:
        print(f"Q={T}: none")
