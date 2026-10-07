#!/usr/bin/env python3
"""climber.py -- random-restart greedy climbs of Q in the accelerated reverse tree.
State: odd y with class Q; steps y -> (2^a y - 1)/3, Q -> Q + (2-a).
Policy: take a=1 (Q+1) when y%3==2 and the result stays alive (not 0 mod 3);
        take a=2 (Q unchanged, y*=4/3) when y%3==1 and the result stays alive;
        otherwise escape with a larger a (costing Q) or stop.
Records the smallest j found for every Q level reached (upper bounds on minj(Q)).
"""
import random, sys, csv, time

def info(j):
    x=j;h=0;o=0
    while x!=1:
        if x&1: x=(3*x+1)//2; o+=1
        else: x//=2
        h+=1
    return h,o

def Q_of(j):
    h,o=info(j); return 2*o-h

def seeds():
    S=set()
    for line in open('/home/user/collatz/data/b006877.txt'):
        p=line.split()
        if len(p)>=2 and p[1].isdigit(): S.add(int(p[1]))
    for r in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
        for col in ('B','streak'):
            try:
                v=int(r[col])
                if v>1: S.add(v); S.add(v+1)
            except Exception: pass
    return sorted(S)

def climb(y, Q, target_hi, rng, maxsteps=30000, best=None, jmax=1<<600):
    if best is None: best={}
    for _ in range(maxsteps):
        if Q>target_hi: return best
        if Q not in best or y<best[Q]: best[Q]=y
        r=y%3
        if r==0: return best                      # dead
        opts=[]
        if r==2:
            opts=[1,3,5,7,9,11]
        else:
            opts=[2,4,6,8,10,12]
        alive=[]; dead=[]
        for a in opts:
            y2=((1<<a)*y-1)//3
            if y2<3 or y2>jmax: continue
            (dead if y2%3==0 else alive).append((a,y2))
        pick=None
        if alive:
            pick=min(alive, key=lambda c:c[0]) if rng.random()<0.85 else rng.choice(alive)
        elif dead:
            Qn=Q+2-dead[0][0]
            pick = dead[0] if Qn==target_hi else (rng.choice(dead) if rng.random()<0.3 else None)
            if pick is None: return best
        if pick is None: return best
        a,y=pick; Q=Q+2-a
        if Q>target_hi: return best
    return best

if __name__=='__main__':
    target_hi=int(sys.argv[1]) if len(sys.argv)>1 else 330
    seconds=float(sys.argv[2]) if len(sys.argv)>2 else 60.0
    S=seeds()
    byQ=sorted(((Q_of(j),j) for j in S), reverse=True)
    rng=random.Random(int(sys.argv[3]) if len(sys.argv)>3 else 7)
    best={}
    t0=time.time(); attempts=0
    pool=[j for q,j in byQ[:60]]
    while time.time()-t0 < seconds:
        attempts+=1
        y0=rng.choice(pool)
        Q0=Q_of(y0)
        climb(y0,Q0,target_hi,rng,maxsteps=20000,best=best)
        if attempts%200==0:
            print(f"# attempts={attempts} levels={len(best)} maxQ={max(best)} t={time.time()-t0:.0f}s",file=sys.stderr)
    out='/tmp/climb_%d.txt'%target_hi
    with open(out,'w') as f:
        for Q in sorted(best): f.write(f"{Q} {best[Q]}\n")
    print(f"wrote {out}: levels {min(best)}..{max(best)}, {attempts} attempts")
    for Q in sorted(best):
        if Q>=max(best)-15: print(f"  Q={Q}: j<={best[Q]} (bits {best[Q].bit_length()})")
