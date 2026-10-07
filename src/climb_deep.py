#!/usr/bin/env python3
"""climb_deep.py -- iterative-deepening climber: repeatedly feeds its own best finds back
as seeds, recording the smallest j found per class Q (upper bounds on minj(Q))."""
import random, sys, time, csv

def info(j):
    x=j;h=0;o=0
    while x!=1:
        if x&1: x=(3*x+1)//2; o+=1
        else: x//=2
        h+=1
    return h,o
def Q_of(j):
    h,o=info(j); return 2*o-h

def load_seeds():
    S=set([5,21,85,341,1365])
    for line in open('/home/user/collatz/data/b006877.txt'):
        p=line.split()
        if len(p)>=2 and p[1].isdigit(): S.add(int(p[1]))
    try:
        for r in csv.DictReader(open('/home/user/collatz/data/collatz_streaks.csv')):
            for col in ('B','streak'):
                try:
                    v=int(r[col])
                    if v>1: S.add(v+v-1 if False else v); S.add(v+1)
                except Exception: pass
    except IOError: pass
    return S

def climb(y,Q,hi,rng,best,maxsteps=40000,jmax=1<<700):
    for _ in range(maxsteps):
        if Q>hi: return
        if Q not in best or y<best[Q][0]: best[Q]=(y,Q_of(y))
        r=y%3
        if r==0: return
        alive=[];dead=[]
        for a in range(1 if r==2 else 2, 16, 2):
            y2=((1<<a)*y-1)//3
            if y2<3 or y2>jmax: continue
            (dead if y2%3==0 else alive).append((a,y2))
        pick=None
        if alive:
            alive.sort()
            pick = alive[0] if rng.random()<0.8 else rng.choice(alive[:3])
        elif dead:
            pick = rng.choice(dead)
        if pick is None: return
        a,y2=pick
        y=y2; Q=Q+2-a

def main():
    hi=int(sys.argv[1]); seconds=float(sys.argv[2]); seed=int(sys.argv[3]) if len(sys.argv)>3 else 1
    rng=random.Random(seed)
    S=load_seeds()
    best={}
    t0=time.time(); rounds=0
    frontier=sorted(S,key=lambda j:-Q_of(j))[:200]
    while time.time()-t0<seconds:
        rounds+=1
        for y0 in frontier:
            climb(y0,Q_of(y0),hi,rng,best,maxsteps=3000)
        # feed back the best finds as new frontier
        mq=max(best) if best else 0
        items=sorted(((Q_of(j),j) for j in [v[0] for v in best.values()] if Q_of(j)>=mq-60), reverse=True)
        frontier=[j for q,j in items[:150]]
        frontier+=sorted(S,key=lambda j:-Q_of(j))[:50]
        if not frontier: frontier=sorted(S,key=lambda j:-Q_of(j))[:200]
        if best:
            mq=max(best)
            print(f"# round={rounds} levels={len(best)} maxQ={mq} t={time.time()-t0:.0f}s",file=sys.stderr)
    out=sys.argv[4] if len(sys.argv)>4 else '/tmp/climb_deep.txt'
    with open(out,'w') as f:
        for Q in sorted(best): f.write(f"{Q} {best[Q][0]}\n")
    print("wrote",out,"levels",min(best),"..",max(best))

if __name__=='__main__': main()
