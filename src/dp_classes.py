#!/usr/bin/env python3
"""dp_classes.py -- exact minimum j per signature class (h,o) / (Q,p) by reverse-tree DP.

class(j): h = #accelerated steps to reach 1, o = #odd values before 1, Q = 2o-h, p = h mod 2.
Inverse moves from 1:  E: y -> 2y  (undo an even step);  O: y -> (2y-1)/3 (undo an odd step,
allowed iff y % 3 == 2).  Both maps are strictly increasing, and operability depends only on
y mod 3, so for each (steps taken, o used, y mod 3) only the minimum y can ever be useful.
Paths that re-enter 1 are discarded (they belong to a shorter class).  Hence

    minj(Q,p) = min over h of DPmin[h][o=(h+Q)/2]   is EXACT.

Writes a text table:  Q p minj h o   (sorted by Q,p) and prints stats.
"""
import sys, time

def main():
    H=int(sys.argv[1]) if len(sys.argv)>1 else 2100
    outpath=sys.argv[2] if len(sys.argv)>2 else '/home/user/collatz/data/classmin.txt'
    t0=time.time()
    cur=[[None,1,None]]          # o -> [min value with v%3==0, ==1, ==2]
    best={}                      # (Q,p) -> (minj,h,o)
    # h=0: value 1, o=0 -> Q=0, p=0
    best[(0,0)]=(1,0,0)
    for h in range(0,H):
        nxt=[[None,None,None] for _ in range(len(cur)+1)]
        for o,row in enumerate(cur):
            for res3,v in enumerate(row):
                if v is None: continue
                # E: 2v
                v2=v+v; r=v2%3
                if nxt[o][r] is None or v2<nxt[o][r]: nxt[o][r]=v2
                # O: (2v-1)/3, only if v%3==2 ; skip if result is 1 (re-entering 1)
                if res3==2:
                    v2=(2*v-1)//3
                    assert (2*v-1)%3==0 and v2>=1
                    if v2!=1:
                        r=v2%3
                        if nxt[o+1][r] is None or v2<nxt[o+1][r]: nxt[o+1][r]=v2
        cur=nxt
        hh=h+1; p=hh&1
        for o,row in enumerate(cur):
            Q=2*o-hh
            M=None
            for v in row:
                if v is not None and (M is None or v<M): M=v
            if M is None: continue
            k=(Q,p)
            e=best.get(k)
            if e is None or M<e[0]: best[k]=(M,hh,o)
        if hh%200==0: print(f"  h={hh} classes={len(best)} t={time.time()-t0:.0f}s",file=sys.stderr)
    with open(outpath,'w') as f:
        f.write("# Q p minj h o   (exact class minimum, DP to h<=%d)\n"%H)
        for (Q,p),(M,hh,o) in sorted(best.items()):
            f.write(f"{Q} {p} {M} {hh} {o}\n")
    print(f"wrote {outpath}: {len(best)} classes, {time.time()-t0:.1f}s",file=sys.stderr)

if __name__=='__main__': main()
