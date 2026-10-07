#!/usr/bin/env python3
"""verify_tops.py -- for the top of each long claimed chain, take the corrected upper-bound
breaker j from the boundary2 table and verify directly that tst(2^n+1) != tst(2^n+j)."""
import csv, sys, time
def tst(x):
    s=0
    while x!=1:
        if x&1: x=(3*x+1)>>1
        else: x>>=1
        s+=1
    return s
old={int(r['n']):r for r in csv.DictReader(open('/tmp/b_ext.csv'))}
new={int(r['n']):r for r in csv.DictReader(open('/tmp/b2_ext.csv'))}
# chain tops: n where old a(n) is a power of two minus 1 and its successor's claim is smaller
tops=[]
for n in sorted(old):
    a=old[n]['a_n']
    if not a: continue
    A=int(a)
    if A>0 and (A & (A+1))==0 and A>=2**30:      # a = 2^k - 1
        nxt=old.get(n+1)
        if nxt and nxt['a_n'] and int(nxt['a_n'])<A: tops.append(n)
print("claimed chain tops with a(n)=2^k-1 >= 2^30:",len(tops))
sel=[]
for n in tops:
    try: k=(int(old[n]['a_n'])+1).bit_length()-1
    except Exception: k=0
    sel.append((k,n))
sel.sort(reverse=True)
print("top 12 by k:",[ (k,n) for k,n in sel[:12]])
ok=ref=0
print("\n#  n      old claim 2^k-1        new a_hi           direct tst check")
for k,n in sel[:12]:
    r=new[n]
    hi=r['a_hi']; j=int(hi)+1 if hi else None
    if j is None: print(f"  n={n}: no new bound"); continue
    t0=time.time()
    t1=tst((1<<n)+1); t2=tst((1<<n)+j)
    verdict = "REFUTED (breaker found)" if t1!=t2 else "no breaker(!)"
    if t1!=t2: ok+=1
    print(f"  n={n}  2^{k}-1   {hi}   tst(2^n+1)={t1} tst(2^n+j)={t2} -> {verdict}  [{time.time()-t0:.1f}s]")
