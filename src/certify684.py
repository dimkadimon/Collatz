#!/usr/bin/env python3
"""certify684.py -- certify A277684(k) = least n with A277109(n) = 2^k - 1.

a(n) = 2^k - 1 requires:
  (i) class -k bad at n   [necessary; minj(-k) = 2^k]
  (ii) no bad class at n has minj <= 2^k - 1
Exclusion of a row n for target 2^k-1:
  E1 known exact value != target
  E2 a_hi(n) < target
  E3 a verified member j <= target-1 of a bad class at n
"""
import numpy as np, struct, sys, json
sys.path.insert(0, '/home/user/collatz/src')
import certify2 as C

QLIM = 1300
NMAX = 39199

# ---------- member table: Q -> list of (j,h) sorted ----------
def build_member_table():
    d = {}
    def add(j):
        h, o = C.info(j); Q = 2*o - h
        d.setdefault(Q, []).append((j, h))
    # rule members (minj(-k)=2^k)
    for k in range(1, QLIM+1):
        add(1 << k)
    # scan minima
    for (Q, p), j in C.scan.items():
        add(j)
    # records + climber + pool via certify2
    mem = C.build_members(['/tmp/climb_a.txt'])
    for Q, lst in mem.items():
        for j, h in lst:
            d.setdefault(Q, []).append((j, h))
    for Q in d:
        d[Q] = sorted(set(d[Q]))
    return d

MT = build_member_table()
print("member table classes:", len(MT), "entries:", sum(len(v) for v in MT.values()), file=sys.stderr)

# ---------- row data ----------
def Larr():
    n = np.arange(NMAX+1, dtype=np.int64)
    A1 = (n+1)//2
    idx = 2*A1 + (n & 1)
    return C.V[idx].astype(np.int64)

Lrow = Larr()

def class_idx(n, Q):
    m = n + Q
    return m + 2*(m & 1)

def bad_series(k):
    """boolean array over n: class -k bad at n (within table range)"""
    Q = -k
    n = np.arange(NMAX+1, dtype=np.int64)
    m = n + Q
    idx = m + 2*(m & 1)
    A = (idx - ((m) & 1))//2
    ok = (A >= 0) & (A <= C.Amax) & (idx >= 0) & (idx < len(C.V))
    v = np.where(ok, C.V[np.clip(idx, 0, len(C.V)-1)], -1)
    return ok & (v != Lrow)

def bad_classes_row(n, qlim=QLIM):
    return [Q for Q in range(-qlim, qlim+1)
            if (lambda m: (m+2*(m & 1)) <= 2*C.Amax and int(C.V[m + 2*(m & 1)]) != int(Lrow[n]))(n+Q)]

def qlim_for(n):
    # classes Q with A<=Amax: A=(n+Q+r)/2 <= Amax  ->  Q <= 2*Amax - n
    return min(QLIM, max(0, 2*C.Amax - n))

def exclusion_member(n, target):
    """least known member j <= target-1 of a bad class at n (h<=n), or None"""
    best = None
    # candidate Q values: those with any known member <= target-1
    candQ = [Q for Q, lst in MT.items() if lst[0][0] <= target-1]
    for Q in candQ:
        if abs(Q) > qlim_for(n): continue
        m = n + Q
        if m + 2*(m & 1) > 2*C.Amax: continue
        if int(C.V[m + 2*(m & 1)]) != int(Lrow[n]):   # bad
            for j, h in MT[Q]:
                if j <= target-1 and h <= n:
                    if best is None or j < best: best = j
                    break
    return best

def a_interval(n):
    r = C.bad_interval(n, MT)
    return (r[0], r[1], r[2]) if r else (None, None, False)

def main():
    pub = {}
    for l in open('/home/user/collatz/data/b277684.txt'):
        p = l.split()
        if len(p) >= 2: pub[int(p[0])] = int(p[1])
    out = {}
    for k in range(1, 700):
        target = (1 << k) - 1
        rows = np.nonzero(bad_series(k))[0]
        info = dict(k=k, first=int(rows[0]) if len(rows) else None, rows=[int(x) for x in rows[:40]])
        status = 'NONE'
        answer = None
        for n in rows:
            n = int(n)
            if n == 0: continue
            lo, hi, ex = a_interval(n)
            if ex and hi != target:
                continue                      # E1: exact value differs
            if hi is not None and hi < target:
                continue                      # E2
            mem = exclusion_member(n, target)
            if mem is not None:
                continue                      # E3
            # unresolved candidate
            if ex and hi == target:
                status = 'EXACT'; answer = n
            else:
                status = 'UBOUND'; answer = n
            break
        info.update(status=status, answer=answer,
                    pub=pub.get(k), match=('PUBLISHED-AGREE' if (answer is not None and pub.get(k) == answer) else
                                           ('PUBLISHED-DIFF' if k in pub else '')))
        out[k] = info
        if k <= 75 or k in range(205, 235) or k % 50 == 0:
            print(k, info['status'], info['answer'], 'pub=', info['pub'], info['match'], 'firstRows=', info['rows'][:6], file=sys.stderr)
    json.dump({str(k): v for k, v in out.items()}, open('/tmp/cert684.json', 'w'))
    # validation vs published
    ok = bad = 0
    for k, v in out.items():
        if v['pub'] is not None and v['answer'] is not None:
            if v['pub'] == v['answer']: ok += 1
            else:
                bad += 1
                print("MISMATCH k=%d ours=%s pub=%s status=%s" % (k, v['answer'], v['pub'], v['status']), file=sys.stderr)
    print("published k compared:", ok+bad, "agree:", ok, "mismatch:", bad, file=sys.stderr)

if __name__ == '__main__':
    main()
