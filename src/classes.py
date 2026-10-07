#!/usr/bin/env python3
"""classes.py -- shared machinery: phi table, class minima, class model of A277109.

Model (validated in boundary2.py / paper):
  a(n) = min{ m>=1 : tst(2^n+m) != tst(2^n+1) } - 1
  class of j: (h,o) = (#accelerated T-steps, #odd values in accelerated orbit),
              Q = 2o - h,  p = h mod 2
  tst(2^n+j) = (n+h-r)/2 + phi(A,c),  r = (n-h)&1, A = (n+Q+r)//2, c = 1+r
  (valid when h <= n; otherwise direct computation)
  class (Q,p) is "bad at n" iff phi(A,c) != phi(A1,c1) where (A1,c1) is the class of j=1.
  a(n) = least member of a bad class, minus 1.
"""
import struct, numpy as np, csv, sys

_ic = {}

def T(x): return (3*x+1)//2 if x & 1 else x >> 1

def info(j):
    r = _ic.get(j)
    if r is not None: return r
    x = j; h = 0; o = 0
    while x != 1:
        if x & 1: o += 1
        x = T(x); h += 1
    _ic[j] = (h, o); return h, o

def Q_of(j):
    h, o = info(j); return 2*o - h

def tst(x):
    s = 0
    while x != 1:
        x = 3*x+1 if x & 1 else x >> 1
        s += 1
    return s

def tst_2npj(n, j):
    """exact standard-step count of 2^n+j via n accelerated steps on j (needs h(j)>=n-ish)."""
    x = j; o = 0
    for _ in range(n):
        if x & 1: o += 1
        x = T(x)
    return n + o + tst(3**o + x)

def load_phi(path):
    d = open(path, 'rb').read(); Amax, = struct.unpack('<I', d[:4])
    return int(Amax), np.frombuffer(d, dtype='<u8', offset=4, count=2*(Amax+1))

def load_scans(paths):
    """exact[(Q,p)] = smallest member found; jmax = largest fully scanned limit"""
    exact = {}; jmax = 0
    for p in paths:
        try: d = open(p, 'rb').read()
        except IOError: continue
        if len(d) < 8 + 16*16002: continue
        N = 8001*2
        JM, = struct.unpack('<q', d[:8]); jmax = max(jmax, JM)
        MJ = np.frombuffer(d, dtype='<q', offset=8, count=N)
        for i in np.nonzero(MJ > 0)[0]:
            Q = int(i//2 - 4000); par = int(i % 2); j = int(MJ[i])
            k = (Q, par)
            if k not in exact or j < exact[k]: exact[k] = j
    return exact, jmax

def load_climb(paths):
    ub = {}
    for p in paths:
        try: f = open(p)
        except IOError: continue
        for line in f:
            if line.startswith('#'): continue
            t = line.split()
            if len(t) != 2: continue
            Q, j = int(t[0]), int(t[1])
            h, o = info(j)
            if 2*o - h != Q:
                print(f"WARN {p}: {j} has Q={2*o-h} != {Q}", file=sys.stderr); continue
            k = (Q, h & 1)
            if k not in ub or j < ub[k][0]: ub[k] = (j, h, 'climb')
    return ub

def build_pool(qlim, streak_csv='/home/user/collatz/data/collatz_streaks.csv',
               records='/home/user/collatz/data/b006877.txt'):
    pool = {}
    def add(j, src):
        h, o = info(j); Q = 2*o - h
        if abs(Q) > qlim: return
        k = (Q, h & 1)
        cur = pool.get(k)
        if cur is None or j < cur[0]: pool[k] = (j, h, src)
    try:
        for line in open(records):
            p = line.split()
            if len(p) >= 2 and p[1].isdigit(): add(int(p[1]), 'record')
    except IOError: pass
    try:
        for r in csv.DictReader(open(streak_csv)):
            for col in ('B', 'streak'):
                try:
                    v = int(r[col])
                    if v > 1: add(v+1 if col == 'streak' else v, 'frinkle')
                except Exception: pass
    except IOError: pass
    return pool

class Model:
    def __init__(self, phi_path, scan_paths=(), climb_paths=(), qlim=800):
        self.Amax, self.V = load_phi(phi_path)
        self.scan, self.jmax = load_scans(list(scan_paths))
        self.climb = load_climb(list(climb_paths))
        self.pool = build_pool(qlim)
        self.qlim = qlim

    def bad_classes(self, n):
        """list of dicts for every phi-bad class in Q range at this n, sorted by least-member info"""
        A1 = (n+1)//2; c1 = 1 + (n & 1); L = self.V[2*A1 + c1 - 1]
        out = []
        for Q in range(-self.qlim, self.qlim+1):
            p = abs(Q) % 2
            r = (n - Q) & 1
            A = (n + Q + r)//2
            if A < 0 or A > self.Amax: continue
            if int(self.V[2*A + r]) != int(L):
                k = (Q, p)
                mem = []
                lb = None
                if Q < 0:
                    mem.append((1 << (-Q), 'rule2^Q')); lb = 1 << (-Q)
                if k in self.scan:
                    j = self.scan[k]; mem.append((j, 'scan'))
                    lb = j if lb is None else min(lb, j)
                if k in self.climb: mem.append((self.climb[k][0], 'climb'))
                if k in self.pool: mem.append((self.pool[k][0], 'pool'))
                if lb is None: lb = self.jmax + 1 if self.jmax > 0 else None
                out.append(dict(Q=Q, p=p, r=r, A=A, members=sorted(mem), lb=lb, nocand=not mem))
        out.sort(key=lambda d: (d['lb'] if d['lb'] is not None else 1 << 62))
        return out

    def a_interval(self, n, maxQ=None):
        """(a_lo, a_hi, info) using phi model + direct checks where h>n."""
        bads = self.bad_classes(n)
        a_hi = None; a_lo = None; notes = []
        for d in bads:
            h_of_mem = info(d['members'][0][0])[0] if d['members'] else None
            if d['members'] and h_of_mem <= n:
                cand = d['members'][0][0]
                a_hi = cand - 1 if a_hi is None else min(a_hi, cand - 1)
                continue
            # no member, or member too tall: no valid upper bound here
            if d['lb'] is not None:
                a_lo = d['lb'] - 1 if a_lo is None else min(a_lo, d['lb'] - 1)
        return a_lo, a_hi, bads

if __name__ == '__main__':
    pass
