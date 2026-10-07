# Collatz — streaks of consecutive numbers with the same height

Certified computations on runs of consecutive integers that all take the same
number of steps to reach 1 under the Collatz (3x+1) map.  This repository
contains the tooling and the data behind our work on

* **OEIS A277109** — starting from `2^n+1`, the length of the longest run of
  consecutive numbers with the same Collatz height;
* **OEIS A277684** — the least `k` with `A277109(k) = 2^n - 1` (exponent
  ladders);
* **OEIS A346775** — the same run-length question starting from `n!+1`
  (factorial rows), authored by D. Kamenetsky.

and the accompanying paper draft (`paper/main.tex`) comparing our method with
K. Frinkle et al., *Computing streaks of consecutive numbers with the same
Collatz height*, J. Supercomputing **81**:1338 (2025).

## The model

For a row `n`, put `a(n) = A277109(n) - 1`, so `a(n) = (least break) - 1`
where the break of row `n` is the least `j >= 1` with
`H(2^n + j) != H(2^n + 1)` (`H` = total stopping time).

Each `j` carries a *class* `(h,o)`, the number of accelerated steps `T(x) =
(3x+1)/2` (x odd), `x/2` (x even) needed to reach 1 from `j` and the number of
odd values seen, with `Q = 2o - h`.  For `h(j) <= n` the exact identity

```
H(2^n + j) = h + o + 3*floor((n-h)/2) + phi(A, 1+r) - A + 2r ,
r = (n-h) & 1 ,  A = (n + Q + r)//2
```

holds, where `phi` is a precomputed table (the "phi table").  Hence the
break-set of row `n` is the union of the *bad* classes (classes whose `phi`
value differs from that of the class of `j = 1`), and

```
a(n) = min over bad classes Q of minj(Q) - 1
```

where `minj(Q)` is the least member of class `Q`.

**Rigorous certification.** `minj(Q) = 2^|Q|` for `Q < 0` (proved), and the
class minima for all classes with `minj <= 5*10^9` are known exactly from an
exhaustive scan.  Every other class has the valid lower bound `minj > 5*10^9`.
Therefore `a(n) = v` is *certified* whenever the least bad class of row `n`
is one of those (i.e. its exact minimum is known) — the program prints
`a_lo == a_hi` exactly in that situation.  Values whose least break lies
beyond the scan limit are reported as bounds `[a_lo, a_hi]`, never as exact.

## Sources

| file | contents |
|------|----------|
| `src/phi_table.c`   | builds the phi table (mpz, `phi_table AMAX OUT`) |
| `src/minj.c`        | exhaustive class-minimum scan, `minj JMAX OUT` |
| `src/classes.py`    | `Model`: phi table + class minima + break search |
| `src/boundary2.py`  | row-by-row certified intervals `[a_lo,a_hi]` |
| `src/verify_rows.c` | model-free direct verification of `(n, a)` claims via GMP `tst` |
| `src/climb_deep.py` | class-minimum climber (upper bounds for large `Q`) |
| `src/certify*.py`   | A277109 / A277684 certification drivers |
| `src/tst2npj.c`, `src/phi.c`, `src/runlen.c`, `src/chaincheck.py`, … | auxiliary tools |
| `src/factstreak.c`, `src/factverify.c`, `src/factcheck.c` | A346775 (factorial rows) |

Build (needs GMP):

```
gcc -O2 -o phi_table phi_table.c -lgmp     # from src/
gcc -O2 -o minj minj.c -lgmp
gcc -O2 -o verify_rows verify_rows.c -lgmp
./phi_table 22000 /tmp/phi22k.bin
./minj 5000000000 /tmp/minj5e9.bin
```

## Data

* `data/b277109.txt`, `data/b277684.txt` — OEIS b-files (published terms).
* `data/framework_n*.csv` — certified intervals `[a_lo,a_hi]` per row.
* `data/a277109_*`, `data/a277684_*` — extension tables and QA vs the
  published/workbook values; `data/frinkle_*.csv` — cross-checks against
  Frinkle et al.
* `data/collatz_streaks.csv` — streak data extracted from the supplementary
  workbook of Frinkle et al. (figshare), used for cross-checks.

Data files retain the licences of their sources (OEIS, figshare); the MIT
licence below applies to the code.

## Reference

K. Frinkle, S. Rowland, J. Pack, M. Gordon, J. Eskue, S. Egbe, M. Morris,
*Computing streaks of consecutive numbers with the same Collatz height*,
J. Supercomput. **81**:1338 (2025).

## Licence

MIT — see `LICENSE`.
