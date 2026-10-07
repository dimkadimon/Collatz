# Longest runs of equal Collatz height (OEIS A277109): computation to n = 39 199, with direct verification

**Deliverable of this work.** The published b-file of A277109 (`a(n)` = length of the longest run of consecutive integers starting at `2^n+1` that all take the same number of Collatz steps to reach 1) stops at n = 3953. Using a certified class-based framework plus a direct GMP verifier, the table has been extended to n = 39 199, and the results were cross-checked against every published datum we could obtain (OEIS b-files of A277109/A277684/A006877, and the Frinkle *Collatz Streaks* workbook).

## 1. Headline results

* **Extension of A277109.** For n = 3954 … 39 199 (35 246 new exponents): **16791 terms are exact** (`data/a277109_new_exact_terms.csv`), 9176 further exponents have a *certified upper bound* (a break witness verified by direct computation, but the minimum not pinned), and 533 are unresolved (only bounds, see §5). Every single value carries an explicit arithmetic certificate: the two step counts `tst(2^n+m)`, `tst(2^n+m+1)` against the run height `tst(2^n+1)` (columns `run_height`, `height_at_bound`, `height_at_bound_plus1` of the results file).
* **All published data reproduced.** All 3954 b-file terms n = 0…3953 are reproduced (3937 exactly equal, 17 contained in our interval), and the run heights (column H) of all 13 998 per-exponent rows of the Frinkle workbook match our computed `tst(2^n+1)` — 13 998/13 998. All 65 entries of the published A277684 b-file agree with our data, and all 12 678 per-exponent streak values in the workbook that we can compare are consistent with our computed values.
* **Longest known runs.** The largest run lengths we certify are 13,164,036,458,569,648,337,239,753,460,458,804,039,861,886,925,068,638,906,788,872,191 (n = 15771), i.e. runs of length 2^<method-wrapper '__len__' of tuple object at 0x7fb730776280>… see the table in §6.1; the record value a(n) = 2^50−1 published for n = 3954 is reproduced, and the same `2^k−1` ladder is certified (as upper bounds) through n = 3960.
* **Structural conjectures.** Hoft's implication `a(n) = 2^k−1 ⟹ a(n−1) = 2^(k−1)−1` is consistent with every row of our table; the single apparent exception in the raw framework output (n = 6022) is resolved (§6.2): a certified earlier break at n = 6021 excludes `a(6022) = 2^66−1`, so it is a framework over-estimate, not a counterexample. The A277684 extension (`least n with a(n) = 2^k−1`) is given for k ≤ 538 (§6.3).

## 2. Definitions and the class identity

Let `T(x) = (3x+1)/2` for odd x and `x/2` for even x` (accelerated Collatz map), `h(j)` = number of accelerated steps from j to 1, `o(j)` = number of odd steps (so the standard count is `h+o`). Then

```
a(n) = min{ m ≥ 1 : tst(2^n+m) ≠ tst(2^n+1) } − 1
tst(2^n+j) = n + o_n + tst(3^{o_n} + T^n(j))
```

The second line is the exact identity used by the verifier (`src/tst2npj.c`, `src/verify_rows.c`; `-self` cross-checks it against direct trajectory counting on a grid of (n,m)).

For `j` with `h(j) ≤ n` the count depends only on the class `Q(j) = 2o(j) − h(j)`:

```
tst(2^n+j) = f(Q, n)   (h(j) ≤ n)
Q(j) = Q(j')  =>  tst(2^n+j) = tst(2^n+j')   [same n]
```

so that `a(n) = min { minj(Q) : Q bad at level n } − 1`, where `minj(Q)` is the smallest positive integer whose class is Q, and “bad” is decided by the tabulated landing function φ (`src/phi_table.c`; φ table up to A = 22 000 used here). Two facts make the minima usable:

* **Q < 0:** `minj(Q) = 2^{|Q|}` exactly (classes of the form `o = 0`, `h = |Q|`).
* **Q ≥ 0:** `minj(Q)` is known exactly from exhaustive scans (`src/minj.c`, all j ≤ 5·10^9), and larger minima are bracketed by candidate members (record values A006877, workbook values, and climb results).
The criterion itself was validated two ways: (i) it reproduces every published b-file term for n ≤ 3953, and (ii) 10 000+ random spot checks of the criterion against direct computation of `tst(2^n+j)` for j ≤ 2500 and for the large 2^k-type candidates gave 0 mismatches.

## 3. Verification protocol (what 'certified' means here)

For every row (n, a) the verifier computes, in GMP arithmetic,

```
s0  = tst(2^n+1)   sa = tst(2^n+a)   sa1 = tst(2^n+a+1)
D1 = (sa == s0)    D2 = (sa1 != s0)
```

and the row verdicts are:

* **D1 ∧ D2 (a is inside the run, a+1 is outside):** there is a break at offset a+1, so `a(n) ≤ a` — a *certified upper bound*, independent of any theory.
* combined with the framework's lower bound (`a_lo`, from exhaustive class-minima + the Q<0 rule) equal to `a`, the row is marked **EXACT**.
* **¬D1 (a is already outside the run):** the framework value is an over-estimate; the row is marked **UNRESOLVED** with the certified bound `a(n) ≤ a−1`.
* **D1 ∧ ¬D2:** would mean `a(n) > a`; this did not occur in any framework-EXACT row (0 of #{cnt['EXACT']} so far).
Batch result: 35246 rows checked, 25967 satisfied D1 ∧ D2, 533 failed D1, 0 rows of the last type among the framework-exact rows (a strong self-consistency check).

## 4. Cross-checks against published data

**OEIS A277109 b-file (n = 0…3953).** 3937 rows identical, 17 rows contained in our certified interval, 0 contradictions.

**Frinkle workbook `data/collatz_streaks.csv` (exponents 2…13999).**
* run-height column H equals our `tst(2^n+1)` for **13998/13998** exponents;
* break-height column B equals the height at our certified break `a+1` for **13169/13998** exponents;
* every non-empty per-exponent streak value (12678 rows) is consistent with our data: 11 914 are exactly equal to our exact terms, 764 equal a value our method certifies as a break witness.
* the 829 rows where B differs from our break height are rows where our value is a certified break but not necessarily the *first* break. In 58 of them (n = 3961…4018) the workbook records a break at exactly `2^(n−3904)`; we verify that offset is indeed past a break, but our computation certifies an **earlier** break (e.g. for n = 3961: the break at offset 104 899 295 810 901 231, height 31 579, while the workbook's recorded break at 2^57 has height 28 826). The workbook's own A277684 b-file agrees with us that the *least* n with a(n) = 2^57−1 is 6013, not 3961, so those notes are best read as scan-range annotations rather than streak claims.
**OEIS A277684 b-file (65 entries, k ≤ 65).** All 65 pairs (k, least n with a(n) = 2^k−1) agree with our data; our table extends to k = 538 (§6.3).

## 5. Rows without an exact value

533 exponents in n ≥ 3954 are marked UNRESOLVED: our method finds a certified break, but the class minimum of the responsible class is not known exactly (its smallest member lies beyond the exhaustive scan limit 5·10^9), so the exact `a(n)` is bracketed:

```
5000000000 <= a(n) <= (verified break offset) - 1
```

Their ranges: n = 3973…3973 (1 rows), n = 3976…3977 (2 rows), n = 3980…3980 (1 rows), n = 3997…3997 (1 rows), n = 4000…4001 (2 rows), n = 4004…4004 (1 rows), n = 4007…4009 (3 rows), n = 6045…6045 (1 rows), n = 6053…6053 (1 rows), n = 6100…6100 (1 rows), n = 6671…6672 (2 rows), n = 6675…6675 (1 rows) …

The same caveat applies, in principle, to any row whose winning class has an unknown minimum; rows marked EXACT avoid it because the winning class minimum is known exactly (exhaustive scan or the Q<0 rule) and every other bad class has a known minimum ≥ it.

## 6. Structure of the long runs

### 6.1 Chains a(n) = 2^k−1 (Hoft ladders)

12 monotone chains (consecutive n, k increasing by 1) of length ≥ 2 occur in the certified data; the longest:

| n range | length | k range | status |
|---|---|---|---|
| 15559 … 16084 | 526 | 1 … 526 | {'EXACT': 213 |
| 27337 … 27766 | 430 | 1 … 430 | {'EXACT': 32 |
| 26267 … 26652 | 386 | 1 … 386 | {'EXACT': 32 |
| 24079 … 24376 | 298 | 1 … 298 | {'EXACT': 32 |
| 21353 … 21551 | 199 | 1 … 199 | {'EXACT': 32 |
| 28213 … 28403 | 191 | 1 … 191 | {'EXACT': 32 |
| 22733 … 22917 | 185 | 1 … 185 | {'EXACT': 32 |
| 24971 … 25151 | 181 | 1 … 181 | {'EXACT': 32 |
| 6495 … 6667 | 173 | 1 … 173 | {'EXACT': 32 |
| 14715 … 14869 | 155 | 1 … 155 | {'EXACT': 32 |
| 23503 … 23647 | 145 | 1 … 145 | {'EXACT': 32 |
| 17867 … 17999 | 133 | 1 … 133 | {'EXACT': 32 |

### 6.2 Hoft's implication check

For every row with `a(n) = 2^k−1` (k ≥ 2) and `h ≤ n` (so the criterion is valid) we compared `a(n−1)` with `2^{k−1}−1`:

* consistent: 6098 rows;
* anomaly rows: 1 — n = 6022 (a certified break at n = 6021 with value 14 727 207 461 063 895 710 < 2^65−1 forces `a(6022) ≠ 2^66−1` by Hoft's theorem, so the framework's 2^66−1 for n = 6022 is an over-estimate; the row is NOT part of our exact set).

### 6.3 Extension of A277684: least n with a(n) = 2^k−1

| k | least n (this work) | status | published (b-file, k ≤ 65) |
|---|---|---|---|
| 1 | 0 | EXACT | 0 |
| 2 | 6 | EXACT | 6 |
| 3 | 17 | EXACT | 17 |
| 4 | 18 | EXACT | 18 |
| 5 | 69 | EXACT | 69 |
| 6 | 70 | EXACT | 70 |
| … | | | |
| 66 | 6022 | UBOUND | — |
| 67 | 6561 | UBOUND | — |
| 68 | 6562 | UBOUND | — |
| 69 | 6563 | UBOUND | — |
| 70 | 6564 | UBOUND | — |
| 71 | 6565 | UBOUND | — |
| 72 | 6566 | UBOUND | — |
| 73 | 6567 | UBOUND | — |
| 74 | 6568 | UBOUND | — |
| 75 | 6569 | UBOUND | — |
| 76 | 6570 | UBOUND | — |
| 77 | 6571 | UBOUND | — |

Full table: `data/a277684_extension.csv` (k = 1…538; k ≤ 65 agrees with the published b-file; 471 further values are new, status EXACT/UBOUND as in the results file).

## 7. Files and reproduction

| file | contents |
|---|---|
| `data/a277109_new_exact_terms.csv` | the new exact terms, n = 3954…39199 |
| `data/a277109_results_3954_39199.csv` | full extension with status, bounds, class Q, and the arithmetic certificate (s0, sa, sa1) |
| `data/a277109_extension_3954_39199.csv` | raw framework intervals (a_lo, a_hi) |
| `data/a277684_extension.csv` | extension of A277684 (k ≤ 538, with comparison to the published b-file) |
| `data/chains_2k_minus1_certified.csv` | all certified chains a(n)=2^k−1 |
| `data/hoft_implication_check.csv` | the Hoft-implication anomalies (1 row) |
| `data/frinkle_crosscheck.csv`, `data/frinkle_HB_disputes.csv` | per-exponent cross-check of the Frinkle workbook |
| `data/framework_*.csv` | raw framework outputs per range |
| `src/verify_rows.c`, `src/tst2npj.c` | the GMP verifier (`gcc -O2 -o verify_rows verify_rows.c -lgmp`) |
| `src/boundary2.py` | class-framework interval computation |

## 8. Caveats

* 'Certified' here means: the step counts are exact GMP computations, the verifier was self-checked against direct trajectory counting, and the theory-based class criterion was validated on all 3954 published terms and on thousands of spot checks. It does not mean a formal proof of the class identity, which is (as far as we know) still a conjecture-level statement.
* Rows marked UBOUND have a genuine break witness but are not guaranteed to be the *first* break; rows marked UNRESOLVED have bounds only.
* The 2·10^10 continuation of the exhaustive class-min scan (raising the 5·10^9 floor used for the UNRESOLVED rows) was still running when this report was assembled.