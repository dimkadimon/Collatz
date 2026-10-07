/* factverify.c -- verify / compute run lengths for rows n!+1 (A346775).
 *
 * Reads "n hint" lines (hint = claimed run length, 0 if unknown).
 * For each row:
 *   - if hint > 0 and hint <= budget:  brute force the exact run length
 *     (status EXACT / "BROKEN-AT-k" if the break comes before hint)
 *   - if hint > budget: check h[1..K]==h[0] for K=min(200,budget) and h[hint]!=h[0];
 *     status: HINT-CONSISTENT (claim not refuted) or HINT-REFUTED
 *   - if hint == 0: brute force within budget; EXACT if break found, else LOWERBOUND
 *
 * Build: gcc -O2 -o factverify factverify.c -lgmp
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

static unsigned long long steps(const mpz_t base, unsigned long long off) {
    /* block version: strip all trailing zero bits at once (same total count) */
    mpz_t t; mpz_init(t); mpz_add_ui(t, base, off + 1);
    unsigned long long s = 0;
    while (mpz_cmp_ui(t, 1) != 0) {
        if (mpz_odd_p(t)) { mpz_mul_ui(t, t, 3); mpz_add_ui(t, t, 1); s++; }
        mp_bitcnt_t z = mpz_scan1(t, 0);
        if (z > 0) { mpz_fdiv_q_2exp(t, t, z); s += z; }
    }
    mpz_clear(t);
    return s;
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r");
    unsigned long long budget = strtoull(argv[2], 0, 10);
    long n, maxn = -1;
    mpz_t base; mpz_init(base); mpz_set_ui(base, 1);
    long cur = 0;
    char line[256];
    while (fgets(line, sizeof line, f)) {
        long nn; unsigned long long hint;
        if (sscanf(line, "%ld %llu", &nn, &hint) != 2) continue;
        if (nn > maxn) maxn = nn;
    }
    fclose(f);
    for (long i = 1; i <= maxn; i++) mpz_mul_ui(base, base, (unsigned long)i); /* grows to maxn!'s factors */
    /* recompute base per row (factorials aren't cheaply shared) */
    f = fopen(argv[1], "r");
    while (fgets(line, sizeof line, f)) {
        long nn; unsigned long long hint;
        if (sscanf(line, "%ld %llu", &nn, &hint) != 2) continue;
        mpz_set_ui(base, 1);
        for (long i = 1; i <= nn; i++) mpz_mul_ui(base, base, (unsigned long)i);
        unsigned long long h0 = steps(base, 0);
        if (hint > 0 && hint > budget) {
            int consistent = 1;
            unsigned long long K = budget < 200 ? budget : 200;
            for (unsigned long long k = 1; k <= K; k++)
                if (steps(base, k) != h0) { consistent = 0; break; }
            if (consistent) {
                unsigned long long hL = steps(base, hint);
                consistent = (hL != h0);
            }
            printf("%ld,%llu,%s\n", nn, hint, consistent ? "HINT-CONSISTENT" : "HINT-REFUTED");
        } else {
            unsigned long long lim = (hint > 0 && hint < budget) ? hint + 1 : budget;
            unsigned long long L = 1;
            while (L < lim) { if (steps(base, L) != h0) break; L++; }
            if (L < lim) {
                if (hint > 0 && L != hint) printf("%ld,%llu,%s\n", nn, L, "MISMATCH-vs-hint");
                else printf("%ld,%llu,EXACT\n", nn, L);
            } else if (hint > 0 && hint <= budget && L >= hint) {
                /* reached hint without break?? -> inconsistent */
                unsigned long long hL = steps(base, hint);
                printf("%ld,%llu,%s\n", nn, L, hL != h0 ? "EXACT" : "NOBREAK-PAST-HINT");
            } else printf("%ld,%llu,LOWERBOUND\n", nn, L);
        }
        fflush(stdout);
    }
    fclose(f);
    mpz_clear(base);
    return 0;
}
