/* factstreak.c -- Collatz run lengths starting at n!+1  (A346775)
 *
 * a(n) = length of the longest sequence of consecutive numbers that all take
 * the same number of steps to reach 1, starting from n!+1.
 *
 * Usage: factstreak nlo nhi [maxscan]
 *   prints: n,value,status  (status: EXACT or LOWERBOUND)
 *
 * Build: gcc -O2 -o factstreak factstreak.c -lgmp
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

static unsigned long long steps(mpz_t x) {
    unsigned long long s = 0;
    mpz_t t; mpz_init(t); mpz_set(t, x);
    while (mpz_cmp_ui(t, 1) != 0) {
        if (mpz_odd_p(t)) { mpz_mul_ui(t, t, 3); mpz_add_ui(t, t, 1); }
        else mpz_fdiv_q_2exp(t, t, 1);
        s++;
    }
    mpz_clear(t);
    return s;
}

int main(int argc, char **argv) {
    long nlo = atol(argv[1]), nhi = atol(argv[2]);
    unsigned long long maxscan = (argc > 3) ? strtoull(argv[3], 0, 10) : 2000000ULL;
    mpz_t base, x;
    mpz_init(base); mpz_init(x);
    mpz_set_ui(base, 1);
    for (long n = 0; n <= nhi; n++) {
        if (n > 0) mpz_mul_ui(base, base, (unsigned long)n);   /* base = n! */
        if (n < nlo) continue;
        mpz_add_ui(x, base, 1);                                 /* x = n!+1 */
        unsigned long long h0 = steps(x);
        unsigned long long L = 1;
        int exact = 1;
        while (L < maxscan) {
            mpz_add_ui(x, base, L + 1);                         /* n!+1+L */
            if (steps(x) != h0) break;
            L++;
        }
        if (L >= maxscan) exact = 0;
        printf("%ld,%llu,%s\n", n, L, exact ? "EXACT" : "LOWERBOUND");
        fflush(stdout);
    }
    mpz_clear(base); mpz_clear(x);
    return 0;
}
