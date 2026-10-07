/* factcheck.c -- verify a claimed run length L for the row n!+1 (A346775).
 *
 * Prints the heights of n!+1+j for j = 0, 1, L/8, L/4, L/2, 3L/4, L-1, L,
 * so a claimed break at offset L can be checked (H at L must differ from H at 0,
 * and sampled interior heights should match).
 *
 * Build: gcc -O2 -o factcheck factcheck.c -lgmp
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

static unsigned long long steps(mpz_t x) {
    unsigned long long s = 0; mpz_t t; mpz_init(t); mpz_set(t, x);
    while (mpz_cmp_ui(t, 1) != 0) {
        if (mpz_odd_p(t)) { mpz_mul_ui(t, t, 3); mpz_add_ui(t, t, 1); }
        else mpz_fdiv_q_2exp(t, t, 1);
        s++;
    }
    mpz_clear(t);
    return s;
}

int main(int argc, char **argv) {
    long n = atol(argv[1]);
    unsigned long long L = strtoull(argv[2], 0, 10);
    mpz_t base, x; mpz_init(base); mpz_init(x);
    mpz_set_ui(base, 1);
    for (long i = 1; i <= n; i++) mpz_mul_ui(base, base, (unsigned long)i);
    unsigned long long offs[8];
    offs[0] = 0; offs[1] = 1; offs[2] = L/8; offs[3] = L/4; offs[4] = L/2;
    offs[5] = (3*L)/4; offs[6] = L-1; offs[7] = L;
    printf("n=%ld claimedL=%llu", n, L);
    unsigned long long h0 = 0;
    for (int i = 0; i < 8; i++) {
        mpz_add_ui(x, base, offs[i]+1);
        unsigned long long h = steps(x);
        if (i == 0) h0 = h;
        printf(" h[%llu]=%llu%s", offs[i], h, (i>0 && h != h0) ? "*" : "");
    }
    printf("  %s\n", "(* = differs from h[0])");
    mpz_clear(base); mpz_clear(x);
    return 0;
}
