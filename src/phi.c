// phi.c -- compute phi(A,c) = A + tst(3^A + c), A = 0..Amax, c in {1,2}
// where tst(x) = total stopping time of x under 3x+1 (odd step x->(3x+1)/2 counts 2).
// Uses a coalescence memo: hash(trajectory value) -> remaining steps to 1.
// Usage: ./phi Amax outfile
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

#define SMALL  22
#define SMALLN (1u<<SMALL)
static uint16_t *tst_small;             /* tst for x < 2^SMALL (fits in 16 bits) */

#define HBITS 22
#define HSIZE (1ULL<<HBITS)
typedef struct { uint64_t h1,h2; uint32_t rem; uint32_t used; } E;
static E *H;

static inline void hash128(const mpz_t x, uint64_t *o1, uint64_t *o2){
    size_t n = mpz_size(x);
    const mp_limb_t *lim = (const mp_limb_t*)mpz_limbs_read((mpz_ptr)x);
    uint64_t a = 0x1234567890ABCDEFULL, b = 0xFEDCBA0987654321ULL;
    for(size_t i=0;i<n;i++){
        a = (a ^ lim[i]) * 0x100000001B3ULL;
        b = (b + lim[i]) * 0x9E3779B97F4A7C15ULL + i;
    }
    a += (uint64_t)n * 0x9E3779B97F4A7C15ULL;
    b ^= a>>29;
    a ^= b>>31; a *= 0xBF58476D1CE4E5B9ULL;
    b ^= a>>27; b *= 0x94D049BB133111EBULL;
    *o1 = a; *o2 = b;
}
static int hm_get(uint64_t h1,uint64_t h2,uint32_t *rem){
    uint64_t i = (h1*0x9E3779B97F4A7C15ULL)>>(64-HBITS);
    for(int k=0;k<32;k++, i=(i+1)&(HSIZE-1)){
        E *e=&H[i];
        if(!e->used) return 0;
        if(e->h1==h1&&e->h2==h2){ *rem=e->rem; return 1; }
    }
    return 0;
}
static void hm_put(uint64_t h1,uint64_t h2,uint32_t rem){
    uint64_t i=(h1*0x9E3779B97F4A7C15ULL)>>(64-HBITS);
    for(int k=0;k<32;k++, i=(i+1)&(HSIZE-1)){
        E *e=&H[i];
        if(!e->used){ e->used=1; e->h1=h1; e->h2=h2; e->rem=rem; return; }
        if(e->h1==h1&&e->h2==h2) return;
    }
}

/* simple direct tst for small x */
static uint32_t tst_small_one(uint64_t x){
    uint64_t s=0;
    while(x!=1){
        if(x&1){ x=(3*x+1)>>1; s+=2; } else { x>>=1; s+=1; }
    }
    return (uint32_t)s;
}
static void build_small(void){
    tst_small = malloc((size_t)SMALLN*sizeof(uint16_t));
    for(uint64_t x=1;x<SMALLN;x++) tst_small[x]=(uint16_t)tst_small_one(x);
}

/* tst(x) with memo; buffer (hash,stepoffset) of visited big values */
#define BCAP (1u<<21)
static uint64_t *b_h1,*b_h2,*b_off;

static uint64_t tst_memo(mpz_t x0){
    mpz_t x; mpz_init_set(x,x0);
    uint64_t steps=0; size_t nb=0;
    for(;;){
        if(mpz_cmp_ui(x,1)==0) break;
        if(mpz_sizeinbase(x,2) <= SMALL){
            uint64_t v = mpz_get_ui(x);
            steps += tst_small[v];
            break;
        }
        uint64_t h1,h2; hash128(x,&h1,&h2);
        uint32_t rem;
        if(hm_get(h1,h2,&rem)){
            steps += rem;
            for(size_t i=0;i<nb;i++) hm_put(b_h1[i],b_h2[i],(uint32_t)(steps-b_off[i]));
            mpz_clear(x); return steps;
        }
        if(nb<BCAP){ b_h1[nb]=h1; b_h2[nb]=h2; b_off[nb]=steps; nb++; }
        if(mpz_odd_p(x)){ mpz_mul_ui(x,x,3); mpz_add_ui(x,x,1); mpz_tdiv_q_2exp(x,x,1); steps+=2; }
        else { mpz_tdiv_q_2exp(x,x,1); steps+=1; }
    }
    for(size_t i=0;i<nb;i++) hm_put(b_h1[i],b_h2[i],(uint32_t)(steps-b_off[i]));
    mpz_clear(x);
    return steps;
}

int main(int argc,char**argv){
    if(argc<3){ fprintf(stderr,"usage: %s Amax outfile\n",argv[0]); return 1; }
    long long Amax = atoll(argv[1]);
    FILE *f=fopen(argv[2],"wb"); if(!f){perror("fopen");return 1;}
    build_small();
    b_h1=malloc((size_t)BCAP*8); b_h2=malloc((size_t)BCAP*8); b_off=malloc((size_t)BCAP*8);
    mpz_t x; mpz_init(x);
    for(int c=1;c<=2;c++){
        H = calloc(HSIZE,sizeof(E));
        if(!H){fprintf(stderr,"oom\n");return 1;}
        for(long long A=Amax;A>=0;A--){
            mpz_ui_pow_ui(x,3,(unsigned long)A);
            mpz_add_ui(x,x,(unsigned long)c);
            uint64_t t=tst_memo(x);
            int32_t a32=(int32_t)A; int64_t ph=(int64_t)((uint64_t)A+t);
            fwrite(&a32,4,1,f); fwrite(&c,1,1,f); fwrite(&ph,8,1,f);
            if(((Amax-A)%20000)==0) fprintf(stderr,"c=%d A=%lld phi=%lld\n",c,A,(long long)ph);
        }
        free(H); H=NULL;
    }
    fclose(f); fprintf(stderr,"done\n"); return 0;
}
