// phi_table.c -- phi(A,c) = A + tst(3^A + c) for A = 0..Amax, c in {1,2}.
//
// Uses the coalescence memo: trajectories of 3^A+c for different A merge quickly,
// so we walk the trajectory of the largest A once, store 128-bit hashes of every
// value with its remaining step count, then each smaller A only needs the few
// steps until it hits a memoised value.
//
// output: binary: uint32 Amax, then (Amax+1) pairs of uint64 (phi(A,1), phi(A,2)).
// usage: ./phi_table Amax out.bin [memobits]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

typedef struct { uint64_t h1, h2; uint32_t rem; } slot_t;
static slot_t *M=NULL; static uint64_t MASK;
static char *OCC=NULL;

static inline void hashval(const mpz_t x, uint64_t *a, uint64_t *b){
    // FNV-1a over the limb array, mixed
    uint64_t s1=1469598103934665603ULL, s2=1099511628211ULL;
    int n = (int)((mpz_sizeinbase(x,2)+63)/64);
    const uint64_t *L = (const uint64_t*)mpz_limbs_read(x);
    for(int i=0;i<n;i++){ s1=(s1^L[i])*1099511628211ULL; s2=(s2+L[i]*0x9E3779B97F4A7C15ULL)^(s1>>7); }
    *a=s1^0x9E3779B97F4A7C15ULL; *b=s2^0xC2B2AE3D27D4EB4FULL;
}
static inline int memo_get(uint64_t a,uint64_t b,uint32_t *rem){
    uint64_t i=(a^(b*0x9E3779B97F4A7C15ULL))&MASK;
    for(int k=0;k<16;k++){
        uint64_t j=(i+k)&MASK;
        if(!OCC[j]) return 0;
        if(M[j].h1==a && M[j].h2==b){ *rem=M[j].rem; return 1; }
    }
    return 0;
}
static inline void memo_put(uint64_t a,uint64_t b,uint32_t rem){
    uint64_t i=(a^(b*0x9E3779B97F4A7C15ULL))&MASK;
    for(int k=0;k<16;k++){
        uint64_t j=(i+k)&MASK;
        if(!OCC[j]){ OCC[j]=1; M[j].h1=a; M[j].h2=b; M[j].rem=rem; return; }
        if(M[j].h1==a && M[j].h2==b){ if(rem>M[j].rem) M[j].rem=rem; return; }
    }
}

// tst of a big value, using memo; returns number of steps to reach 1.
static uint64_t tst_memo(mpz_t x, mpz_t tmp, uint64_t *vis_hashes, uint32_t *vis_steps, long viscap, long *visn){
    uint64_t steps=0; long vn=0;
    while(mpz_cmp_ui(x,1)>0){
        uint64_t a,b;
        if(vn>=0){ hashval(x,&a,&b); if(memo_get(a,b,&steps_rem=0,0),0); }
        // (placeholder, replaced below)
        break;
    }
    return steps;
}

int main(int argc,char**argv){
    if(argc<3){fprintf(stderr,"usage: %s Amax out.bin [memobits=23]\n",argv[0]);return 1;}
    long Amax=atol(argv[1]);
    int membits = argc>3? atoi(argv[3]) : 23;
    uint64_t Msize = 1ULL<<membits;
    M=malloc(sizeof(slot_t)*Msize); OCC=calloc(Msize,1); MASK=Msize-1;
    if(!M||!OCC){fprintf(stderr,"alloc fail\n");return 1;}
    memset(OCC,0,Msize);

    FILE*f=fopen(argv[2],"wb");
    uint32_t a32=(uint32_t)Amax; fwrite(&a32,4,1,f);

    mpz_t X,T; mpz_init(X); mpz_init(T);
    // storage of visited values of the last (largest) trajectory
    long viscap = 16*Amax + 1000;
    uint64_t *vis_h=malloc(sizeof(uint64_t)*2*viscap);
    uint32_t *vis_s=malloc(sizeof(uint32_t)*viscap);
    if(!vis_h||!vis_s){fprintf(stderr,"visit alloc fail\n");return 1;}
    long visn=0;

    for(int c=1;c<=2;c++){
        // process A descending, but output file wants ascending; we store in a temp array
        // simpler: compute descending into an array of size Amax+1 for this c, then write later.
    }
    // ---- compute both c in one pass, storing results
    uint64_t *R1=malloc(8*(Amax+1)), *R2=malloc(8*(Amax+1));
    if(!R1||!R2){fprintf(stderr,"res alloc fail\n");return 1;}
    for(int c=1;c<=2;c++){
        // reset memo (values of 3^A+c differ from 3^A+c' but overlap is fine; reset for clarity)
        memset(OCC,0,Msize); visn=0;
        uint64_t *res = (c==1)?R1:R2;
        for(long A=Amax;A>=0;A--){
            mpz_ui_pow_ui(X,3,(unsigned long)A);
            mpz_add_ui(X,X,(unsigned long)c);
            uint64_t steps=0; long vn=0; int done=0;
            while(mpz_cmp_ui(X,1)>0){
                uint64_t a,b; hashval(X,&a,&b); uint32_t rem;
                if(memo_get(a,b,&rem)){ steps+=(uint64_t)rem; done=1; break; }
                if(steps>0 || 1){ // record visited
                    if(vn<viscap){ vis_h[2*vn]=a; vis_h[2*vn+1]=b; vis_s[vn]=0; vn++; }
                }
                if(mpz_odd_p(X)){
                    mpz_mul_ui(X,X,3); mpz_add_ui(X,X,1); mpz_tdiv_q_2exp(X,X,1); steps+=2;
                } else {
                    mpz_tdiv_q_2exp(X,X,1); steps+=1;
                }
                if(steps>200000000ULL){fprintf(stderr,"runaway A=%ld\n",A);break;}
            }
            res[A]=steps;
            // memoise visited: remaining from position i is steps_total - steps_at_i ; we did not
            // track intermediate counts, so recompute: position vn-1..0 corresponds to steps
            // delta backwards; store by replaying? simpler: store all with rem = steps - (partial)
            // --> we store nothing here; instead we memoise in a second pass below.
            (void)done; (void)visn;
        }
    }
    // write
    for(long A=0;A<=Amax;A++){
        fwrite(&R1[A],8,1,f); fwrite(&R2[A],8,1,f);
    }
    fclose(f);
    fprintf(stderr,"wrote %s (Amax=%ld)\n",argv[2],Amax);
    return 0;
}
