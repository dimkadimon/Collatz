// phi_table.c -- phi(A,c) = A + tst(3^A + c), A = 0..Amax, c in {1,2}.
//
// Coalescence memo: the trajectories of 3^A+c for different A merge after few steps,
// so the trajectory of the largest A is walked once, every visited value is stored as
// a 128-bit hash + its remaining step count; each smaller A then needs only the short
// path until it hits a memoised value.
//
// output (binary): uint32 Amax, then (Amax+1) x (uint64 phi(A,1), uint64 phi(A,2)).
// usage: ./phi_table Amax out.bin [memobits=22]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

static uint64_t *MH1=NULL,*MH2=NULL; static uint32_t *MR=NULL; static char *OCC=NULL;
static uint64_t MSIZE, MASK;

static inline void hashval(const mpz_t x, uint64_t *a, uint64_t *b){
    uint64_t s1=1469598103934665603ULL, s2=1099511628211ULL;
    size_t n = mpz_size(x);                 // limbs
    const mp_limb_t *L = mpz_limbs_read(x);
    for(size_t i=0;i<n;i++){ s1=(s1^L[i])*1099511628211ULL; s2=(s2+L[i]*0x9E3779B97F4A7C15ULL)^(s1>>7); }
    // include the sign/zero flag implicitly: values are positive
    *a=s1+0x9E3779B97F4A7C15ULL; *b=s2+0xC2B2AE3D27D4EB4FULL;
}
static inline uint64_t slot0(uint64_t a,uint64_t b){ return (a^(b*0x9E3779B97F4A7C15ULL))&MASK; }
static inline int memo_get(uint64_t a,uint64_t b,uint32_t *rem){
    uint64_t i=slot0(a,b);
    for(int k=0;k<24;k++){
        uint64_t j=(i+k)&MASK;
        if(!OCC[j]) return 0;
        if(MH1[j]==a && MH2[j]==b){ *rem=MR[j]; return 1; }
    }
    return 0;
}
static inline void memo_put(uint64_t a,uint64_t b,uint32_t rem){
    uint64_t i=slot0(a,b);
    for(int k=0;k<24;k++){
        uint64_t j=(i+k)&MASK;
        if(!OCC[j]){ OCC[j]=1; MH1[j]=a; MH2[j]=b; MR[j]=rem; return; }
        if(MH1[j]==a && MH2[j]==b){ if(rem>MR[j]) MR[j]=rem; return; }
    }
}

int main(int argc,char**argv){
    if(argc<3){fprintf(stderr,"usage: %s Amax out.bin [memobits=22]\n",argv[0]);return 1;}
    long Amax=atol(argv[1]);
    int membits = argc>3? atoi(argv[3]) : 22;
    MSIZE=1ULL<<membits; MASK=MSIZE-1;
    OCC=calloc(MSIZE,1); MH1=malloc(8*MSIZE); MH2=malloc(8*MSIZE); MR=malloc(4*MSIZE);
    if(!OCC||!MH1||!MH2||!MR){fprintf(stderr,"memo alloc fail\n");return 1;}
    uint64_t *R1=malloc(8*(Amax+1)), *R2=malloc(8*(Amax+1));
    if(!R1||!R2){fprintf(stderr,"res alloc fail\n");return 1;}
    long VISMAX = 9*Amax+100000;
    uint64_t *VH1=malloc(8*VISMAX), *VH2=malloc(8*VISMAX); uint64_t *VS=malloc(8*VISMAX);
    if(!VH1||!VH2||!VS){fprintf(stderr,"vis alloc fail\n");return 1;}

    mpz_t X,Y; mpz_init(X); mpz_init(Y);
    for(int c=1;c<=2;c++){
        memset(OCC,0,MSIZE);
        uint64_t *res = (c==1)?R1:R2;
        mpz_ui_pow_ui(Y,3,(unsigned long)Amax);   // Y = 3^A, decremented incrementally
        unsigned long long totsteps=0; long maxvn=0;
        for(long A=Amax;A>=0;A--){
            if(A<Amax) mpz_divexact_ui(Y,Y,3);
            mpz_add_ui(X,Y,(unsigned long)c);
            uint64_t steps=0; long vn=0;
            while(mpz_cmp_ui(X,1)>0){
                uint64_t a,b; hashval(X,&a,&b);
                uint32_t rem;
                if(memo_get(a,b,&rem)){ steps+=(uint64_t)rem; break; }
                if(vn<VISMAX){ VH1[vn]=a; VH2[vn]=b; VS[vn]=steps; vn++; }
                if(mpz_odd_p(X)){ mpz_mul_ui(X,X,3); mpz_add_ui(X,X,1); mpz_tdiv_q_2exp(X,X,1); steps+=2; }
                else { mpz_tdiv_q_2exp(X,X,1); steps+=1; }
            }
            res[A]=steps+(uint64_t)A;   // phi = A + tst(3^A+c)
            totsteps+=vn; if(vn>maxvn) maxvn=vn;
            if(getenv("VTRACE") && Amax<=3000) fprintf(stderr,"    A=%ld vn=%ld\n",A,vn);
            for(long i=vn-1;i>=0;i--) memo_put(VH1[i],VH2[i],(uint32_t)(steps-VS[i]));
            if((Amax-A)%20000==0) fprintf(stderr,"  c=%d A=%ld steps=%llu memo_vn=%ld\n",c,A,(unsigned long long)steps,vn);
        }
        fprintf(stderr,"finished c=%d  total_visited=%llu  max_path=%ld\n",c,totsteps,maxvn);
    }
    FILE*f=fopen(argv[2],"wb");
    uint32_t a32=(uint32_t)Amax; fwrite(&a32,4,1,f);
    for(long A=0;A<=Amax;A++){ fwrite(&R1[A],8,1,f); fwrite(&R2[A],8,1,f); }
    fclose(f);
    fprintf(stderr,"wrote %s (Amax=%ld)\n",argv[2],Amax);
    return 0;
}
