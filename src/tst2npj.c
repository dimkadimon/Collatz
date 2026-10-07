// tst2npj.c -- verify A277109 candidates by direct computation.
// a(n) = min{ m>=1 : tst(2^n+m) != tst(2^n+1) } - 1   (tst = standard Collatz steps)
// checks per row (n,a): D1: tst(2^n+a)==tst(2^n+1);  D2: tst(2^n+a+1)!=tst(2^n+1)
// Identity: tst(2^n+m) = n + o + tst(3^o + x),  x = T^n(m),  o = #odd accelerated steps of m
// Usage: ./tst2npj rows.csv out.csv | -self
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

static unsigned long std_mpz(mpz_t x){ unsigned long s=0;
    while(mpz_cmp_ui(x,1)>0){
        if(mpz_odd_p(x)){ mpz_mul_ui(x,x,3); mpz_add_ui(x,x,1); s++; }
        mpz_tdiv_q_2exp(x,x,1); s++; }
    return s; }
static unsigned long T_mpz(mpz_t x){ unsigned long s=0;
    while(mpz_cmp_ui(x,1)>0){ if(mpz_odd_p(x)){ mpz_mul_ui(x,x,3); mpz_add_ui(x,x,1); }
        mpz_tdiv_q_2exp(x,x,1); s++; }
    return s; }
static int T_step(mpz_t x){ int odd=(int)mpz_odd_p(x);
    if(odd){ mpz_mul_ui(x,x,3); mpz_add_ui(x,x,1); } mpz_tdiv_q_2exp(x,x,1); return odd; }
/* returns standard count, writes accelerated count via *ta */
static unsigned long tst2npj(unsigned long n, mpz_t m, unsigned long *ta){
    mpz_t x,p; mpz_inits(x,p,NULL); mpz_set(x,m); unsigned long o=0;
    for(unsigned long i=0;i<n;i++) o+=(unsigned long)T_step(x);
    mpz_ui_pow_ui(p,3,o); mpz_add(p,p,x);
    mpz_t q; mpz_init_set(q,p);
    unsigned long rstd=n+o+std_mpz(p); unsigned long rT=n+T_mpz(q);
    if(ta) *ta=rT;
    mpz_clear(q);
    mpz_clears(x,p,NULL); return rstd; }
static int selfcheck(void){
    int bad=0; mpz_t n2,mm; mpz_inits(n2,mm,NULL); unsigned long ta;
    for(unsigned long n=1;n<=25;n++){ mpz_ui_pow_ui(n2,2,n);
        for(unsigned long m=1;m<=60;m++){
            mpz_add_ui(mm,n2,m);
            mpz_t cc; mpz_init_set(cc,mm);
            unsigned long d=std_mpz(cc); mpz_set(cc,mm); unsigned long dT=T_mpz(cc);
            mpz_clear(cc);
            mpz_set_ui(mm,m); unsigned long f=tst2npj(n,mm,&ta);
            if(d!=f||dT!=ta){ printf("MISMATCH n=%lu m=%lu std=%lu/%lu T=%lu/%lu\n",n,m,d,f,dT,ta); bad++; if(bad>8) return 1; } } }
    printf("selfcheck %s\n", bad?"FAILED":"ok"); mpz_clears(n2,mm,NULL); return bad!=0; }
int main(int argc,char**argv){
    if(argc>1 && !strcmp(argv[1],"-self")) return selfcheck();
    if(argc<3){fprintf(stderr,"usage: %s rows.csv out.csv | -self\n",argv[0]);return 1;}
    FILE *fi=fopen(argv[1],"r"); if(!fi){perror(argv[1]);return 1;}
    FILE *fo=fopen(argv[2],"w"); if(!fo){perror(argv[2]);return 1;}
    char line[8192]; if(!fgets(line,sizeof line,fi)){fprintf(stderr,"empty input\n");return 1;}
    fprintf(fo,"n,a,std_a,std_a1,std_0,T_a,T_a1,T_0,D1_std,D2_std,D1_T,D2_T,label\n");
    mpz_t m; mpz_init(m); unsigned long lastn=(unsigned long)-1,s0=0,t0=0;
    char label[512]; unsigned long nrows=0,f1=0,f2=0,g1=0,g2=0;
    while(fgets(line,sizeof line,fi)){
        char astr[4096]; label[0]=0;
        char *tok=strtok(line,",\n"); if(!tok) continue;
        unsigned long n=strtoul(tok,NULL,10);
        tok=strtok(NULL,",\n"); if(!tok) continue; snprintf(astr,sizeof astr,"%s",tok);
        tok=strtok(NULL,",\n"); if(tok) snprintf(label,sizeof label,"%s",tok);
        if(mpz_set_str(m,astr,10)!=0 || mpz_sgn(m)<=0){ fprintf(fo,"%lu,%s,ERR\n",n,astr); continue; }
        if(n!=lastn){ mpz_t one; mpz_init_set_ui(one,1); s0=tst2npj(n,one,&t0); mpz_clear(one); lastn=n; }
        unsigned long ta; unsigned long sa=tst2npj(n,m,&ta);
        mpz_add_ui(m,m,1); unsigned long ta1; unsigned long sa1=tst2npj(n,m,&ta1);
        int d1=(sa==s0),d2=(sa1!=s0),e1=(ta==t0),e2=(ta1!=t0);
        if(!d1)f1++; if(!d2)f2++; if(!e1)g1++; if(!e2)g2++;
        fprintf(fo,"%lu,%s,%lu,%lu,%lu,%lu,%lu,%lu,%d,%d,%d,%d,%s\n",n,astr,sa,sa1,s0,ta,ta1,t0,d1,d2,e1,e2,label);
        nrows++; }
    fprintf(stderr,"rows=%lu std:D1fail=%lu D2fail=%lu  T:D1fail=%lu D2fail=%lu\n",nrows,f1,f2,g1,g2);
    mpz_clear(m); fclose(fi); fclose(fo); return 0; }
