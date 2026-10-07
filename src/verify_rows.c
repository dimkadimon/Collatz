// verify_rows.c -- batch D1/D2 verification for A277109 rows (n,a).
// a(n) = min{m>=1 : tst(2^n+m) != tst(2^n+1)} - 1.
// D1: tst(2^n+a)  == tst(2^n+1)      D2: tst(2^n+a+1) != tst(2^n+1)
// D1&&D2  =>  a(n) == a  (proved by direct computation).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
static unsigned long std_mpz(mpz_t x){ unsigned long s=0;
    while(mpz_cmp_ui(x,1)>0){ if(mpz_odd_p(x)){ mpz_mul_ui(x,x,3); mpz_add_ui(x,x,1); s++; } mpz_tdiv_q_2exp(x,x,1); s++; } return s; }
static int T_step(mpz_t x){ int odd=(int)mpz_odd_p(x); if(odd){ mpz_mul_ui(x,x,3); mpz_add_ui(x,x,1);} mpz_tdiv_q_2exp(x,x,1); return odd; }
static unsigned long tst2npj(unsigned long n, mpz_t m){ mpz_t x,p; mpz_inits(x,p,NULL); mpz_set(x,m); unsigned long o=0;
    for(unsigned long i=0;i<n;i++) o+=(unsigned long)T_step(x);
    mpz_ui_pow_ui(p,3,o); mpz_add(p,p,x);
    unsigned long r=n+o+std_mpz(p); mpz_clears(x,p,NULL); return r; }
int main(int argc,char**argv){
    if(argc<3){ fprintf(stderr,"usage: %s rows.csv out.csv\n",argv[0]); return 1; }
    FILE*fi=fopen(argv[1],"r"); if(!fi){perror("in");return 1;}
    FILE*fo=fopen(argv[2],"w"); if(!fo){perror("out");return 1;}
    fprintf(fo,"n,a,s0,sa,sa1,D1,D2,label\n");
    char line[65536]; if(!fgets(line,sizeof line,fi)){return 1;}
    mpz_t a,a1; mpz_inits(a,a1,NULL);
    unsigned long lastn=(unsigned long)-1, s0=0, nrows=0,f1=0,f2=0;
    while(fgets(line,sizeof line,fi)){
        char *tok=strtok(line,",\n"); if(!tok) continue; unsigned long n=strtoul(tok,NULL,10);
        tok=strtok(NULL,",\n"); if(!tok) continue; 
        char lab[256]=""; char *t3=strtok(NULL,",\n"); if(t3) snprintf(lab,sizeof lab,"%s",t3);
        if(mpz_set_str(a,tok,10)!=0||mpz_sgn(a)<=0){ fprintf(fo,"%lu,ERR\n",n); continue; }
        if(n!=lastn){ mpz_t one; mpz_init_set_ui(one,1); s0=tst2npj(n,one); mpz_clear(one); lastn=n; }
        unsigned long sa=tst2npj(n,a);
        mpz_set(a1,a); mpz_add_ui(a1,a1,1);
        unsigned long sa1=tst2npj(n,a1);
        int d1=(sa==s0), d2=(sa1!=s0);
        if(!d1)f1++; if(!d2)f2++; nrows++;
        fprintf(fo,"%lu,%s,%lu,%lu,%lu,%d,%d,%s\n",n,tok,s0,sa,sa1,d1,d2,lab);
        fflush(fo);
    }
    fprintf(stderr,"rows=%lu D1fail=%lu D2fail=%lu\n",nrows,f1,f2);
    mpz_clears(a,a1,NULL); fclose(fi); fclose(fo); return 0;
}
