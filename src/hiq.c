// hiq.c -- search for numbers with large Q = o(j) - e(j) (odd-step surplus) by random
// biased walks in the reverse Collatz tree: from 1, inverse-even y->2y (Q-1) or
// inverse-odd y->(2y-1)/3 when y%3==2 (Q+1).  Records the smallest j found per Q.
// usage: ./hiq nwalks capbits targetQmax out.txt [seed]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef unsigned __int128 u128;
static inline int odd_inv_ok(u128 y){ return (int)(y % 3) == 2; }

static u128 best[4096];       // best[Q] = smallest j found with that Q
static int have[4096];

int main(int argc,char**argv){
    long nwalks = argc>1? atol(argv[1]) : 200000;
    int capbits = argc>2? atoi(argv[2]) : 80;
    int Qmax = argc>3? atoi(argv[3]) : 600;
    const char*outf = argc>4? argv[4] : "/tmp/hiq.txt";
    unsigned seed = argc>5? (unsigned)atoi(argv[5]) : 12345u;
    srand(seed);
    memset(have,0,sizeof have);
    u128 cap = ((u128)1)<<capbits;
    for(long w=0; w<nwalks; w++){
        double p = 0.5 + 0.5*((double)(rand()%1000)/1000.0);   // bias to take O steps
        u128 y=1; int o=0,e=0;
        for(int step=0; step<4000; step++){
            int tookO=0;
            if(odd_inv_ok(y)){
                u128 z=(2*y-1)/3;
                if(z!=1 && (double)rand()/RAND_MAX < p){ y=z; o++; tookO=1; }
            }
            if(!tookO){ y*=2; e++; }
            if(y>=cap) break;
            int Q=o-e;
            if(Q>=1 && Q<=Qmax && (!have[Q] || y<best[Q])){ best[Q]=y; have[Q]=1; }
        }
    }
    FILE*f=fopen(outf,"w");
    for(int Q=1;Q<=Qmax;Q++) if(have[Q]){
        char buf[64]; int n=0; u128 v=best[Q];
        if(v==0) buf[n++]='0';
        while(v){ buf[n++]=(char)('0'+(int)(v%10)); v/=10; }
        for(int i=0;i<n/2;i++){ char t=buf[i]; buf[i]=buf[n-1-i]; buf[n-1-i]=t; }
        buf[n]=0;
        fprintf(f,"%d %s\n",Q,buf);
    }
    fclose(f);
    fprintf(stderr,"wrote %s\n",outf);
    return 0;
}
