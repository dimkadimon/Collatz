// runlen.c -- compute A277109(n) for n in [nlo,nhi].
// Method: exact identity  tst(2^n+j) = n + o_j + tst(3^{o_j} + T^n(j)),
// plus, for h(j) <= n, closed form in terms of (h(j),o(j));  j is "good" iff
//    phi((n+Q+r)/2, 1+r) == phi(ceil(n/2), 1+n%2),  Q=2o-h, r=(n-h)&1.
// phi(A,c)=A+tst(3^A+c) from table; smallest j per (Q,h%2)-class from table;
// far-field candidates j=m*2^k (m odd) handled by a precomputed table.
// usage: ./runlen phi.bin minj.bin nlo nhi out.csv [Mmax]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int64_t *phi1,*phi2; static long long AMAX=-1;
static int64_t *MJ; static uint64_t JMAX=0;
#define QSCAN 1600
#define QR 4000
static inline int phi_get(long long A,int c){ if(A<0||A>AMAX) return -1; return c==1?phi1[A]:phi2[A]; }

#define SMALL 24
#define SMALLN (1u<<SMALL)
static uint16_t *hs,*os;

static int64_t farMin[2*QSCAN+1][2];   /* min j = m*2^k > JMAX with signature (Q,p) */

int main(int argc,char**argv){
    if(argc<6){fprintf(stderr,"usage: %s phi.bin minj.bin nlo nhi out.csv [Mmax]\n",argv[0]);return 1;}
    /* --- load phi --- */
    FILE *f=fopen(argv[1],"rb"); if(!f){perror("phi");return 1;}
    { int32_t a; int8_t c; int64_t p; long long mn=1LL<<60,mx=-1,nA=0;
      while(fread(&a,4,1,f)==1){ fread(&c,1,1,f); fread(&p,8,1,f); if(a<mn)mn=a; if(a>mx)mx=a; nA++; }
      AMAX=mx; fclose(f);
      f=fopen(argv[1],"rb"); phi1=calloc((size_t)AMAX+1,8); phi2=calloc((size_t)AMAX+1,8);
      for(long long i=0;i<nA;i++){ fread(&a,4,1,f); fread(&c,1,1,f); fread(&p,8,1,f);
          if(c==1) phi1[a]=p; else phi2[a]=p; }
      fclose(f);
      fprintf(stderr,"phi table: A=0..%lld\n",AMAX);
    }
    /* --- load minj --- */
    f=fopen(argv[2],"rb"); if(!f){perror("minj");return 1;}
    fread(&JMAX,8,1,f);
    MJ=malloc((size_t)8001*2*8);
    { int64_t *tmp=malloc((size_t)8001*2*8);
      fread(tmp,8,8001*2,f); memcpy(MJ,tmp,8001*2*8);
      fread(tmp,8,8001*2,f); /* minjodd, unused */ free(tmp); }
    fclose(f);
    fprintf(stderr,"minj table: Jmax=%llu\n",(unsigned long long)JMAX);
    /* --- small tables --- */
    hs=malloc((size_t)SMALLN*2); os=malloc((size_t)SMALLN*2);
    hs[1]=0;os[1]=0;
    for(uint64_t x=2;x<SMALLN;x++){ uint64_t y=x; uint32_t h=0,o=0;
        while(y>=x){ if(y&1){y=(3*y+1)>>1;h++;o++;} else {y>>=1;h++;} }
        hs[x]=h+hs[y]; os[x]=o+os[y]; }
    /* --- far field table --- */
    long long MMAX = argc>6? atoll(argv[6]) : 1000000;
    for(int i=0;i<=2*QSCAN;i++){ farMin[i][0]=farMin[i][1]=0; }
    for(long long m=1;m<=MMAX;m+=2){
        int hm,om;
        if(m<SMALLN){ hm=hs[m]; om=os[m]; } else { uint64_t y=m; uint32_t h=0,o=0;
            while(y>=SMALLN){ if(y&1){y=(3*y+1)>>1;h++;o++;} else {y>>=1;h++;} }
            hm=h+hs[y]; om=o+os[y]; }
        int Qm=2*om-hm;
        for(int k=0;k<=120;k++){
            unsigned __int128 j128=((unsigned __int128)m)<<k;
            if(j128<=JMAX) continue;
            if(j128 > (unsigned __int128)4e18) break;
            int64_t j=(int64_t)j128;
            int Q=Qm-k; if(Q<-QSCAN||Q>QSCAN) continue;
            int p=(hm+k)&1;
            if(farMin[Q+QSCAN][p]==0 || j<farMin[Q+QSCAN][p]) farMin[Q+QSCAN][p]=j;
        }
    }
    long long nlo=atoll(argv[3]), nhi=atoll(argv[4]);
    FILE *o=fopen(argv[5],"w"); fprintf(o,"n,a_n,first_bad_j,note\n");
    char *gq[2]; gq[0]=calloc(2*QSCAN+3,1); gq[1]=calloc(2*QSCAN+3,1);
    for(long long n=nlo;n<=nhi;n++){
        long long A1=(n+1)/2; int c1=1+(int)(n&1);
        int64_t L=phi_get(A1,c1);
        if(L<0){ fprintf(o,"%lld,NA,NA,phi table too small\n",n); continue; }
        memset(gq[0],0,2*QSCAN+3); memset(gq[1],0,2*QSCAN+3);
        int loQ[2]={QSCAN+1,QSCAN+1}, hiQ[2]={-QSCAN-1,-QSCAN-1}; int contig[2]={1,1}; int nseg[2]={0,0};
        int prevgood[2]={0,0};
        for(int p=0;p<2;p++){
            int r=((n-p)&1); int cc=1+r;
            for(int Q=-QSCAN;Q<=QSCAN;Q++){
                if(((Q&1)+2)%2 != p) continue;
                long long A=(n+Q+r)/2;
                int g = (phi_get(A,cc)==L);
                if(g){ gq[p][Q+QSCAN]=1;
                    if(Q<loQ[p]) loQ[p]=Q; if(Q>hiQ[p]) hiQ[p]=Q;
                    if(!prevgood[p]){ nseg[p]++; } prevgood[p]=1;
                } else if(prevgood[p] && Q<hiQ[p]) contig[p]=0;   /* gap inside range */
            }
        }
        int64_t best=0; const char *why="";
        for(int p=0;p<2;p++){
            int r=((n-p)&1);
            for(int Q=-QSCAN;Q<=QSCAN;Q++){
                if(((Q&1)+2)%2 != p) continue;
                if(Q>=-QR && Q<=QR && MJ[(Q+4000)*2+p]>0){
                    long long A=(n+Q+r)/2;
                    if(!(A>=0&&A<=AMAX)) continue;
                    if(gq[p][Q+QSCAN]) continue;
                    int64_t j=MJ[(Q+4000)*2+p];
                    if(best==0||j<best){ best=j; why="scan"; }
                } else if(gq[p][Q+QSCAN]==0){
                    long long A=(n+Q+r)/2;
                    if(!(A>=0&&A<=AMAX)) continue;
                    int64_t j=farMin[Q+QSCAN][p];
                    if(j>0 && (best==0||j<best)){ best=j; why=(Q<=0)?"farfield(2^k-ish)":"farfield(hi-Q)"; }
                }
            }
        }
        if(best==0) fprintf(o,"%lld,NA,NA,no violator\n",n);
        else fprintf(o,"%lld,%lld,%lld,%s[goodQ p0:%d..%d p1:%d..%d,segs %d/%d]\n",
                     n,best-1,best,why,loQ[0],hiQ[0],loQ[1],hiQ[1],nseg[0],nseg[1]);
    }
    fclose(o); return 0;
}
