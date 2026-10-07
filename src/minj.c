// minj.c -- scan j=1..Jmax, compute (h(j),o(j)) = (accelerated height, #odd values
// in accelerated orbit), and record for Q = 2o-h and parity p = h mod 2:
//    minj[Q][p] = smallest j with that signature
// Also record the largest Q seen per parity, max h, etc.
// Usage: ./minj Jmax outfile
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define SMALL 24
#define SMALLN (1u<<SMALL)
static uint16_t *hs, *os;

static void build_small(void){
    hs = malloc((size_t)SMALLN*2); os = malloc((size_t)SMALLN*2);
    hs[1]=0; os[1]=0;
    for(uint64_t x=2;x<SMALLN;x++){
        uint64_t y=x; uint32_t h=0,o=0;
        while(y>=x){   /* iterate until below x (table already filled) */
            if(y&1){ y=(3*y+1)>>1; h+=1; o+=1; } else { y>>=1; h+=1; }
        }
        hs[x]=h+hs[y]; os[x]=o+os[y];
    }
}

#define QLO (-4000)
#define QHI (4000)
#define QN  (QHI-QLO+1)
static int64_t minj[QN][2];
static int64_t minjodd[QN][2];
static int16_t minH[2], maxH[2];
static int64_t argminH[2];

int main(int argc,char**argv){
    if(argc<3){fprintf(stderr,"usage: %s Jmax outfile\n",argv[0]);return 1;}
    uint64_t Jmax = strtoull(argv[1],0,10);
    uint64_t Jmin = (argc>3)? strtoull(argv[3],0,10) : 1;
    if(Jmin<1) Jmin=1;
    build_small();
    for(int i=0;i<QN;i++){ minj[i][0]=minj[i][1]=-1; minjodd[i][0]=minjodd[i][1]=-1; }
    maxH[0]=maxH[1]=-30000; minH[0]=minH[1]=30000;
    int64_t maxQ[2]={-100000,-100000}, minQ[2]={100000,100000};
    for(uint64_t j=Jmin;j<=Jmax;j++){
        uint64_t y=j; uint32_t h=0,o=0;
        while(y>=SMALLN){
            if(y&1){ y=(3*y+1)>>1; h+=1; o+=1; } else { y>>=1; h+=1; }
        }
        h+=hs[y]; o+=os[y];
        int Q = (int)2*(int)o - (int)h; int p = h&1;
        if(Q>=QLO&&Q<=QHI){
            if(minj[Q-QLO][p]<0) minj[Q-QLO][p]=(int64_t)j;
            if((j&1)&&minjodd[Q-QLO][p]<0) minjodd[Q-QLO][p]=(int64_t)j;
        }
        if((int)h>maxH[p]){ maxH[p]=(int16_t)h; argminH[p]=(int64_t)j; }
        if((int)h<minH[p]&&j>1) minH[p]=(int16_t)h;
        if(Q>maxQ[p]) maxQ[p]=Q;
        if(Q<minQ[p]) minQ[p]=Q;
    }
    FILE *f=fopen(argv[2],"wb");
    fwrite(&Jmax,8,1,f);
    fwrite(minj,sizeof(int64_t),QN*2,f);
    fwrite(minjodd,sizeof(int64_t),QN*2,f);
    fwrite(maxQ,sizeof(int64_t),2,f); fwrite(minQ,sizeof(int64_t),2,f);
    fwrite(maxH,sizeof(int16_t),2,f); fwrite(minH,sizeof(int16_t),2,f);
    fwrite(argminH,sizeof(int64_t),2,f);
    fclose(f);
    fprintf(stderr,"Jmax=%llu maxQ(odd/even h)=%lld/%lld minQ=%lld/%lld maxH=%d/%d (at %lld/%lld) minH=%d/%d\n",
      (unsigned long long)Jmax,(long long)maxQ[0],(long long)maxQ[1],(long long)minQ[0],(long long)minQ[1],
      maxH[0],maxH[1],(long long)argminH[0],(long long)argminH[1],minH[0],minH[1]);
    return 0;
}
