// minj_beam.c -- beam search for small j with given class Q = o(j)-e(j) in accelerated terms.
// Reverse accelerated tree on odd values: from odd y, predecessors are y'=(2^a y-1)/3 for a with
// 2^a y == 1 (mod 3); a odd iff y%3==2, a even iff y%3==1; Q changes by (2-a), height by a.
// Seeds: successors of 1 in that tree: y=(2^a-1)/3 (a even) -> 5,21,85,341,...
// Keeps for each level s the SMALLEST values found (upper bounds on minj).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
#define CAP 512
#define SLO (-40)
#define SHI 800
#define AMAX 9
typedef struct { u128 y; int32_t s,h,depth; } st;
static u128 *tab; static int32_t *tabh; static int32_t cnt[SHI-SLO+1];
static st *next; static int nnext;
static int cmp_y(const void*a,const void*b){ const st*x=a,*y=b; return x->y<y->y?-1:(x->y>y->y?1:0); }
int main(int argc,char**argv){
    int Qmax = argc>1?atoi(argv[1]):340;
    int32_t maxdepth = argc>2?atoi(argv[2]):4000;
    int usebig = argc>3?atoi(argv[3]):1;         // allow a up to AMAX
    tab=malloc(sizeof(u128)*(SHI-SLO+1)*(size_t)CAP); tabh=malloc(sizeof(int32_t)*(SHI-SLO+1)*(size_t)CAP);
    next=malloc(sizeof(st)*(size_t)4*CAP*(SHI-SLO+1));
    memset(cnt,0,sizeof cnt);
    #define ADD(S,Y,H) do{ int _s=(S); if(_s>=SLO&&_s<=SHI){ int *c=&cnt[_s-SLO]; u128 _y=(Y); int i; \
        for(i=0;i<*c;i++) if(tab[(_s-SLO)*CAP+i]==_y) break; \
        if(i==*c && *c<CAP){ tab[(_s-SLO)*CAP+i]=_y; tabh[(_s-SLO)*CAP+i]=(H); (*c)++; } \
        else if(i<*c && (H)<tabh[(_s-SLO)*CAP+i]) tabh[(_s-SLO)*CAP+i]=(H); \
      } }while(0)
    // seeds: y=(2^a-1)/3 for even a>=4
    u128 p4=1;
    for(int a=4;a<=40;a+=2){
        p4<<=2;                       // now 2^a for a=4,6,...
        u128 y=(((u128)1<<a)-1)/3;
        int s=2-a, h=a;
        if(s>=SLO) ADD(s,y,h);
    }
    u128 CAPV=((u128)1)<<126;
    for(int depth=0; depth<maxdepth; depth++){
        nnext=0;
        for(int s=SLO;s<=SHI;s++){
            int c=cnt[s-SLO];
            for(int i=0;i<c;i++){
                u128 y=tab[(s-SLO)*CAP+i]; int h=tabh[(s-SLO)*CAP+i];
                int r=(int)(y%3);
                if(r==0) continue;
                int amax = usebig?AMAX:4;
                for(int a= (r==2?1:2); a<=amax; a+=2){
                    u128 num=(((u128)1<<a)*y-1);
                    if(num%3!=0) continue;
                    u128 y2=num/3;
                    if(y2<3||y2>CAPV) continue;
                    int s2=s+2-a;
                    if(s2>SHI) continue;
                    if(nnext < 4*CAP*(SHI-SLO+1)){ next[nnext].y=y2; next[nnext].s=s2; next[nnext].h=h+a; next[nnext].depth=depth; nnext++; }
                }
            }
        }
        if(!nnext) break;
        // merge: sort by y ascending, insert into tables keeping CAP smallest per level
        qsort(next,nnext,sizeof(st),cmp_y);
        for(int i=0;i<nnext;i++){
            int s2=next[i].s;
            u128 y2=next[i].y; int h2=next[i].h;
            int *c=&cnt[s2-SLO];
            u128 *row=&tab[(s2-SLO)*CAP]; int32_t*rh=&tabh[(s2-SLO)*CAP];
            int j;
            for(j=0;j<*c;j++) if(row[j]==y2) break;
            if(j==*c){
                if(*c<CAP){ // append then bubble up
                    row[*c]=y2; rh[*c]=h2; int k=*c; (*c)++;
                    while(k>0 && row[k]<row[k-1]){ u128 t=row[k];row[k]=row[k-1];row[k-1]=t; int32_t th=rh[k];rh[k]=rh[k-1];rh[k-1]=th; k--; }
                }
            } else if(h2<rh[j]) rh[j]=h2;
        }
        if(depth%50==0) fprintf(stderr,"depth=%d states=%d maxQ=%d\n",depth,nnext,SHI);
    }
    // report smallest y per level
    printf("# Q j h\n");
    for(int s=SLO;s<=SHI;s++){
        int c=cnt[s-SLO];
        if(!c) continue;
        if(s<1) continue;
        // smallest y in level
        u128 best=0; int bh=0;
        for(int i=0;i<c;i++){ if(best==0||tab[(s-SLO)*CAP+i]<best){ best=tab[(s-SLO)*CAP+i]; bh=tabh[(s-SLO)*CAP+i]; } }
        char buf[64]; int n=0; u128 v=best; while(v){buf[n++]='0'+(int)(v%10);v/=10;}
        for(int i=0;i<n/2;i++){char t=buf[i];buf[i]=buf[n-1-i];buf[n-1-i]=t;} buf[n]=0;
        if(s<=Qmax) printf("%d %s %d\n",s,buf,bh);
    }
    return 0;
}
