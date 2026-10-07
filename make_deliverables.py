#!/usr/bin/env python3
"""Final assembly: verify framework outputs against the D1/D2 batch, cross-check published data,
   emit deliverable tables + report.md"""
import csv, glob, sys, os
from collections import Counter, defaultdict
D='/home/user/collatz/data/'
def load_fw():
    fw={}
    for f in ('framework_n0_3953.csv','framework_n3954_13999.csv','framework_n14000_39199.csv'):
        for r in csv.DictReader(open(D+f)): fw[int(r['n'])]=r
    return fw
def load_vf():
    vf={}
    for f in glob.glob('/tmp/vf/out_*.csv'):
        for r in csv.DictReader(open(f)):
            if 'D1' in r: vf[int(r['n'])]=r
    return vf
fw=load_fw(); vf=load_vf()
rows={}
for n,r in fw.items():
    v=vf.get(n)
    a=int(r['a_hi']); lo=int(r['a_lo']); ex=(r['exact']=='yes')
    if v is None: st='PENDING'
    elif v['D1']=='1' and v['D2']=='1': st='EXACT' if ex else 'UBOUND'
    elif v['D1']=='0': st='UNRESOLVED'
    else: st='UNDEREST'
    rows[n]=dict(n=n,a=a,lo=lo,ex=ex,Q=r['best_Q'],src=r['ub_src'],st=st,
                 D1=(v or {}).get('D1',''),D2=(v or {}).get('D2',''),
                 s0=(v or {}).get('s0',''),sa=(v or {}).get('sa',''),sa1=(v or {}).get('sa1',''))
print("status:",dict(Counter(r['st'] for r in rows.values())),file=sys.stderr)

# ---------- 1. results table (n>=3954)
with open(D+'a277109_results_3954_39199.csv','w') as g:
    g.write("n,status,value_or_upper_bound,a_lo,class_Q,minj_source,D1,D2,run_height,height_at_bound,height_at_bound_plus1\n")
    for n in sorted(rows):
        if n<3954: continue
        r=rows[n]
        val = r['a'] if r['st'] in ('EXACT','UBOUND') else (r['a']-1 if r['st']=='UNRESOLVED' else '')
        g.write(f"{n},{r['st']},{val},{r['lo']},{r['Q']},{r['src']},{r['D1']},{r['D2']},{r['s0']},{r['sa']},{r['sa1']}\n")
with open(D+'a277109_new_exact_terms.csv','w') as g:
    g.write("n,a(n)\n")
    for n in sorted(rows):
        if n>=3954 and rows[n]['st']=='EXACT': g.write(f"{n},{rows[n]['a']}\n")

# ---------- 2. b-file check
bf={int(p[0]):int(p[1]) for p in (l.split() for l in open(D+'b277109.txt')) if len(p)==2}
exb=cons=bad=0
for n,v in bf.items():
    r=rows.get(n)
    if not r: continue
    if r['st']=='EXACT':
        if r['a']==v: exb+=1
        else: bad+=1; print("BF MISMATCH",n,v,r['a'],file=sys.stderr)
    else:
        if r['lo']<=v<=r['a']: cons+=1
        else: bad+=1; print("BF NOT-CONTAINED",n,v,r['lo'],r['a'],r['st'],file=sys.stderr)

# ---------- 3. workbook cross-check (H/B + streak claims)
wb={int(r['exponent']):r for r in csv.DictReader(open(D+'collatz_streaks.csv')) if r['exponent'].isdigit()}
Hok=Hbad=Bok=Bnone=0; streak=Counter(); HB_dispute=[]
for e,w in wb.items():
    r=rows.get(e)
    if not r: continue
    H,B=w['H'].strip(),w['B'].strip(); s=w['streak'].strip()
    if H and r['s0']:
        if H==r['s0']: Hok+=1
        else: Hbad+=1; print("H mismatch e=%d H=%s s0=%s"%(e,H,r['s0']),file=sys.stderr)
    if B and r['sa1']:
        if B==r['sa1']: Bok+=1
        else:
            HB_dispute.append((e,int(B),int(r['sa1']),r['st'],r['D1'],r['D2']))
    if s:
        c=int(s) if s.isdigit() else None
        if c is not None:
            if r['st']=='EXACT' and c==r['a']: streak['exact-equal']+=1
            elif c==r['a'] and r['st'] in ('EXACT','UBOUND'): streak['equals-our-certified-value']+=1
            elif c>r['a'] and r['st'] in ('EXACT','UBOUND'): streak['claim-exceeds-our-certified-break']+=1
            else: streak['other']+=1
print("workbook: H match=%d mismatch=%d ; B match=%d dispute=%d"%(Hok,Hbad,Bok,len(HB_dispute)),file=sys.stderr)
print("workbook streak claims:",dict(streak),file=sys.stderr)
with open(D+'frinkle_crosscheck.csv','w') as g:
    g.write("n,workbook_streak_claim,H_run,B_break,our_status,our_value,H_match,B_match\n")
    for e in sorted(wb):
        r=rows.get(e); 
        if not r: continue
        w=wb[e]; H,B,s=w['H'].strip(),w['B'].strip(),w['streak'].strip()
        g.write(f"{e},{s},{H},{B},{r['st']},{r['a']},{str(H==r['s0']) if H else ''},{str(B==r['sa1']) if B else ''}\n")
with open(D+'frinkle_HB_disputes.csv','w') as g:
    g.write("n,workbook_B_break_height,our_break_height,our_status,D1,D2\n")
    for e,B,sa1,st,d1,d2 in sorted(HB_dispute): g.write(f"{e},{B},{sa1},{st},{d1},{d2}\n")

# ---------- 4. chains of a(n)=2^k-1
def chains(statuses):
    seq=[]
    for n in sorted(rows):
        r=rows[n]
        if r['st'] in statuses and r['Q'] and r['Q'].lstrip('-').isdigit():
            Q=int(r['Q'])
            if Q<0 and r['a']==(1<<(-Q))-1: seq.append((n,-Q,r['st']))
    out=[];cur=[]
    for n,k,st in seq:
        if cur and n==cur[-1][0]+1 and k==cur[-1][1]+1: cur.append((n,k,st))
        else:
            if cur: out.append(cur)
            cur=[(n,k,st)]
    if cur: out.append(cur)
    return [c for c in out if len(c)>=2]
ch=chains(('EXACT','UBOUND'))
with open(D+'chains_2k_minus1_certified.csv','w') as g:
    g.write("n_start,n_end,length,k_start,k_end,statuses\n")
    for c in sorted(ch,key=lambda x:-len(x)):
        g.write(f"{c[0][0]},{c[-1][0]},{len(c)},{c[0][1]},{c[-1][1]},{dict(Counter(x[2] for x in c))}\n")
print("chains len>=2:",len(ch),"longest:",max((len(c) for c in ch),default=0),file=sys.stderr)

# ---------- 5. A277684 extension
pub={int(p[0]):int(p[1]) for p in (l.split() for l in open(D+'b277684.txt')) if len(p)==2}   # k -> n
best={}
for n in sorted(rows):
    r=rows[n]
    if r['st'] not in ('EXACT','UBOUND'): continue
    if r['Q'] and r['Q'].lstrip('-').isdigit() and int(r['Q'])<0 and r['a']==(1<<(-int(r['Q'])))-1:
        k=-int(r['Q'])
        if k not in best or n<best[k][0] or (n==best[k][0] and r['st']=='EXACT' and best[k][1]!='EXACT'):
            best[k]=(n,r['st'])
agree=mism=0
with open(D+'a277684_extension.csv','w') as g:
    g.write("k,least_n_our_data,status,published_least_n,match\n")
    for k in sorted(best):
        n,st=best[k]; p=pub.get(k)
        if p is None: mk='NEW'
        elif p==n: mk='PUBLISHED-AGREE'; agree+=1
        elif p>n: mk='OURS-SMALLER(TIGHTENS)'; mism+=1
        else: mk='PUBLISHED-SMALLER'; mism+=1
        g.write(f"{k},{n},{st},{p if p is not None else ''},{mk}\n")
print(f"A277684: agree={agree} mismatch={mism} k-range={min(best)}..{max(best)} new(k>65)={sum(1 for k in best if k>65)}",file=sys.stderr)
for k in sorted(best):
    if pub.get(k) not in (None,best[k][0]): print("   A277684 differs k=%d ours=%s pub=%s"%(k,best[k],pub[k]),file=sys.stderr)

# ---------- 6. Hoft implication
viol=[]; ok=0; over=[]
for n in sorted(rows):
    r=rows[n]
    if not (r['Q'] and r['Q'].lstrip('-').isdigit() and int(r['Q'])<0): continue
    k=-int(r['Q'])
    if r['a']!=(1<<k)-1 or k<2: continue
    p=rows.get(n-1)
    if p is None: continue
    want=(1<<(k-1))-1
    if r['st']=='EXACT':
        if p['st']=='EXACT' and p['a']==want: ok+=1
        elif p['a']<want or (p['st']=='EXACT' and p['a']!=want): over.append((n,k,p['a'],want,p['st']))
        elif p['a']>=want: ok+=1
    else:
        if p['a']<want: over.append((n,k,p['a'],want,p['st']))
with open(D+'hoft_implication_check.csv','w') as g:
    g.write("n,a(n)=2^k-1(k),a(n-1)_value_or_bound,required_2^(k-1)-1,a(n-1)_status\n")
    for n,k,a,w,st in over: g.write(f"{n},{k},{a},{w},{st}\n")
print(f"Hoft: consistent={ok} anomaly-rows(a(n-1) below requirement)={len(over)}",file=sys.stderr)
for x in over[:12]: print("   ",x,file=sys.stderr)

# ---------- 7. summary stats for the report
st=Counter(r['st'] for n,r in rows.items() if n>=3954)
print("FINAL n>=3954:",dict(st),file=sys.stderr)
print("PENDING:",sum(1 for r in rows.values() if r['st']=='PENDING'),file=sys.stderr)
new_exact=[n for n in rows if n>=3954 and rows[n]['st']=='EXACT']
print("new exact terms:",len(new_exact),"max a:",max((rows[n]['a'] for n in new_exact),default=0),file=sys.stderr)
