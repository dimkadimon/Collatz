#!/usr/bin/env python3
"""Assemble final A277109 deliverables from the framework outputs + the D1/D2 verification batch."""
import csv, glob, sys
from collections import Counter, defaultdict

DATA='/home/user/collatz/data/'
FW=('framework_n0_3953.csv','framework_n3954_13999.csv','framework_n14000_39199.csv')

fw={}
for f in FW:
    for r in csv.DictReader(open(DATA+f)): fw[int(r['n'])]=r
vf={}
for f in glob.glob('/tmp/vf/out_*.csv'):
    for r in csv.DictReader(open(f)):
        if 'D1' in r: vf[int(r['n'])]=r
print(f"framework rows={len(fw)} verified rows={len(vf)}",file=sys.stderr)

rows=[]
for n in sorted(fw):
    r=fw[n]; v=vf.get(n)
    a=int(r['a_hi']); lo=int(r['a_lo']); ex=(r['exact']=='yes')
    if v is None:
        st='PENDING'; ub=a; d1=d2=''
    else:
        d1,d2=v['D1'],v['D2']
        if d1=='1' and d2=='1': st='EXACT' if ex else 'UBOUND'; ub=a
        elif d1=='0': st='UNRESOLVED'; ub=a-1
        else: st='UNDEREST'; ub=None
    rows.append(dict(n=n,a_hi=a,a_lo=lo,exact=ex,best_Q=r['best_Q'],ub_src=r['ub_src'],status=st,ub=ub,D1=d1,D2=d2,
                     s0=(v or {}).get('s0',''),sa=(v or {}).get('sa',''),sa1=(v or {}).get('sa1','')))

with open(DATA+'a277109_results_3954_39199.csv','w') as g:
    g.write("n,a_status,value_or_upper_bound,a_lo,class_Q,minj_source,certified_D1,certified_D2,height_run_s0,height_at_a_sa,height_at_a1_sa1\n")
    for r in rows:
        if r['n']<3954: continue
        val = r['a_hi'] if r['status'] in ('EXACT','UBOUND') else (r['ub'] if r['ub'] else '')
        g.write(f"{r['n']},{r['status']},{val},{r['a_lo']},{r['best_Q']},{r['ub_src']},{r['D1']},{r['D2']},{r['s0']},{r['sa']},{r['sa1']}\n")

with open(DATA+'a277109_new_exact_terms.csv','w') as g:
    g.write("n,a(n)\n")
    for r in rows:
        if r['n']>=3954 and r['status']=='EXACT': g.write(f"{r['n']},{r['a_hi']}\n")

cnt=Counter((r['status'],'n<3954' if r['n']<3954 else 'n>=3954') for r in rows)
print("STATUS COUNTS:",dict(cnt),file=sys.stderr)

# --- b-file cross-check (n<=3953)
bf={}
for line in open(DATA+'b277109.txt'):
    p=line.split()
    if len(p)==2: bf[int(p[0])]=int(p[1])
ag=ex=cons=bad=0
for n,v in bf.items():
    r=next((x for x in rows if x['n']==n),None)
    if not r: continue
    if r['status']=='EXACT':
        if r['a_hi']==v: ex+=1
        else: bad+=1; print("BF MISMATCH",n,v,r['a_hi'],file=sys.stderr)
    else:
        if r['a_lo']<=v<=r['a_hi']: cons+=1
        else: bad+=1; print("BF NOT-CONTAINED",n,v,r['a_lo'],r['a_hi'],r['status'],file=sys.stderr)
print(f"b-file: exact-match={ex} contained={cons} bad={bad}",file=sys.stderr)

# --- workbook cross-check
def claim(r):
    st=r['streak'].strip()
    if st:
        try: return int(st),'streak'
        except Exception: pass
    return None,''
wb=[r for r in csv.DictReader(open(DATA+'collatz_streaks.csv')) if r['exponent'].isdigit()]
res=Counter(); refs=[]
with open(DATA+'frinkle_crosscheck.csv','w') as g:
    g.write("n,workbook_claim,claim_source,our_status,our_value,delta_vs_bound,evidence\n")
    for r in wb:
        n=int(r['exponent']); c,src=claim(r)
        if c is None or n not in (x['n'] for x in rows): continue
        x=next(y for y in rows if y['n']==n)
        if x['status'] in ('EXACT','UBOUND'):
            if c==x['a_hi'] and x['status']=='EXACT': k='CONFIRMED-EXACT'
            elif c==x['a_hi']: k='CONFIRMED-VALUE'
            elif c>x['a_hi']: k='OURS-SMALLER'   # certified break before their claim
            else: k='WORKBOOK-SMALLER'
        else: k='NO-COMPARE'
        res[(k,src)]+=1
        if c>x['ub'] if x['ub'] else False: refs.append((n,c,src,x['a_hi'],x['status']))
        g.write(f"{n},{c},{src},{x['status']},{x['a_hi']},{c-x['a_hi'] if x['a_hi'] else ''},{x['status']}\n")
print("WORKBOOK:",dict(res),file=sys.stderr)
open(DATA+'frinkle_claim_exceeds_our_bound.csv','w').write("n,claim,source,our_bound,status\n"+"\n".join(",".join(map(str,r)) for r in sorted(refs)))

# --- chains: consecutive n with a(n)=2^k-1 (certified rows first)
def chain_rows(statuses):
    seq=[]
    for r in rows:
        if r['status'] in statuses and r['best_Q'] and r['best_Q'].lstrip('-').isdigit():
            Q=int(r['best_Q'])
            if Q<0 and r['a_hi']==(1<<(-Q))-1: seq.append((r['n'],-Q,r['status']))
    chains=[];cur=[]
    for n,k,st in seq:
        if cur and n==cur[-1][0]+1 and k==cur[-1][1]+1: cur.append((n,k,st))
        else:
            if cur: chains.append(cur)
            cur=[(n,k,st)]
    if cur: chains.append(cur)
    return [c for c in chains if len(c)>=2]

ch=chain_rows(('EXACT','UBOUND'))
with open(DATA+'chains_2k_minus1_certified.csv','w') as g:
    g.write("n_start,n_end,length,k_start,k_end,status_mix\n")
    for c in sorted(ch,key=lambda c:-len(c)):
        st=Counter(x[2] for x in c)
        g.write(f"{c[0][0]},{c[-1][0]},{len(c)},{c[0][1]},{c[-1][1]},{dict(st)}\n")

# --- A277684 extension: least n with a(n)=2^k-1 (certified status required)
best={}
for r in rows:
    if r['status'] not in ('EXACT','UBOUND'): continue
    Q=r['best_Q']
    if not Q or not Q.lstrip('-').isdigit(): continue
    Q=int(Q)
    if Q<0 and r['a_hi']==(1<<(-Q))-1:
        k=-Q
        if k not in best or r['n']<best[k][0] or (r['n']==best[k][0] and r['status']=='EXACT' and best[k][1]!='EXACT'):
            best[k]=(r['n'],r['status'])
with open(DATA+'a277684_extension.csv','w') as g:
    g.write("k (a(n)=2^k-1),least_n_certified,status\n")
    for k in sorted(best): g.write(f"{k},{best[k][0]},{best[k][1]}\n")

# --- Hoft implication check over certified rows
byN={r['n']:r for r in rows}
viol=[]; ok=0; unk=0
for r in rows:
    if r['best_Q'] and r['best_Q'].lstrip('-').isdigit() and int(r['best_Q'])<0 and r['a_hi']==(1<<(-int(r['best_Q'])))-1:
        n=r['n']; prev=byN.get(n-1)
        if prev is None: continue
        k=-int(r['best_Q'])
        want=(1<<(k-1))-1
        if prev['status']=='EXACT':
            if prev['a_hi']==want: ok+=1
            else: viol.append((n,k,prev['a_hi'],want))
        else:
            if prev['a_hi']<want: viol.append((n,k,prev['a_hi'],want))
            else: unk+=1
with open(DATA+'hoft_implication_check.csv','w') as g:
    g.write("n,k,a(n-1)_value_or_bound,required_2^(k-1)-1,verdict\n")
    for n,k,av,want in viol: g.write(f"{n},{k},{av},{want},VIOLATION-OR-OVERESTIMATE\n")
print(f"hoft: ok={ok} unresolved={unk} potential-violations={len(viol)}",file=sys.stderr)
if viol: print("   violations:",viol[:10],file=sys.stderr)

summ=Counter(r['status'] for r in rows)
print("SUMMARY:",dict(summ),file=sys.stderr)
