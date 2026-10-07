def T(x): return (3*x+1)//2 if x & 1 else x//2
def acc_info(j):
    x=j;h=0;o=0
    while x!=1:
        if x&1:o+=1
        x=T(x);h+=1
    return h,o
def tst(x):
    s=0
    while x!=1:
        x = 3*x+1 if x&1 else x//2; s+=1
    return s
N=64
info={j:acc_info(j) for j in range(1,2000)}
def landing(n,j):
    h,o=info[j]
    if h>n: return None
    A = o + (n-h+1)//2
    c = 1 + ((n-h)&1)
    return A,c,3**A+c
print("j : h(j) o(j) | for n=17,18,19: A,c,tst(L),tst(2^n+j)")
for j in range(1,36):
    h,o=info[j]
    row=[f"{j:3d}: h={h:3d} o={o:3d} tst(j)={h+o:3d} |"]
    for n in (17,18,19,20):
        A,c,L=landing(n,j)
        row.append(f"n={n}: A={A:4d} c={c} tstL={tst(L):4d} tst={n+A+tst(L):4d}")
    print(" ".join(row))
