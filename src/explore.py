def T(x): return (3*x+1)//2 if x & 1 else x//2
def tst(x):
    s=0
    while x!=1:
        x = 3*x+1 if x&1 else x//2; s+=1
    return s

print("== tst(3^A + c), c=1,2 ==")
A1=[];A2=[]
for A in range(0,64):
    a=tst(3**A+1); b=tst(3**A+2); A1.append(a);A2.append(b)
    print(f"A={A:3d}: tst(3^A+1)={a:5d}   tst(3^A+2)={b:5d}   [diff to A+1: {a+1-A:5d},{b+1-A:5d}]")
print("deltas c=1:", [A1[i+1]-A1[i] for i in range(len(A1)-1)])
print("deltas c=2:", [A2[i+1]-A2[i] for i in range(len(A2)-1)])
