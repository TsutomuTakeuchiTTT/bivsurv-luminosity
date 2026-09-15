"""Exact linear algebra and certified natural-log enclosures (standard library)."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb

def rat(x):
    if isinstance(x,float): raise TypeError('Supply exact integers/strings/Fraction')
    return F(x)

def down(x,D=10**30): return F(x.numerator*D//x.denominator,D)
def up(x,D=10**30): return -down(-x,D)
@lru_cache(None)
def nearlog(x):
    assert 1<=x<=2
    t=(x-1)/(x+1); z=t;s=F(0)
    for j in range(40): s+=2*z/(2*j+1);z*=t*t
    return down(s),up(s+2*z/(81*(1-t*t)))
@lru_cache(None)
def logint(x):
    x=rat(x)
    if x<=0:raise ValueError('log argument must be positive')
    k=x.numerator.bit_length()-x.denominator.bit_length()
    y=x/(F(2**k) if k>=0 else F(1,2**(-k)))
    if y<1:k-=1;y*=2
    if y>=2:k+=1;y/=2
    l,u=nearlog(y);a,b=nearlog(F(2))
    return (l+k*a,u+k*b) if k>=0 else (l+k*b,u+k*a)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def solve(A,b):
    n=len(b);M=[list(map(F,a))+[F(t)] for a,t in zip(A,b)]
    for j in range(n):
        k=next((k for k in range(j,n) if M[k][j]),None)
        if k is None:return None
        M[j],M[k]=M[k],M[j];d=M[j][j];M[j]=[v/d for v in M[j]]
        for k in range(n):
            if k!=j and M[k][j]:
                d=M[k][j];M[k]=[x-d*y for x,y in zip(M[k],M[j])]
    return tuple(row[-1] for row in M)
def constraints(V,box):
    H=len(V[0]);A=[];b=[]
    for h in range(H):A.append(tuple(-int(j==h) for j in range(H)));b.append(F(0))
    for v,(l,u) in zip(V,box):
        A+=[tuple(v),tuple(-x for x in v)];b+=[u,-l]
    return A,b
@lru_cache(None)
def vertices(V,box):
    H=len(V[0]);A,b=constraints(V,box);out=set()
    if comb(len(b),H-1)>200000:raise RuntimeError('Exact vertex-enumeration resource limit')
    for ids in combinations(range(len(b)),H-1):
        w=solve([(1,)*H]+[A[i] for i in ids],[F(1)]+[b[i] for i in ids])
        if w is not None and all(dot(a,w)<=z for a,z in zip(A,b)):out.add(w)
    return tuple(sorted(out))
def one_cut_bound(values,excess):
    """Max linear objective on a vertex-mixture simplex with one halfspace.
    Exact exhaustive vertices of this auxiliary simplex: one or two weights.
    """
    pos=[i for i,e in enumerate(excess) if e>=0];neg=[i for i,e in enumerate(excess) if e<0]
    if not pos:return None
    candidates=[values[i] for i in pos]
    for i in pos:
        for j in neg:
            e,f=excess[i],excess[j];candidates.append(((-f)*values[i]+e*values[j])/(e-f))
    return max(candidates)
def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x

def one_cut_dual(values,excess):
    primal=one_cut_bound(values,excess)
    if primal is None:return None,None
    candidates={F(0)}
    for i in range(len(values)):
        for j in range(i):
            if excess[i]!=excess[j]:
                lam=(values[j]-values[i])/(excess[i]-excess[j])
                if lam>=0:candidates.add(lam)
    pairs=[(max(v+lam*e for v,e in zip(values,excess)),lam) for lam in candidates]
    dual,lam=min(pairs)
    if dual!=primal:raise ArithmeticError('Exact auxiliary LP primal-dual mismatch')
    return dual,lam
