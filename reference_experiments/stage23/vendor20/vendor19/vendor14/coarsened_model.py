"""Deterministic finite-resolution recording of a continuous latent pair.

All real-space cells, ties, overflow bins, and all non-D outcomes are retained.
Cell masses are coordinates of an arbitrary Borel measure, not a finite true
support or a within-cell uniformity assumption. No noise model is implicit.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import factorial
from typing import Sequence
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'prior'))
from split_inference import dot, conditional_probabilities, em_predictor, conditional_likelihood, count_probability, encode


def rat(x):
    if isinstance(x,float): raise TypeError('Supply exact integers, Fraction, or rational strings')
    return F(x)

@dataclass(frozen=True)
class Interval:
    lo: F|None
    hi: F|None
    lc: bool=False
    rc: bool=False

    def contains(self,x):
        x=rat(x)
        return (self.lo is None or x>self.lo or (self.lc and x==self.lo)) and (self.hi is None or x<self.hi or (self.rc and x==self.hi))

    def intersect(self,other):
        lows=[v for v in (self.lo,other.lo) if v is not None]
        highs=[v for v in (self.hi,other.hi) if v is not None]
        lo=max(lows) if lows else None; hi=min(highs) if highs else None
        lc=False if lo is None else self.contains(lo) and other.contains(lo)
        rc=False if hi is None else self.contains(hi) and other.contains(hi)
        if lo is not None and hi is not None:
            if lo>hi or (lo==hi and not (lc and rc)): return None
        return Interval(lo,hi,lc,rc)

    def representative(self):
        if self.lo is None and self.hi is None:return F(0)
        if self.lo is None:return self.hi-1
        if self.hi is None:return self.lo+1
        return (self.lo+self.hi)/2

    def test_points(self):
        r=self.representative()
        if self.lo==self.hi and self.lo is not None:return [r]
        if self.lo is None and self.hi is None:return [F(-7),F(0),F(7)]
        if self.lo is None:return [self.hi-F(1,7),self.hi-1,self.hi-9]
        if self.hi is None:return [self.lo+F(1,7),self.lo+1,self.lo+9]
        return [(6*self.lo+self.hi)/7,r,(self.lo+6*self.hi)/7]

    def encoded(self):
        return {'lo':None if self.lo is None else str(self.lo),'hi':None if self.hi is None else str(self.hi),'left_closed':self.lc,'right_closed':self.rc}


def make_bins(edges):
    edges=tuple(map(rat,edges))
    if tuple(sorted(set(edges)))!=edges:raise ValueError('Bin edges must be strictly increasing')
    return tuple(Interval(a,b,False,b is not None) for a,b in zip((None,)+edges,edges+(None,)))


def partition(cuts):
    cuts=sorted(set(map(rat,cuts)))
    if not cuts:return (Interval(None,None),)
    out=[Interval(None,cuts[0])]
    for j,x in enumerate(cuts):
        out.append(Interval(x,x,True,True))
        out.append(Interval(x,cuts[j+1] if j+1<len(cuts) else None))
    return tuple(out)


def record(z,c,bins):
    """Detection uses the exact threshold BEFORE reporting a bin index."""
    delta=tuple(int(z[b]<=c[b]) for b in (0,1))
    if not any(delta):return None
    labels=[]
    for b in (0,1):
        labels.append(next(k for k,I in enumerate(bins[b]) if I.contains(z[b])) if delta[b] else -1)
    return delta+tuple(labels)


def record_rectangle(label,c,bins):
    out=[]
    for b in (0,1):
        if label[b]:
            I=bins[b][label[2+b]].intersect(Interval(None,c[b],False,True))
        else:I=Interval(c[b],None)
        if I is None:raise ValueError('Impossible record')
        out.append(I)
    return tuple(out)


def build_model(edges,limits,targets=()):
    edges=tuple(tuple(map(rat,e)) for e in edges)
    cc=tuple(tuple(map(rat,c)) for c in limits)
    tt=tuple(tuple(map(rat,t)) for t in targets)
    if len(edges)!=2 or not cc or any(len(c)!=2 for c in cc+tt):raise ValueError('Bivariate finite design required')
    bins=tuple(make_bins(e) for e in edges)
    axes=tuple(partition(list(edges[b])+[c[b] for c in cc]+[t[b] for t in tt]) for b in (0,1))
    cells=tuple(product(*axes));reps=tuple(tuple(I.representative() for I in cell) for cell in cells)
    groups=[]
    for c in cc:
        labels=tuple(record(z,c,bins) for z in reps)
        records=[]
        for r in sorted({x for x in labels if x is not None}):
            records.append({'label':r,'delta':r[:2],'A':tuple(int(x==r) for x in labels)})
        V=tuple(int(x is not None) for x in labels)
        if not records:raise AssertionError('Every finite limit admits a nonempty real-space region')
        assert all(sum(r['A'][h] for r in records)==V[h] for h in range(len(cells)))
        groups.append({'limit':c,'V':V,'records':records})
    coeff=tuple(tuple(int(z[0]>t[0] and z[1]>t[1]) for z in reps) for t in tt)
    return {'edges':edges,'bins':bins,'limits':cc,'targets':tt,'axes':axes,'cells':cells,'representatives':reps,'groups':groups,'targets_coefficients':coeff}


def quotient(model,extra_coeffs=()):
    """Only EXACT column equivalence. No dominance or zero-truth pruning."""
    groups=model['groups'];H=len(groups[0]['V']);rows=[]
    for g in groups:rows.extend([g['V']]+[r['A'] for r in g['records']])
    rows.extend(extra_coeffs)
    classes={}
    for h in range(H):classes.setdefault(tuple(row[h] for row in rows),[]).append(h)
    classes=list(classes.values());reps=[x[0] for x in classes]
    gg=[]
    for g in groups:
        gg.append({'limit':g['limit'],'V':tuple(g['V'][h] for h in reps),
                   'records':[{**r,'A':tuple(r['A'][h] for h in reps)} for r in g['records']]})
    return gg,classes


def aggregate_masses(m,classes):return tuple(sum((m[h] for h in cc),F(0)) for cc in classes)

class FGMLaw:
    """Absolutely continuous FGM law on (0,4)^2, uniform margins.

    f(x,y)=1/16 [1+lambda(1-x/2)(1-y/2)] there, zero elsewhere.
    It is only used by the evaluator; the estimator never receives it.
    """
    def __init__(self,lam):
        self.lam=rat(lam)
        if abs(self.lam)>1:raise ValueError('FGM lambda must be in [-1,1]')
    def cdf(self,x,y):
        def normalized(v):
            if v is None:return F(1)
            return min(F(1),max(F(0),v/4))
        u,v=normalized(x),normalized(y)
        return u*v*(1+self.lam*(1-u)*(1-v))
    def rectangle(self,rect):
        x,y=rect
        if x.lo==x.hi and x.lo is not None or y.lo==y.hi and y.lo is not None:return F(0)
        def C(a,b):
            if a=='minus' or b=='minus':return F(0)
            return self.cdf(a,b)
        a=x.lo if x.lo is not None else 'minus';b=y.lo if y.lo is not None else 'minus'
        return C(x.hi,y.hi)-C(a,y.hi)-C(x.hi,b)+C(a,b)
    def survival(self,t):
        u=min(F(1),max(F(0),rat(t[0])/4));v=min(F(1),max(F(0),rat(t[1])/4))
        return (1-u)*(1-v)*(1+self.lam*u*v)
    def direct_law(self,c,bins):
        # Independently enumerate products of the coordinate reporting regions.
        coordinate=[]
        for b in (0,1):
            possibilities=[(0,-1,Interval(c[b],None))]
            for j,I in enumerate(bins[b]):
                B=I.intersect(Interval(None,c[b],False,True))
                if B is not None:possibilities.append((1,j,B))
            coordinate.append(possibilities)
        result={};missing=F(0)
        for p,q in product(*coordinate):
            value=self.rectangle((p[2],q[2]))
            if p[0]==q[0]==0:missing=value
            else:result[(p[0],q[0],p[1],q[1])]=value
        return result,1-missing

class UniformRectangle:
    def __init__(self,lo,hi):
        self.lo=tuple(map(rat,lo));self.hi=tuple(map(rat,hi))
        if any(a>=b for a,b in zip(self.lo,self.hi)):raise ValueError('Nonempty rectangle')
    def rectangle(self,rect):
        ans=F(1)
        for i,I in enumerate(rect):
            a=max(self.lo[i],I.lo) if I.lo is not None else self.lo[i]
            b=min(self.hi[i],I.hi) if I.hi is not None else self.hi[i]
            ans*=max(F(0),b-a)/(self.hi[i]-self.lo[i])
        return ans
    def survival(self,t):return self.rectangle(tuple(Interval(rat(x),None) for x in t))


def measure_masses(model,law):return tuple(law.rectangle(Q) for Q in model['cells'])


def compositions(n,k):
    if k==1:yield (n,);return
    for a in range(n+1):
        for rest in compositions(n-a,k-1):yield (a,)+rest


def positive_count_tables(laws,per):
    """Enumeration uses truth only to skip events of probability exactly zero.

    Count vectors retain every possible declared outcome, including all zeros.
    The training predictor receives neither the truth nor the positive indices.
    """
    possibilities=[]
    for law in laws:
        pos=[j for j,p in enumerate(law) if p>0];rows=[]
        for counts in compositions(per,len(pos)):
            out=[0]*len(law)
            for j,k in zip(pos,counts):out[j]=k
            rows.append(tuple(out))
        possibilities.append(rows)
    return list(product(*possibilities))


def exact_coverage(groups,true_mass,per,alpha=F(1,20),steps=2):
    law0=conditional_probabilities(groups,true_mass)
    hist=positive_count_tables(law0,per)
    tests=[(va,count_probability(law0,va),conditional_likelihood(law0,va)) for va in hist]
    assert sum(v[1] for v in tests)==1
    rejection=F(0);zero=F(0);worst=F(0);mean=F(0);maxE=F(0);pairs=0;details=[]
    for tr,pr,_ in tests:
        pred=em_predictor(groups,tr,steps)
        law=conditional_probabilities(groups,pred)
        fail=F(0);ezero=F(0);E=F(0)
        # Powers are cached to make exhaustive exact enumeration practical.
        powers=[[[p**k for k in range(per+1)] for p in row] for row in law]
        for va,pv,L0 in tests:
            D=F(1)
            for gg,nn in enumerate(va):
                for j,k in enumerate(nn):
                    if k:D*=powers[gg][j][k]
            if L0<alpha*D:fail+=pv
            if D==0:ezero+=pv
            E+=pv*(D/L0);pairs+=1
        assert fail<=alpha and E<=1
        rejection+=pr*fail;zero+=pr*ezero;mean+=pr*E
        worst=max(worst,fail);maxE=max(maxE,E)
        details.append({'train_counts':tr,'predictor':pred,'conditional_rejection':fail,'conditional_e_expectation':E})
    return {'per_limit_per_half':per,'n_train':per*len(groups),'n_validation':per*len(groups),
            'alpha':alpha,'coverage':1-rejection,'rejection':rejection,'worst_conditional_rejection':worst,
            'zero_prediction_probability':zero,'mean_e_expectation':mean,'maximum_conditional_e':maxE,
            'positive_count_pairs':pairs,'training_tables':len(hist),'training_details':details}
