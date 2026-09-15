"""Conditional split-likelihood inference on a declared finite observation model.

The statistical confidence set is NOT the Stage 12 kappa-near-optimal set.
Training-only predictions and the *full* conditional record law are essential.
This module does not calibrate an outcome-dependent continuous point-mass
likelihood, select latent support from the validation data, or modify the EM map.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product
from math import comb
from typing import Sequence


def dot(a, b):
    if len(a) != len(b):
        raise ValueError('Dimension mismatch')
    return sum((F(x)*F(y) for x,y in zip(a,b)), F(0))


def observation_model(support, limits):
    """All possible non-D records, including unobserved outcomes, by limit."""
    z = tuple(tuple(F(x) for x in p) for p in support)
    cc = tuple(tuple(F(x) for x in p) for p in limits)
    if not z or len(set(z))!=len(z) or any(len(p)!=2 for p in z+cc):
        raise ValueError('Distinct finite bivariate support and limits required')
    groups=[]
    for c in cc:
        records={}; vis=[]
        for j,p in enumerate(z):
            delta=tuple(int(p[b]<=c[b]) for b in range(2))
            t=tuple(min(p[b],c[b]) for b in range(2))
            visible=int(any(delta));vis.append(visible)
            if visible:
                records.setdefault((t,delta),[0]*len(z))[j]=1
        if not records: raise ValueError('A design stratum has no possible record')
        ordered=sorted(records.items())
        group={'limit':c,'V':tuple(vis),
               'records':[{'t':t,'delta':d,'A':tuple(row)} for (t,d),row in ordered]}
        assert all(sum(r['A'][j] for r in group['records'])==vis[j] for j in range(len(z)))
        groups.append(group)
    return z,groups


def conditional_probabilities(groups, mass):
    m=tuple(map(F,mass))
    if min(m)<0 or sum(m)!=1: raise ValueError('Not a probability mass')
    laws=[]
    for g in groups:
        q=dot(g['V'],m)
        if q<=0: raise ValueError('Zero selection probability for a design stratum')
        pp=tuple(dot(r['A'],m)/q for r in g['records'])
        if sum(pp)!=1 or min(pp)<0: raise AssertionError('Conditional law not normalized')
        laws.append(pp)
    return tuple(laws)


def aggregate(groups, count_rows):
    if len(groups)!=len(count_rows): raise ValueError('Wrong number of strata')
    A=[];V=[];counts=[]
    for g,row in zip(groups,count_rows):
        if len(row)!=len(g['records']) or any(type(n)!=int or n<0 for n in row):
            raise ValueError('Invalid counts')
        for r,n in zip(g['records'],row):
            if n:
                A.append(r['A']);V.append(g['V']);counts.append(n)
    if not counts: raise ValueError('Nonempty training or validation sample required')
    return tuple(A),tuple(V),tuple(counts)


def em_predictor(groups, train_counts, steps=2):
    """Exactly two (by default) unchanged nonnegative EM updates from uniform.

    No held-out outcomes enter this function. Zero predicted outcome probabilities
    are allowed: their realized held-out likelihood is zero and yields no test.
    """
    if type(steps)!=int or steps<0: raise ValueError('Invalid iteration count')
    A,V,counts=aggregate(groups,train_counts); H=len(A[0]);m=(F(1,H),)*H
    for _ in range(steps):
        aa=[dot(a,m) for a in A];qq=[dot(v,m) for v in V]
        if min(aa)<=0 or min(qq)<=0: raise AssertionError('Training record inadmissible')
        K=sum((F(n)/q for n,q in zip(counts,qq)),F(0))
        m=tuple(m[h]*sum((F(n)*(F(a[h])/ap+F(1-v[h])/qp)
                for a,v,n,ap,qp in zip(A,V,counts,aa,qq)),F(0))/K for h in range(H))
        if min(m)<0 or sum(m)!=1: raise AssertionError('EM normalization failed')
    conditional_probabilities(groups,m)
    return m


def conditional_likelihood(laws, counts):
    result=F(1)
    for pp,nn in zip(laws,counts):
        for p,n in zip(pp,nn):
            if n: result*=p**n
    return result


def count_probability(laws, counts):
    """Independent stratum multinomials; currently audit uses binary strata."""
    from math import factorial
    c=1
    for nn in counts:
        x=factorial(sum(nn))
        for k in nn:x//=factorial(k)
        c*=x
    return c*conditional_likelihood(laws,counts)


def binary_histories(per_stratum, n_strata):
    if type(per_stratum)!=int or per_stratum<1: raise ValueError('Positive size')
    return [tuple((k,per_stratum-k) for k in kk)
            for kk in product(range(per_stratum+1),repeat=n_strata)]


def exact_split_audit(groups, true_mass, per_train, per_test, alpha=F(1,20), steps=2):
    """Enumerate EVERY stratum-count table; no Monte Carlo or asymptotic cutoff."""
    if any(len(g['records'])!=2 for g in groups): raise ValueError('Binary audit only')
    alpha=F(alpha)
    if not 0<alpha<1: raise ValueError('alpha in (0,1)')
    law0=conditional_probabilities(groups,true_mass)
    trs=binary_histories(per_train,len(groups));vas=binary_histories(per_test,len(groups))
    trains=[(t,count_probability(law0,t),em_predictor(groups,t,steps)) for t in trs]
    tests=[(v,count_probability(law0,v),conditional_likelihood(law0,v)) for v in vas]
    assert sum(p for _,p,_ in trains)==sum(p for _,p,_ in tests)==1
    failure=F(0);zero=F(0);meanE=F(0);maxcond=F(0);mincond=F(1)
    maxE=F(0);minE=F(1);tested=0;details=[]
    for tr,pr,m in trains:
        if not pr:continue
        pred=conditional_probabilities(groups,m)
        condfail=F(0);condzero=F(0);condE=F(0)
        for va,pv,L0 in tests:
            if not pv:continue
            assert L0>0
            D=conditional_likelihood(pred,va)
            E=D/L0
            condE+=pv*E
            if L0<alpha*D:condfail+=pv
            if D==0:condzero+=pv
            tested+=1
        assert 0<=condE<=1 and 0<=condfail<=alpha
        failure+=pr*condfail;zero+=pr*condzero;meanE+=pr*condE
        maxcond=max(maxcond,condfail);mincond=min(mincond,condfail)
        maxE=max(maxE,condE);minE=min(minE,condE)
        details.append({'training_counts':tr,'predictive_mass':m,
                        'conditional_rejection_probability':condfail,
                        'conditional_e_expectation':condE})
    assert failure<=alpha and 0<=meanE<=1
    return {'true_mass':tuple(map(F,true_mass)), 'train_per_limit':per_train,
            'test_per_limit':per_test, 'n_train':per_train*len(groups),
            'n_validation':per_test*len(groups),'alpha':alpha,'em_steps':steps,
            'coverage':1-failure,'rejection_probability':failure,
            'worst_conditional_rejection':maxcond,
            'mean_e_expectation':meanE,'min_conditional_e':minE,'max_conditional_e':maxE,
            'zero_predictive_likelihood_probability':zero,
            'positive_probability_count_pairs':tested,
            'all_conditional_markov_checks_pass':True,'training_details':details}


def fixed_kappa_binomial_coverage(n, p=F(1,2), kappa=F(99,100)):
    p=F(p);kappa=F(kappa)
    if not 0<p<1 or not 0<kappa<1:raise ValueError('Interior parameters')
    coverage=F(0); accepted=[]
    for k in range(n+1):
        L0=p**k*(1-p)**(n-k)
        q=F(k,n)
        Lmax=q**k*(1-q)**(n-k)
        if L0>=kappa*Lmax:
            coverage+=comb(n,k)*L0;accepted.append(k)
    return {'n':n,'true_p':p,'kappa':kappa,'coverage':coverage,'accepted_counts':accepted}


def encode(obj):
    if isinstance(obj,F): return str(obj)
    if isinstance(obj,dict): return {str(k):encode(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)): return [encode(v) for v in obj]
    return obj
