"""Stage 19: one declared coarsened observation model for fit and inference.

Full real-space cells, singleton boundaries and overflow intervals are retained.
Only exact (all-record A, V, all-target c) column equivalence is used.
The training predictor is defined on observation equivalence classes first;
target refinement lifts it without resetting its initial aggregate masses.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import sys, json
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(200000)
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0,str(Path(__file__).resolve().parent/'vendor14'))
from coarsened_model import build_model,quotient,record,aggregate_masses
from split_inference import aggregate,conditional_probabilities,conditional_likelihood
from exact import dot,encode,logint

def em_step(groups,counts,m):
    A,V,nn=aggregate(groups,counts);m=tuple(map(F,m));H=len(m)
    aa=[dot(a,m) for a in A];qq=[dot(v,m) for v in V]
    if min(aa)<=0 or min(qq)<=0: raise ValueError('Inadmissible EM input')
    K=sum((F(n)/q for n,q in zip(nn,qq)),F(0))
    nu=tuple(m[h]*sum((n*(F(a[h])/a0+F(1-v[h])/q0)
         for a,v,n,a0,q0 in zip(A,V,nn,aa,qq)),F(0)) for h in range(H))
    out=tuple(x/K for x in nu)
    assert sum(out)==1 and min(out)>=0
    return out,nu,K

def observation_initial(groups):
    return (F(1,len(groups[0]['V'])),)*len(groups[0]['V'])

def compile_input(spec,target_index=0):
    model=build_model(spec['edges'],spec['limits'],spec['targets'])
    obs,oc=quotient(model)
    groups,classes=quotient(model,model['targets_coefficients'])
    H=len(classes);original_index={h:k for k,C in enumerate(oc) for h in C}
    parent=[original_index[C[0]] for C in classes]
    for C,k in zip(classes,parent): assert all(original_index[h]==k for h in C)
    multiplicity=[parent.count(k) for k in range(len(oc))]
    base=observation_initial(obs);predictor=base
    history=[]
    for _ in range(spec.get('training_steps',2)):
        predictor,nu,K=em_step(obs,spec['train'],predictor)
        history.append(predictor)
    lifted=tuple(predictor[k]/multiplicity[k] for k in parent)
    initial_lift=tuple(base[k]/multiplicity[k] for k in parent)
    laws=conditional_probabilities(groups,lifted)
    assert laws==conditional_probabilities(obs,predictor)
    eps=tuple(F(x) for x in spec['epsilon']);alpha=F(spec['alpha'])
    if len(eps)!=len(groups) or any(not 0<=x<1 for x in eps): raise ValueError('Explicit finite tail bounds required')
    seed=lifted
    seed_policy=spec.get('initial_seed_policy','prediction')
    if seed_policy=='remove_completely_invisible_mass':
        invisible=[h for h in range(H) if all(g['V'][h]==0 for g in groups)]
        theta=sum((lifted[h] for h in invisible),F(0))
        if theta>=1:raise ValueError('No visible predictive mass')
        seed=tuple(F(0) if h in invisible else lifted[h]/(1-theta) for h in range(H))
        assert conditional_probabilities(groups,seed)==laws
    elif seed_policy!='prediction':raise ValueError('Unknown explicit seed policy')
    if any(dot(g['V'],seed)<1-e for g,e in zip(groups,eps)):
        raise ValueError('Initial witness is outside the declared tail class; no implicit repair')
    tau=alpha*conditional_likelihood(laws,spec['validation'])
    if tau<=0:raise ValueError('Zero prediction requires the noninformative route, not this positive-threshold pilot')
    A=[];ns=[];js=[]
    for j,(g,counts) in enumerate(zip(groups,spec['validation'])):
        if sum(counts)<=0:raise ValueError('Every validation stratum must be nonempty')
        for r,n in zip(g['records'],counts):
            if n:A.append(r['A']);ns.append(n);js.append(j)
    c=tuple(model['targets_coefficients'][target_index][C[0]] for C in classes)
    data={'name':spec['name']+'_target'+str(target_index),'design':spec,'target_index':target_index,
      'points':tuple(model['representatives'][C[0]] for C in classes),
      'limits':model['limits'],'target':model['targets'][target_index],
      'all_records':[[r['label'] for r in g['records']] for g in groups],
      'all_record_operators':[[r['A'] for r in g['records']] for g in groups],
      'train':spec['train'],'validation':spec['validation'],
      'A':A,'V':[g['V'] for g in groups],'counts':ns,'row_stratum':js,'c':c,
      'epsilon':eps,'alpha':alpha,'prediction':laws,'prediction_mass':lifted,'initial_mass':seed,'tau':tau,
      'initial_seed_policy':seed_policy,
      'support_scope':'full_R2_deterministic_interval_records','predictor':'two_EM_on_observation_quotient_then_lift',
      'cell_geometry':[[I.encoded() for I in cell] for cell in model['cells']],
      'classes':classes,'observation_classes':oc,'class_parent':parent,
      'class_multiplicity':multiplicity,'initial_lift':initial_lift,
      'observation_predictor':predictor,'observation_initial':base}
    return data,model,obs,groups

def numeric_em(groups,counts,initial,max_iter=20000,tol=1e-10):
    """Original unconstrained map. No tail enforcement or post-normalization."""
    aa,vv,nn=aggregate(groups,counts);A=np.asarray(aa,float);V=np.asarray(vv,float);n=np.asarray(nn,float)
    m=np.asarray(initial,float);history=[]
    for k in range(max_iter+1):
        a=A@m;q=V@m
        if min(a)<=0 or min(q)<=0:raise ArithmeticError('Positive observation probabilities lost')
        score=(n/a)@A-(n/q)@V
        ell=float(n@(np.log(a)-np.log(q)))
        residual=max(0.,float(score.max()))/n.sum()
        history.append({'iteration':k,'ell':ell,'positive_score':residual,
          'complementarity':float(np.max(np.abs(m*score)))/n.sum(),'mass_sum':float(m.sum()),
          'mass':m.tolist(),'selection':[float(np.dot(g['V'],m)) for g in groups]})
        if residual<=tol:break
        K=float((n/q).sum());m=m*((n/a)@A+(n/q)@(1-V))/K
    return m,history

def constrained_gem_step(groups,counts,m,eps):
    """Explicit NEW constrained generalized M-step; not the old closed update.
    Search proposes a rational mass; exact constraints, surrogate ascent and
    original-likelihood ascent must all pass before it is accepted.
    """
    eps=tuple(map(F,eps));V=tuple(g['V'] for g in groups);m=tuple(map(F,m));H=len(m)
    if any(dot(v,m)<1-e for v,e in zip(V,eps)):raise ValueError('Initial mass violates tail model')
    closed,nu,K=em_step(groups,counts,m);w=np.asarray([float(x/K) for x in nu]);VV=np.asarray(V,float)
    lower=np.asarray([float(1-e) for e in eps]);x=np.asarray(m,float)
    con=[{'type':'eq','fun':lambda z:z.sum()-1,'jac':lambda z:np.ones(H)},
         {'type':'ineq','fun':lambda z:VV@z-lower,'jac':lambda z:VV}]
    def fun(z):
        if np.any(z<=0):return 1e90
        return -float(w@np.log(z))
    res=minimize(fun,x,jac=lambda z:-w/z,method='SLSQP',bounds=[(1e-14,1)]*H,
                 constraints=con,options={'ftol':1e-13,'maxiter':500})
    z=tuple(F(round(float(a)*10**12),10**12) for a in res.x)
    if min(z)<=0: return m,{'accepted':False,'reason':'nonpositive proposal'}
    z=tuple(a/sum(z) for a in z)
    old=conditional_likelihood(conditional_probabilities(groups,m),counts)
    # A fixed, declared strict interior point is used only to repair a NEW
    # proposed rational mass before exact feasibility/ascent checks.
    common=[h for h in range(H) if all(v[h] for v in V)]
    anchor=list(m)
    if common:
        positive=[e for v,e in zip(V,eps) if not all(v)]
        weight=min(positive+[F(1,2)])/2
        anchor=[weight/F(H)]*H;anchor[common[0]]+=1-weight
    anchor=tuple(anchor)
    for mix in (F(0),F(1,10**10),F(1,10**8),F(1,10**6),F(1,10000),F(1,100),F(1,10)):
        u=tuple((1-mix)*a+mix*b for a,b in zip(z,anchor))
        if any(dot(v,u)<1-e for v,e in zip(V,eps)):continue
        logs=[logint(u0/m0) for u0,m0 in zip(u,m)]
        lo=sum((nh*l[0] for nh,l in zip(nu,logs)),F(0))
        hi=sum((nh*l[1] for nh,l in zip(nu,logs)),F(0))
        new=conditional_likelihood(conditional_probabilities(groups,u),counts)
        if lo>=0 and new>=old:
            return u,{'accepted':True,'surrogate_gain_interval':(lo,hi),'old_likelihood':old,
                      'new_likelihood':new,'new_mass':u,'old_mass':m,'nu':nu,'K':K,
                      'q':tuple(dot(v,u) for v in V),'mix':mix,'proposal_success':bool(res.success)}
    return m,{'accepted':False,'reason':'no certified ascent proposal'}

def nested_profile(counts):
    """Independent analytic factorization for THIS two-limit coarsened design.
    Counts order: shallow C,B,A; deeper C,B_low,B_mid,A_low,A_mid.
    No truth masses or true support enter it.
    """
    s,d=counts
    if len(s)!=3 or len(d)!=5 or min(s+d)<=0:raise ValueError('Positive-count nested diagnostic only')
    N1=sum(s);N2=sum(d);k=d[2];M=N1+N2-k
    r=(F(s[0]+d[0]+d[4],M),F(s[1]+d[1],M),F(s[2]+d[3],M))
    a=F(d[0],d[0]+d[4]);mode=F(k,N2)
    C=(r[0]**(s[0]+d[0]+d[4]))*(r[1]**(s[1]+d[1]))*(r[2]**(s[2]+d[3]))
    C*=a**d[0]*(1-a)**d[4]
    def table(t,theta=F(0)):
        w=(r[0]*a*(1-t),r[1]*(1-t),t,r[2]*(1-t),r[0]*(1-a)*(1-t))
        return ((1-theta)*w[3],(1-theta)*w[1],(1-theta)*w[4],(1-theta)*w[2],(1-theta)*w[0],theta)
    return {'mode':mode,'r':r,'a':a,'C':C,'k':k,'N2':N2,'M':M,'table':table,
            'maximum':C*mode**k*(1-mode)**(N2-k)}

def profile_roots(profile,tau,tol=F(1,10**10)):
    mode=profile['mode'];k=profile['k'];N=profile['N2'];C=profile['C']
    f=lambda t:C*t**k*(1-t)**(N-k)
    assert f(mode)>=tau
    lo=F(0);hi=mode
    while hi-lo>tol:
        mid=(lo+hi)/2
        if f(mid)>=tau:hi=mid
        else:lo=mid
    lower=(lo,hi)
    lo=mode;hi=F(1)
    while hi-lo>tol:
        mid=(lo+hi)/2
        if f(mid)>=tau:lo=mid
        else:hi=mid
    return lower,(lo,hi)
