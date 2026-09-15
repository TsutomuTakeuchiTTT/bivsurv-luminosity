"""Stage24 numerical proposals and exact rational constrained-M-step audit.
Two crossing exclusion sets have the DECLARED bounds 1/2,1/2;
the third exclusion set is contained in both. Not a general J-constraint solver.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,json
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import numpy as np
from scipy.optimize import minimize,linprog
from scipy.linalg import null_space
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'vendor'))
from coarsened_model import build_model,quotient,aggregate_masses,measure_masses
from split_inference import aggregate,conditional_probabilities
from exact import dot,logint,encode

class ReferenceLaw:
    @staticmethod
    def cdf(x):
        if x is None:return F(1)
        if x<=1:return 1/(3*(2-x))
        if x<=2:return x/3
        return 1-1/(3*(x-1))
    def rectangle(self,cell):
        p=F(1)
        for I in cell:p*=self.cdf(I.hi)-(self.cdf(I.lo) if I.lo is not None else 0)
        return p

def compile_model(inp):
    mod=build_model(inp['edges'],inp['limits'],inp['targets']);groups,C=quotient(mod)
    raw=measure_masses(mod,ReferenceLaw());initial=aggregate_masses(raw,C)
    coeff=tuple(tuple(sum((raw[h]*c[h] for h in cc),F(0))/initial[k] for k,cc in enumerate(C)) for c in mod['targets_coefficients'])
    counts=[[a+b for a,b in zip(x,y)] for x,y in zip(inp['train'],inp['validation'])]
    A,VR,n=aggregate(groups,counts);V=tuple(g['V'] for g in groups);D=tuple(tuple(1-v for v in row) for row in V)
    if tuple(map(F,inp['epsilon']))!=(F(1,2),)*3:raise ValueError('This audit requires the three declared half bounds')
    assert all(D[2][h]<=D[0][h] and D[2][h]<=D[1][h] for h in range(len(initial)))
    sig=tuple(2*D[0][h]+D[1][h] for h in range(len(initial)))
    if set(sig)!={0,1,2,3}:raise ValueError('All four exclusion signatures required for this closed formula')
    return {'A':A,'VR':VR,'counts':n,'V':V,'D':D,'groups':groups,'initial':initial,'sig':sig,'coeff':coeff,'N':sum(n),'model':mod,'classes':C}

def expected(d,m):
    A,VR,n=d['A'],d['VR'],d['counts'];a=[dot(x,m) for x in A];q=[dot(x,m) for x in VR]
    if min(a+q)<=0:raise ValueError('Nonpositive observed probability')
    K=sum((F(k)/y for k,y in zip(n,q)),F(0))
    nu=tuple(m[h]*sum((k*(F(x[h])/xx+F(1-v[h])/yy) for x,v,k,xx,yy in zip(A,VR,n,a,q)),F(0)) for h in range(len(m)))
    assert sum(nu)==K
    return tuple(x/K for x in nu),K

def closed_mstep(w,sig):
    W=tuple(sum((x for x,s in zip(w,sig) if s==j),F(0)) for j in range(4))
    if min(W)<=0:raise ValueError('Positive four-group expected masses required; degenerate case not implemented')
    candidates=[('none',W,(F(0),F(0)),F(1))]
    for j in (0,1):
        bits=[s//2 if j==0 else s%2 for s in range(4)]
        t=sum((W[s] for s in range(4) if bits[s]),F(0))
        P=tuple(W[s]/(2*(t if bits[s] else 1-t)) for s in range(4))
        beta=[F(0),F(0)];beta[j]=4*t-2
        candidates.append(('tail'+str(j+1),P,tuple(beta),2*(1-t)))
    x=(W[0]+W[3])/2;y=F(1,2)-x
    eta=W[0]/x;beta=(W[2]/y-eta,W[1]/y-eta)
    candidates.append(('both',(x,y,y,x),beta,eta))
    for active,P,beta,eta in candidates:
        if min(beta)<0 or P[2]+P[3]>F(1,2) or P[1]+P[3]>F(1,2):continue
        if any(beta[j]*(sum((P[s] for s in range(4) if (s//2 if j==0 else s%2)),F(0))-F(1,2))!=0 for j in (0,1)):continue
        assert all(W[s]/P[s]==eta+beta[0]*(s//2)+beta[1]*(s%2) for s in range(4))
        u=tuple(P[s]*x/W[s] for x,s in zip(w,sig))
        assert sum(u)==1 and min(u)>=0
        return u,{'active':active,'group_weights':W,'group_solution':P,'multipliers':beta,'eta':eta}
    raise ArithmeticError('No exact KKT active set')

def qvals(d,m):return tuple(dot(v,m) for v in d['V'])
def feasible(d,m):
    return min(m)>=0 and sum(m)==1 and min(qvals(d,m))>=F(1,2) and all(dot(a,m)>0 for a in d['A'])

def rounded_candidate(x,d,anchor,digits=22):
    # This is a NEW numerical proposal, not a projection of a previously reported fit.
    # Its feasibility and increase are tested before acceptance.
    with localcontext() as ctx:
        ctx.prec=digits
        u=tuple(F(str(Decimal(z.numerator)/Decimal(z.denominator))) if z else F(0) for z in x)
    u=tuple(z/sum(u) for z in u);mix=F(0)
    for q,q0 in zip(qvals(d,u),qvals(d,anchor)):
        if q<F(1,2):
            if q0<=F(1,2):raise ValueError('Anchor must be strictly feasible')
            mix=max(mix,(F(1,2)-q)/(q0-q))
    # Exact mixing is only to make the rounded trial mass feasible.
    if mix:u=tuple((1-mix)*a+mix*b for a,b in zip(u,anchor))
    assert feasible(d,u)
    return u,mix

def gain_interval(d,old,new):
    lo=hi=F(0)
    for A,V,n in zip(d['A'],d['VR'],d['counts']):
        x=dot(A,new)*dot(V,old)/(dot(A,old)*dot(V,new))
        l,h=logint(x);lo+=n*l;hi+=n*h
    return lo/d['N'],hi/d['N']

def surrogate_interval(w,old,new):
    lo=hi=F(0)
    for a,b,c in zip(w,old,new):
        if a:
            if b<=0 or c<=0:raise ValueError('Undefined surrogate')
            l,h=logint(c/b);lo+=a*l;hi+=a*h
    return lo,hi

def grad_exact(d,m):
    aa=[dot(a,m) for a in d['A']];qq=[dot(v,m) for v in d['VR']]
    return tuple(sum((n*(F(a[h])/x-F(v[h])/y) for a,v,n,x,y in zip(d['A'],d['VR'],d['counts'],aa,qq)),F(0))/d['N'] for h in range(len(m)))

def stationary_certificate(d,m):
    g=grad_exact(d,m);D=d['D'][:2];H=len(m)
    res=linprog(-np.array(g,float),A_ub=np.array(D,float),b_ub=[.5,.5],A_eq=np.ones((1,H)),b_eq=[1],bounds=(0,None),method='highs')
    if not res.success:raise ArithmeticError(res.message)
    beta=tuple(F(format(max(0.,-float(x)),'.16g')) for x in res.ineqlin.marginals)
    z=max(g[h]-sum((beta[j]*D[j][h] for j in (0,1)),F(0)) for h in range(H))
    upper=z+sum(beta,F(0))/2-dot(g,m)
    assert upper>=0 and dot(g,m)==0
    return {'multipliers':beta,'simplex_multiplier':z,'directional_gap_upper':upper,'maximum_raw_score':max(g),'kind':'first_order_stationarity_not_global_likelihood_gap'}

def numeric_problem(d):
    A=np.array(d['A'],float);V=np.array(d['VR'],float);w=np.array(d['counts'],float)/d['N'];DV=np.array(d['D'][:2],float)
    def objective(x):
        a=A@x;q=V@x
        if min(a)<=0 or min(q)<=0:return 1e100
        return -float(w@(np.log(a)-np.log(q)))
    def gradient(x):return -(w/(A@x))@A+(w/(V@x))@V
    def hessian(x):return A.T@((w/(A@x)**2)[:,None]*A)-V.T@((w/(V@x)**2)[:,None]*V)
    return objective,gradient,hessian,DV

def polish(d,z):
    fun,jac,hes,D=numeric_problem(d)
    support=np.flatnonzero(z>1e-9)
    active=np.flatnonzero(abs(D@z-.5)<1e-7)
    E=np.vstack([np.ones(len(support)),D[active][:,support]]);b=np.array([1]+[.5]*len(active))
    v=z[support];v=v+E.T@np.linalg.lstsq(E@E.T,b-E@v,rcond=None)[0]
    if min(v)<0:return z
    B=null_space(E);H=len(z)
    zz=np.zeros(H);zz[support]=v
    for k in range(30):
        gg=B.T@jac(zz)[support]
        if max(abs(gg),default=0)<1e-13:break
        HH=B.T@hes(zz)[np.ix_(support,support)]@B
        step=-B@np.linalg.lstsq(HH,gg,rcond=1e-10)[0]
        if float(jac(zz)[support]@step)>=0:break
        t=1.
        accepted=False
        for _ in range(50):
            vv=zz[support]+t*step;trial=np.zeros(H);trial[support]=vv
            if min(vv)>=0 and max(D@trial)<=.5+1e-12 and fun(trial)<=fun(zz)+1e-14:
                zz=trial;accepted=True;break
            t/=2
        if not accepted:break
    return zz

def independent_candidates(d,gem):
    fun,jac,hes,D=numeric_problem(d);H=len(gem)
    starts=[np.array(gem,float),np.array(d['initial'],float),np.ones(H)/H]
    results=[];best=None
    for start in starts:
        res=minimize(fun,start,jac=jac,method='SLSQP',bounds=[(0,1)]*H,
            constraints=[{'type':'eq','fun':lambda x:x.sum()-1,'jac':lambda x:np.ones(H)},
                         {'type':'ineq','fun':lambda x:.5-D@x,'jac':lambda x:-D}],
            options={'ftol':1e-14,'maxiter':2000})
        x=polish(d,res.x)
        raw=tuple(F(format(max(0.,float(a)),'.17g')) for a in x)
        raw=tuple(a/sum(raw) for a in raw)
        u,mix=rounded_candidate(raw,d,d['initial'])
        cert=stationary_certificate(d,u)
        results.append({'success_proposal_only':bool(res.success),'iterations':int(res.nit),'objective':float(fun(np.array(u,float))),'certificate':cert})
        key=(results[-1]['objective'],float(cert['directional_gap_upper']))
        if best is None or key<best[0]:best=(key,u,cert)
    # Prioritize a certified-small stationary residual among effectively tied values.
    feasible_results=[]
    # best from primary objective, all starts are separately diagnosed.
    return best[1],best[2],results
