"""Certified mean log-likelihood global upper bounds for saved Stage24 cases.
The only numerical roles are proposing tangent points and LP multipliers.
Every upper bound is completed and evaluated with exact rational arithmetic.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import sys,json,time,heapq,hashlib
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize,linprog
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'reference_stage24'))
from core import compile_model,feasible
from exact import logint,dot,encode


def prepare(record):
    d=compile_model(record['public_input'])
    assert encode({k:d[k] for k in record['compiled']})==record['compiled']
    candidate=tuple(map(F,record['candidate']['mass']))
    assert feasible(d,candidate)
    # Deepest visible set is the union of all visible sets.
    eta=d['V'][2]
    assert all(max(v[h] for v in d['V'])==eta[h] for h in range(len(eta)))
    ids=[h for h,e in enumerate(eta) if e]
    assert all(not A[h] for h,e in enumerate(eta) if not e for A in d['A'])
    r=dot(eta,candidate)
    ref=tuple(candidate[h]/r for h in ids)
    A=tuple(tuple(a[h] for h in ids) for a in d['A'])
    V=tuple(tuple(v[h] for h in ids) for v in d['V'])
    VR=tuple(tuple(v[h] for h in ids) for v in d['VR'])
    assert all(x==1 for x in V[2])
    n=d['counts'];N=d['N'];omega=tuple(F(x,N) for x in n)
    layers=tuple(next(j for j,v in enumerate(V) if v==vv) for vv in VR)
    lam=tuple(sum((w for w,jj in zip(omega,layers) if jj==j),F(0)) for j in range(3))
    assert all(dot(a,candidate)/dot(v,candidate)==dot(aa,ref)/dot(vv,ref) for a,v,aa,vv in zip(d['A'],d['VR'],A,VR))
    assert min(dot(v,ref) for v in V)>=F(1,2)
    return {'A':A,'V':V,'VR':VR,'n':n,'N':N,'omega':omega,'layers':layers,'lam':lam,'ref':ref,'candidate':candidate,'visible_ids':ids}


def loglik_interval(d,m):
    low=high=F(0)
    for a,v,w in zip(d['A'],d['VR'],d['omega']):
        l,h=logint(dot(a,m)/dot(v,m));low+=w*l;high+=w*h
    return low,high


def tangent_coeff(d,ref,box):
    # Rational affine upper bound valid on the box. Logarithm upper endpoints
    # appear with nonnegative interpolation weights on the actual q-domain.
    c=[F(0)]*len(ref);base=F(0)
    for a,w in zip(d['A'],d['omega']):
        t=dot(a,ref);assert t>0
        base+=w*(logint(t)[1]-1)
        for h,ah in enumerate(a):c[h]+=w*ah/t
    for j,(l,u) in enumerate(box):
        L=-logint(l)[0];U=-logint(u)[0]
        if l==u:base+=d['lam'][j]*L
        else:
            slope=(U-L)/(u-l);intercept=L-slope*l
            base+=d['lam'][j]*intercept
            for h,vh in enumerate(d['V'][j]):c[h]+=d['lam'][j]*slope*vh
    return tuple(x+base for x in c)


def polytope(d,box):
    G=[];b=[]
    for j,(l,u) in enumerate(box):
        G.extend([d['V'][j],tuple(-v for v in d['V'][j])]);b.extend([u,-l])
    return G,b


def certify_tangent(d,ref,box):
    t=tangent_coeff(d,ref,box);G,b=polytope(d,box)
    sol=linprog(-np.array(t,float),A_ub=np.array(G,float),b_ub=np.array(b,float),A_eq=np.ones((1,len(ref))),b_eq=[1],bounds=(0,None),method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
    if sol.success:y=tuple(F(format(max(0.,-float(v)),'.16g')) for v in sol.ineqlin.marginals)
    else:y=(F(0),)*len(b)
    z=max(t[h]-sum((y[j]*G[j][h] for j in range(len(b))),F(0)) for h in range(len(ref)))
    upper=z+dot(y,b)
    return {'reference':ref,'multipliers':y,'simplex_multiplier':z,'upper':upper,'proposal_LP_success':bool(sol.success)}


def propose(d,box,start):
    A=np.array(d['A'],float);V=np.array(d['V'][:2],float);w=np.array(d['omega'],float)
    l=np.array([float(x[0]) for x in box]);u=np.array([float(x[1]) for x in box]);lam=np.array(d['lam'][:2],float)
    slopes=np.array([0. if a==b else -np.log1p(float((b-a)/a))/float(b-a) for a,b in box])
    intercept=-np.log(l)-slopes*l
    linear=lam*slopes
    def fun(x):
        aa=A@x
        if np.min(aa)<=0:return 1e10
        return -float(w@np.log(aa)+linear@(V@x)+lam@intercept)
    def jac(x):return -(w/np.maximum(A@x,1e-100))@A-linear@V
    constraints=[{'type':'eq','fun':lambda x:x.sum()-1,'jac':lambda x:np.ones(len(x))},
                 {'type':'ineq','fun':lambda x:np.r_[V@x-l,u-V@x,A@x-1e-12],
                  'jac':lambda x:np.vstack([V,-V,A])}]
    result=minimize(fun,np.array(start,float),jac=jac,bounds=[(0.,1.)]*len(start),constraints=constraints,method='SLSQP',options={'ftol':3e-13,'maxiter':200})
    x=np.maximum(result.x,0.);x/=sum(x)
    q=tuple(F(format(float(v),'.12g')) for v in x)
    if sum(q)>0:q=tuple(x/sum(q) for x in q)
    if min(dot(a,q) for a in d['A'])<=0:q=d['ref']
    return q,{'success':bool(result.success),'nit':int(result.nit)}


def run_case(path,out,tol=F(1,10**6),max_nodes=4001):
    tic=time.monotonic();record=json.loads(Path(path).read_text());d=prepare(record)
    lower,lh=loglik_interval(d,d['ref']);tree=[];active={};heap=[];proposals=0
    def node(box,start):
        nonlocal proposals
        ref,pinfo=propose(d,box,start);proposals+=1
        cert=certify_tangent(d,ref,box)
        # Previous candidate tangent can be preferable if proposal did not settle.
        alt=certify_tangent(d,d['ref'],box)
        if alt['upper']<cert['upper']:cert=alt;ref=d['ref']
        k=len(tree);tree.append({'id':k,'box':box,'proof':cert,'proposal':pinfo})
        active[k]=True;heapq.heappush(heap,(-cert['upper'],k))
        return k,ref
    node(((F(1,2),F(1)),(F(1,2),F(1))),d['ref'])
    while heap and len(tree)+2<=max_nodes:
        ub,k=heap[0];ub=-ub
        if ub-lower<=tol:break
        heapq.heappop(heap);active.pop(k);p=tree[k];box=p['box'];start=p['proof']['reference']
        j=max(range(2),key=lambda j:d['lam'][j]*(box[j][1]-box[j][0])**2/box[j][0]**2)
        l,u=box[j];mid=(l+u)/2;b1=list(box);b2=list(box);b1[j]=(l,mid);b2[j]=(mid,u)
        a,_=node(tuple(b1),start);b,_=node(tuple(b2),start);p['children']=[a,b];p['split_axis']=j
        if len(tree)%100==1:print(record['case'],'nodes',len(tree),'gap',float(-heap[0][0]-lower),flush=True)
    upper=-heap[0][0]
    cert={'case':record['case'],'source_input_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
       'visible_indices':d['visible_ids'],'normalized_candidate':d['ref'],
       'mean_lower_interval':(lower,lh),'mean_upper':upper,'gap_upper':upper-lower,
       'tolerance':tol,'precision_achieved':upper-lower<=tol,'max_nodes':max_nodes,'tree':tree,
       'seconds':time.monotonic()-tic,'new_candidate_optimization':False}
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/(record['case']+'.json')).write_text(json.dumps(encode(cert),sort_keys=True,indent=2)+'\n')
    print(record['case'],'DONE',len(tree),float(upper-lower),cert['precision_achieved'],'sec',cert['seconds'],flush=True)
    return cert

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--case');p.add_argument('--outdir',default=str(BASE/'outputs/certificates'));p.add_argument('--tol',default='1/1000000');p.add_argument('--max-nodes',type=int,default=4001);args=p.parse_args()
    paths=[BASE/'inputs'/(args.case+'.json')] if args.case else sorted((BASE/'inputs').glob('coarse_*.json'))
    for path in paths:run_case(path,args.outdir,F(args.tol),args.max_nodes)
