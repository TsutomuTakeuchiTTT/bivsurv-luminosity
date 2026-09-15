"""Common-parent multi-limit confidence projection by denominator boxes.

Numerical optimization proposes tangent points and primal witnesses only.
Final outer bounds use exact polytope vertices and rational log enclosures.
Tail restrictions q_j>=1-epsilon_j>0 are explicitly part of this model.
No support columns are deleted for projection.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json,hashlib,math,warnings
import numpy as np
from scipy.optimize import minimize
from exact import rat,logint,dot,vertices,one_cut_bound,one_cut_dual,encode

def digest(x):return hashlib.sha256(json.dumps(encode(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()
def record(z,c):
    de=tuple(int(x<=y) for x,y in zip(z,c))
    if not any(de):return None
    return tuple(z[b] if de[b] else c[b] for b in range(2))+de

def from_observations(points,limits,train,validation,target,eps,alpha=F(1,20),name='case'):
    points=tuple(tuple(map(rat,z)) for z in points);limits=tuple(tuple(map(rat,z)) for z in limits)
    records=[];A=[];V=[];row_stratum=[];counts=[];full=[]
    for j,lim in enumerate(limits):
        lab=sorted(set(record(z,lim) for z in points)-{None})
        if not lab:raise ValueError('Empty selection stratum')
        AA=[tuple(int(record(z,lim)==r) for z in points) for r in lab]
        vv=tuple(int(record(z,lim) is not None) for z in points)
        if any(sum(a[h] for a in AA)!=vv[h] for h in range(len(points))):raise AssertionError('Partition')
        if len(train[j])!=len(lab) or len(validation[j])!=len(lab):raise ValueError(('Counts do not cover all records',j,lab))
        if any(type(x)!=int or x<0 for x in train[j]+validation[j]) or not sum(validation[j]):raise ValueError('Counts')
        V.append(vv);full.append(AA);records.append(lab)
        for a,n in zip(AA,validation[j]):
            if n:A.append(a);counts.append(n);row_stratum.append(j)
    H=len(points);m=[F(1,H)]*H
    # Two unchanged completion-EM training updates. No tail-constrained fitting is claimed.
    for _ in range(2):
        q=[dot(v,m) for v in V];K=sum(F(sum(t),1)/x for t,x in zip(train,q));new=[F(0)]*H
        for j,(AA,nn) in enumerate(zip(full,train)):
            for a,n in zip(AA,nn):
                if not n:continue
                av=dot(a,m)
                for h in range(H):new[h]+=n*m[h]*(F(a[h],1)/av+F(1-V[j][h],1)/q[j])
        m=[x/K for x in new]
    eps=tuple(map(rat,eps));alpha=rat(alpha)
    if len(eps)!=len(V) or any(not 0<=e<1 for e in eps):raise ValueError('Explicit epsilon<1 required')
    if any(dot(v,m)<1-e for v,e in zip(V,eps)):raise ValueError('Training prediction outside declared tail class')
    prediction=[[dot(a,m)/dot(v,m) for a in AA] for AA,v in zip(full,V)]
    tau=alpha
    for pr,nn in zip(prediction,validation):
        for p,n in zip(pr,nn):tau*=p**n
    if tau<=0:raise ValueError('Positive predictive likelihood required for this pilot')
    t=tuple(map(rat,target));c=tuple(int(all(x>y for x,y in zip(z,t))) for z in points)
    data={'name':name,'points':points,'limits':limits,'target':t,'epsilon':eps,'alpha':alpha,
          'all_records':records,'all_record_operators':full,'train':train,'validation':validation,
          'A':A,'V':V,'counts':counts,'row_stratum':row_stratum,'c':c,
          'prediction':prediction,'initial_mass':m,'tau':tau,
          'support_scope':'declared_fixed_finite_support','predictor':'two_training_completion_EM_updates'}
    return data

def prepare(data):
    d=dict(data)
    for key in ('epsilon','initial_mass','c'):d[key]=tuple(map(rat,d[key]))
    d['A']=tuple(tuple(map(int,a)) for a in d['A']);d['V']=tuple(tuple(map(int,v)) for v in d['V'])
    d['tau']=rat(d['tau']);d['counts']=tuple(d['counts']);d['row_stratum']=tuple(d['row_stratum'])
    d['N']=sum(d['counts']);d['stratum_counts']=tuple(sum(d['counts'][i] for i,jj in enumerate(d['row_stratum']) if jj==j) for j in range(len(d['V'])))
    return d

def likelihood(data,m):
    m=tuple(map(rat,m))
    if any(x<0 for x in m) or sum(m)!=1:return None
    q=tuple(dot(v,m) for v in data['V'])
    if any(x<1-e for x,e in zip(q,data['epsilon'])):return None
    L=F(1)
    for a,n,j in zip(data['A'],data['counts'],data['row_stratum']):
        av=dot(a,m)
        if av<=0:return None
        L*=(av/q[j])**n
    return L

def rational_proposal(w,D=10**12):
    # New proposed rational point; never accepted without exact feasibility testing.
    v=[F(round(float(x)*D),D) for x in w]
    if any(x<0 for x in v) or sum(v)<=0:return None
    return tuple(x/sum(v) for x in v)

def tangent_vertex_bounds(data,box,v,verts):
    """Upper affine tangent to concave numerator log plus convex-denominator secants.
    Computes its value at every polytope vertex with rigorous upper rounding.
    """
    av=tuple(dot(a,v) for a in data['A'])
    if any(x<=0 for x in av):raise ValueError('Nonpositive tangent numerator')
    logs=[logint(x) for x in av];N=data['N'];out=[]
    logq={q:logint(q) for interval in box for q in interval}
    for w in verts:
        T=sum((F(n,N)*(lu[1]+dot(a,w)/x-1) for a,n,x,lu in zip(data['A'],data['counts'],av,logs)),F(0))
        for vv,n,(l,u) in zip(data['V'],data['stratum_counts'],box):
            q=dot(vv,w)
            if l==u:sec=-logq[l][0]
            else:sec=-(u-q)/(u-l)*logq[l][0]-(q-l)/(u-l)*logq[u][0]
            T+=F(n,N)*sec
        out.append(T)
    return tuple(out)

def certify_tangent(data,box,v,verts,sign):
    T=tangent_vertex_bounds(data,box,v,verts);tlo=logint(data['tau'])[0]/data['N']
    values=[sign*dot(data['c'],w) for w in verts]
    return one_cut_dual(values,[t-tlo for t in T])

def project(data,tolerance=F(1,10000),max_nodes=2001):
    d=prepare(data);tol=rat(tolerance);H=len(d['c']);J=len(d['V']);N=d['N']
    if tol<=0 or type(max_nodes)!=int or max_nodes<1:raise ValueError('Positive precision and node budget required')
    if len(d['epsilon'])!=J or any(not 0<=e<1 for e in d['epsilon']):raise ValueError('Tail specification')
    A=np.asarray(d['A'],float);V=np.asarray(d['V'],float);w=np.asarray(d['counts'],float)/N;nj=np.asarray(d['stratum_counts'],float)/N
    c=np.asarray(d['c'],float);tl,tu=logint(d['tau']);threshold=float((tl+tu)/(2*N))
    root=tuple((F(1) if all(v) else 1-e,F(1)) for v,e in zip(d['V'],d['epsilon']))
    m0=d['initial_mass'];L0=likelihood(d,m0)
    if L0 is None or L0<d['tau']:raise ValueError('No certified initial witness')
    incumbent={s:(s*dot(d['c'],m0),m0) for s in (-1,1)}
    numerical_calls=0;nodes=[];active={}

    def eval_f(m,box=None):
        a=A@m
        if np.any(a<=0):return -1e50,np.zeros(H)
        val=float(w@np.log(a));grad=(w/a)@A
        q=V@m
        if np.any(q<=0):return -1e50,np.zeros(H)
        if box is None:
            return val-float(nj@np.log(q)),grad-(nj/q)@V
        for j,(ll,uu) in enumerate(box):
            l,u=float(ll),float(uu)
            if l==u:val-=nj[j]*math.log(l)
            else:
                a0=-math.log(l);a1=-math.log(u);s=(a1-a0)/(u-l)
                val+=nj[j]*(a0+s*(q[j]-l));grad+=nj[j]*s*V[j]
        return val,grad

    def constraints_num(box,relaxed):
        lo=np.asarray([float(a) for a,b in box]);hi=np.asarray([float(b) for a,b in box])
        result=[{'type':'eq','fun':lambda m:float(m.sum()-1),'jac':lambda m:np.ones(H)}]
        # Trivial all-visible rows are already equalities by normalization; avoid duplicates.
        ids=[j for j in range(J) if not np.all(V[j]==1)]
        if ids:
            VV=V[ids];l=lo[ids];u=hi[ids]
            result +=[{'type':'ineq','fun':lambda m:VV@m-l,'jac':lambda m:VV},
                      {'type':'ineq','fun':lambda m:u-VV@m,'jac':lambda m:-VV}]
        result +=[{'type':'ineq','fun':lambda m:A@m-1e-13,'jac':lambda m:A}]
        if relaxed is not False:
            b=box if relaxed else None
            result +=[{'type':'ineq','fun':lambda m,b=b:eval_f(m,b)[0]-threshold,
                       'jac':lambda m,b=b:eval_f(m,b)[1]}]
        return result

    def attempt_witness(z):
        p=rational_proposal(z)
        if p is None:return
        for scale in (F(0),F(1,10**10),F(1,10**8),F(1,10**6),F(1,10000),F(1,100),F(1,2)):
            m=tuple((1-scale)*a+scale*b for a,b in zip(p,m0));L=likelihood(d,m)
            if L is not None and L>=d['tau']:
                for s in (-1,1):
                    val=s*dot(d['c'],m)
                    if val>incumbent[s][0]:incumbent[s]=(val,m)
                return

    # Local original optimizers supply only certified feasible lower bounds.
    for s in (-1,1):
        for seed in (np.array(list(map(float,m0))),np.ones(H)/H):
            con=constraints_num(root,False)
            con.append({'type':'ineq','fun':lambda m:eval_f(m,None)[0]-threshold,'jac':lambda m:eval_f(m,None)[1]})
            res=minimize(lambda m:-s*float(c@m),seed,jac=lambda m:-s*c,bounds=[(0,1)]*H,
                         constraints=con,method='SLSQP',options={'ftol':1e-12,'maxiter':600})
            numerical_calls+=1;attempt_witness(res.x)

    def process(box):
        nonlocal numerical_calls
        idx=len(nodes);verts=vertices(d['V'],box)
        node={'id':idx,'box':box,'vertices':verts,'status':'leaf','certificates':{}}
        nodes.append(node)
        if not verts:node['status']='empty_polytope';return idx
        seed=tuple(sum((v[h] for v in verts),F(0))/len(verts) for h in range(H))
        if any(dot(a,seed)==0 for a in d['A']):node['status']='zero_observed_numerator';return idx
        z=np.array(list(map(float,seed)))
        # Feasibility phase solves a concave max proposal if needed.
        if eval_f(z,box)[0]<threshold:
            res=minimize(lambda m:-eval_f(m,box)[0],z,jac=lambda m:-eval_f(m,box)[1],
                         bounds=[(0,1)]*H,constraints=constraints_num(box,False),method='SLSQP',
                         options={'ftol':1e-12,'maxiter':400})
            numerical_calls+=1;z=res.x
        ts=[];upper={}
        for s in (-1,1):
            wi=incumbent[s][1]
            zi=np.array(list(map(float,wi))) if all(l<=dot(v,wi)<=u for v,(l,u) in zip(d['V'],box)) else z
            res=minimize(lambda m:-s*float(c@m),zi,jac=lambda m:-s*c,bounds=[(0,1)]*H,
                         constraints=constraints_num(box,True),method='SLSQP',
                         options={'ftol':1e-12,'maxiter':400})
            numerical_calls+=1;attempt_witness(res.x)
            ref=rational_proposal(res.x)
            if ref is None or any(dot(a,ref)<=0 for a in d['A']):ref=seed
            bound,lam=certify_tangent(d,box,ref,verts,s)
            node['certificates'][str(s)]={'reference':ref,'upper_bound':bound,'dual_lambda':lam,'proposal_success':bool(res.success)}
            if bound is None:
                node['status']='excluded_by_likelihood';return idx
            upper[s]=bound
        active[idx]=upper
        return idx

    process(root)
    while active:
        gaps={s:max([incumbent[s][0]]+[v[s] for v in active.values()])-incumbent[s][0] for s in (-1,1)}
        if max(gaps.values())<=tol:break
        if len(nodes)+2>max_nodes:break
        idx=max(active,key=lambda i:max(active[i][s]-incumbent[s][0] for s in (-1,1)))
        box=nodes[idx]['box']
        cand=[j for j,(l,u) in enumerate(box) if l<u]
        if not cand:break
        j=max(cand,key=lambda j:float(F(d['stratum_counts'][j],N)*(box[j][1]-box[j][0])**2/box[j][0]**2))
        l,u=box[j];mid=(l+u)/2;b1=list(box);b2=list(box);b1[j]=(l,mid);b2[j]=(mid,u)
        del active[idx];left=process(tuple(b1));right=process(tuple(b2))
        nodes[idx]['status']='split';nodes[idx]['split_axis']=j;nodes[idx]['children']=[left,right]
    upper={s:max([incumbent[s][0]]+[v[s] for v in active.values()]) for s in (-1,1)}
    result={'schema':'stage18_denominator_boxes_v1','input':data,'input_sha256':digest(data),
            'root_box':root,'tolerance':tol,'max_nodes':max_nodes,'nodes':nodes,
            'leaf_ids':sorted(active),'witnesses':{str(s):{'mass':incumbent[s][1],'value':incumbent[s][0],
                                                     'likelihood':likelihood(d,incumbent[s][1])} for s in (-1,1)},
            'lower_bracket':[-upper[-1],-incumbent[-1][0]],
            'upper_bracket':[incumbent[1][0],upper[1]],
            'precision_reached':max(upper[s]-incumbent[s][0] for s in (-1,1))<=tol,
            'numerical_proposals':numerical_calls,
            'proof':'concave-numerator tangents plus denominator secants, exact linear vertex bound'}
    return result
