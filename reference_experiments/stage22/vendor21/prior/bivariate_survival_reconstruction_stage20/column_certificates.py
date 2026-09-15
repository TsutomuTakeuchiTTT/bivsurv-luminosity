"""Stage 20: rational column-dual bounds without polytope vertex enumeration.

The likelihood, tail class and full-space compilation are inherited unchanged.
SLSQP/HiGHS propose points and multipliers. Every retained bound is verified
with Fraction inequalities on ALL original columns, not solver status.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import sys, math, json, time
sys.path.insert(0,str(Path(__file__).resolve().parent/'vendor19'))
import numpy as np
from scipy.optimize import minimize, linprog
from exact import dot, logint, encode
from multilimit import prepare, likelihood, rational_proposal, digest


def affine_upper(d, box, reference):
    """Return rational beta,z: ell(m)/N <= beta+z.m on the box.

    Numerator log tangents use upper log bounds. Denominator secants use
    minus lower log bounds at BOTH endpoints; their interpolation weights
    are nonnegative on the box. Coefficients themselves need not be positive.
    """
    H=len(d['c']);N=d['N'];coef=[F(0)]*H;beta=F(0)
    for a,n in zip(d['A'],d['counts']):
        av=dot(a,reference)
        if av<=0:raise ValueError('Tangent reference misses observed record')
        w=F(n,N);beta+=w*(logint(av)[1]-1)
        for h in range(H):
            if a[h]:coef[h]+=w/av
    for v,n,(l,u) in zip(d['V'],d['stratum_counts'],box):
        w=F(n,N)
        if l==u:beta-=w*logint(l)[0]
        else:
            sl=-logint(l)[0];su=-logint(u)[0];slope=(su-sl)/(u-l)
            beta+=w*(sl-slope*l)
            for h in range(H):
                if v[h]:coef[h]+=w*slope
    return beta,tuple(coef)


def record_thresholds(d):
    """Binary KL consequences of the SAME total likelihood threshold.
    e^{-B}=tau/L_sat is exact; no per-stratum confidence re-calibration.
    """
    sat=F(1)
    for n,j in zip(d['counts'],d['row_stratum']):sat*=F(n,d['stratum_counts'][j])**n
    ratio=d['tau']/sat
    if not 0<ratio<=1:raise ValueError('No nonempty confidence set under saturation bound')
    out=[]
    for k,j in zip(d['counts'],d['row_stratum']):
        n=d['stratum_counts'][j];a=F(k,n)
        out.append((k,n,ratio*(a**k)*((1-a)**(n-k))))
    return out


def make_record_envelopes(d,steps=24):
    result=[]
    for k,n,t in record_thresholds(d):
        f=lambda x:x**k*(1-x)**(n-k)
        mode=F(k,n);lo=F(0);hi=mode
        for _ in range(steps):
            mid=(lo+hi)/2
            if f(mid)>=t:hi=mid
            else:lo=mid
        lower=lo
        if k==n:upper=F(1)
        else:
            lo=mode;hi=F(1)
            for _ in range(steps):
                mid=(lo+hi)/2
                if f(mid)>=t:lo=mid
                else:hi=mid
            upper=hi
        result.append((lower,upper))
    return tuple(result)


def verify_record_envelopes(d,envelopes):
    ee=tuple(tuple(map(F,x)) for x in envelopes);tt=record_thresholds(d)
    assert len(ee)==len(tt)
    for (lo,hi),(k,n,t) in zip(ee,tt):
        mode=F(k,n);assert 0<=lo<=mode<=hi<=1
        assert lo==0 or lo**k*(1-lo)**(n-k)<=t
        assert hi==1 or hi**k*(1-hi)**(n-k)<=t
    return ee


def linear_constraints(d,box,reference=None):
    """G m <= b, m>=0, sum m=1. No numerical probability floors here."""
    G=[];b=[]
    for v,(l,u) in zip(d['V'],box):
        G.extend((tuple(map(F,v)),tuple(-F(x) for x in v)));b.extend((u,-l))
    for a,j,(lo,hi) in zip(d['A'],d['row_stratum'],d.get('record_envelopes',())):
        v=d['V'][j]
        G.append(tuple(lo*v[h]-a[h] for h in range(len(a))));b.append(F(0))
        G.append(tuple(a[h]-hi*v[h] for h in range(len(a))));b.append(F(0))
    affine=None
    if reference is not None:
        beta,t=affine_upper(d,box,reference)
        threshold_lo=logint(d['tau'])[0]/d['N']
        G.append(tuple(-x for x in t));b.append(beta-threshold_lo)
        affine={'constant':beta,'coefficients':t,'threshold_lower':threshold_lo}
    return tuple(G),tuple(b),affine


def rational_dual(G,b,objective,multipliers):
    """Construct an exact feasible dual from any nonnegative proposal.

    Projection to nonnegative multipliers is just a dual proposal, not a
    modification of the parent feasible set. Normalization multiplier is
    chosen by checking ALL columns, automatically repairing dual residuals.
    """
    if any(not math.isfinite(float(t)) for t in multipliers):raise ValueError('Nonfinite dual proposal')
    y=tuple(F(format(max(0.,float(t)),'.15g')) for t in multipliers)
    H=len(objective)
    residual=tuple(F(objective[h])-sum((y[i]*G[i][h] for i in range(len(G))),F(0)) for h in range(H))
    z=max(residual)
    upper=z+dot(b,y)
    return {'multipliers':y,'normalization_multiplier':z,'upper_bound':upper}


def lp_certificate(d,box,reference,sign):
    G,b,aff=linear_constraints(d,box,reference);H=len(d['c']);obj=tuple(sign*x for x in d['c'])
    gf=np.asarray(G,float);bf=np.asarray(b,float)
    # Scaling used ONLY by the proposal solver; map duals back before proof.
    scale=np.maximum(1.,np.maximum(np.max(np.abs(gf),axis=1),np.abs(bf)))
    res=linprog(-np.asarray(obj,float),A_ub=gf/scale[:,None],b_ub=bf/scale,
                A_eq=np.ones((1,H)),b_eq=[1.],bounds=(0,None),method='highs')
    if res.success:
        cert=rational_dual(G,b,obj,-res.ineqlin.marginals/scale)
        cert.update({'reference':reference,'kind':'upper','affine':aff,'proposal_status':int(res.status)})
        return cert
    # A solver's infeasible flag is NOT an exclusion certificate.
    # Search a normalized Farkas combination; verify its negativity exactly.
    rows=len(G)
    aux=linprog(np.r_[bf,1.],A_ub=np.vstack([np.c_[-gf.T,-np.ones(H)],np.r_[np.ones(rows),0.]]),
                b_ub=np.r_[np.zeros(H),1.],bounds=[(0,None)]*rows+[(None,None)],method='highs')
    if aux.success:
        cert=rational_dual(G,b,[0]*H,aux.x[:-1])
        if cert['upper_bound']<0:
            cert.update({'reference':reference,'kind':'empty','affine':aff,'proposal_status':int(res.status)})
            return cert
    cert=rational_dual(G,b,obj,[0]*len(G))
    cert.update({'reference':reference,'kind':'upper','affine':aff,'proposal_status':int(res.status),'fallback':True})
    return cert


def check_dual(d,box,sign,cert):
    ref=None if cert['reference'] is None else tuple(map(F,cert['reference']))
    if ref is not None:
        assert len(ref)==len(d['c']) and min(ref)>=0 and sum(ref)==1
    G,b,aff=linear_constraints(d,box,ref)
    assert encode(aff)==encode(cert['affine'])
    y=tuple(map(F,cert['multipliers']));z=F(cert['normalization_multiplier']);U=F(cert['upper_bound'])
    assert len(y)==len(G) and min(y)>=0
    obj=[0]*len(d['c']) if cert['kind']=='empty' else [sign*x for x in d['c']]
    assert all(z+sum((y[i]*G[i][h] for i in range(len(G))),F(0))>=obj[h] for h in range(len(obj)))
    assert z+dot(y,b)==U
    if cert['kind']=='empty':assert U<0
    else:assert cert['kind']=='upper'
    return len(obj)


def project_columns(data,tolerance=F(1,1000),max_nodes=201,proposal_maxiter=180):
    """Box search. Uncertified numerical failures only leave wider bounds."""
    tic=time.perf_counter();d=prepare(data);H=len(d['c']);J=len(d['V']);N=d['N'];tol=F(tolerance)
    d['record_envelopes']=make_record_envelopes(d)
    if tol<=0 or max_nodes<1:raise ValueError('Positive precision/budget required')
    A=np.asarray(d['A'],float);V=np.asarray(d['V'],float);w=np.asarray(d['counts'],float)/N;nj=np.asarray(d['stratum_counts'],float)/N
    c=np.asarray(d['c'],float);tl,tu=logint(d['tau']);threshold=float((tl+tu)/(2*N))
    root=tuple((F(1) if all(v) else 1-e,F(1)) for v,e in zip(d['V'],d['epsilon']))
    if any(l<=0 or u<l for l,u in root):raise ValueError('Explicit positive tail lower bound needed')
    m0=tuple(d['initial_mass']);L0=likelihood(d,m0)
    if L0 is None or L0<d['tau']:raise ValueError('Initial mass not a certified inner witness')
    incumbent={s:(s*dot(d['c'],m0),m0) for s in (-1,1)};nodes=[];active={};calls=0;phase_calls=0

    def eval_f(m,box=None):
        a=A@m;q=V@m
        if min(a)<=0 or min(q)<=0:return -1e70,np.zeros(H)
        val=float(w@np.log(a));grad=(w/a)@A
        if box is None:return val-float(nj@np.log(q)),grad-(nj/q)@V
        for j,(ll,uu) in enumerate(box):
            l,u=float(ll),float(uu)
            if l==u:val-=nj[j]*math.log(l)
            else:
                sl=-math.log(l);slope=(-math.log(u)-sl)/(u-l)
                val+=nj[j]*(sl+slope*(q[j]-l));grad+=nj[j]*slope*V[j]
        return val,grad

    def cons(box,with_likelihood=False,original=False):
        lo=np.asarray([float(l) for l,u in box]);hi=np.asarray([float(u) for l,u in box]);ids=[j for j,v in enumerate(V) if not np.all(v==1)]
        cc=[{'type':'eq','fun':lambda m:float(m.sum()-1),'jac':lambda m:np.ones(H)},
            {'type':'ineq','fun':lambda m:A@m-1e-13,'jac':lambda m:A}]
        if ids:
            vv=V[ids];l=lo[ids];u=hi[ids]
            cc += [{'type':'ineq','fun':lambda m:vv@m-l,'jac':lambda m:vv},
                   {'type':'ineq','fun':lambda m:u-vv@m,'jac':lambda m:-vv}]
        if with_likelihood:
            b0=None if original else box
            cc += [{'type':'ineq','fun':lambda m:eval_f(m,b0)[0]-threshold,'jac':lambda m:eval_f(m,b0)[1]}]
        return cc

    def attempt(z):
        p=rational_proposal(z)
        if p is None:return
        for e in (F(0),F(1,10**10),F(1,10**8),F(1,10**6),F(1,10000),F(1,100),F(1,2)):
            m=tuple((1-e)*a+e*b for a,b in zip(p,m0));L=likelihood(d,m)
            if L is not None and L>=d['tau']:
                for s in (-1,1):
                    val=s*dot(d['c'],m)
                    if val>incumbent[s][0]:incumbent[s]=(val,m)
                return

    for s in (-1,1):
        res=minimize(lambda m:-s*float(c@m),np.asarray(m0,float),jac=lambda m:-s*c,
                     bounds=[(0,1)]*H,constraints=cons(root,True,True),method='SLSQP',
                     options={'ftol':1e-12,'maxiter':400})
        calls+=1;attempt(res.x)

    def process(box):
        nonlocal calls,phase_calls
        idx=len(nodes);node={'id':idx,'box':box,'status':'leaf','certificates':{}};nodes.append(node)
        # LP maximin observed numerator supplies a numerical starting point.
        G,b,_=linear_constraints(d,box,None)
        Gf=np.asarray(G,float);bf=np.asarray(b,float)
        phase=linprog(np.r_[np.zeros(H),-1.],A_ub=np.vstack([np.c_[Gf,np.zeros(len(G))],np.c_[-A,np.ones(len(A))]]),
                      b_ub=np.r_[bf,np.zeros(len(A))],A_eq=np.r_[np.ones(H),0.][None,:],b_eq=[1.],bounds=[(0,None)]*(H+1),method='highs')
        phase_calls+=1
        if phase.success and phase.x[-1]>1e-14:z=phase.x[:H]
        else:z=np.asarray(m0,float)
        refs=[m0]
        rp=rational_proposal(z)
        if rp is not None and all(dot(a,rp)>0 for a in d['A']):refs.append(rp)
        # Always a globally admissible fallback reference, even if outside box.
        if eval_f(z,box)[0]<threshold:
            res=minimize(lambda m:-eval_f(m,box)[0],z,jac=lambda m:-eval_f(m,box)[1],
                         bounds=[(0,1)]*H,constraints=cons(box),method='SLSQP',
                         options={'ftol':1e-12,'maxiter':proposal_maxiter})
            calls+=1;z=res.x
            rp=rational_proposal(z)
            if rp is not None and all(dot(a,rp)>0 for a in d['A']):refs.append(rp)
        upper={}
        for s in (-1,1):
            wi=incumbent[s][1]
            zi=np.asarray(wi,float) if all(l<=dot(v,wi)<=u for v,(l,u) in zip(d['V'],box)) else z
            res=minimize(lambda m:-s*float(c@m),zi,jac=lambda m:-s*c,
                         bounds=[(0,1)]*H,constraints=cons(box,True),method='SLSQP',
                         options={'ftol':1e-12,'maxiter':proposal_maxiter})
            calls+=1;attempt(res.x)
            ref=rational_proposal(res.x)
            if ref is None or any(dot(a,ref)<=0 for a in d['A']):ref=m0
            refs.append(ref)
            cert=None
            for rr in dict.fromkeys(refs):
                candidate=lp_certificate(d,box,rr,s)
                if candidate['kind']=='empty':cert=candidate;break
                if cert is None or candidate['upper_bound']<cert['upper_bound']:cert=candidate
            check_dual(d,box,s,cert)
            node['certificates'][str(s)]=cert
            if cert['kind']=='empty':node['status']='excluded';return idx
            upper[s]=cert['upper_bound']
        active[idx]=upper
        return idx

    process(root)
    while active:
        ub={s:max([incumbent[s][0]]+[a[s] for a in active.values()]) for s in (-1,1)}
        if max(ub[s]-incumbent[s][0] for s in (-1,1))<=tol or len(nodes)+2>max_nodes:break
        idx=max(active,key=lambda i:max(active[i][s]-incumbent[s][0] for s in (-1,1)))
        box=nodes[idx]['box'];choices=[j for j,(l,u) in enumerate(box) if l<u]
        if not choices:break
        j=max(choices,key=lambda j:float(F(d['stratum_counts'][j],N)*(box[j][1]-box[j][0])**2/box[j][0]**2))
        l,u=box[j];mid=(l+u)/2;left=list(box);right=list(box);left[j]=(l,mid);right[j]=(mid,u)
        del active[idx];i1=process(tuple(left));i2=process(tuple(right))
        nodes[idx].update({'status':'split','axis':j,'children':[i1,i2]})
    ub={s:max([incumbent[s][0]]+[a[s] for a in active.values()]) for s in (-1,1)}
    return {'schema':'stage20_column_dual_v1','input':data,'input_sha256':digest(data),'root_box':root,
            'nodes':nodes,'record_envelopes':d['record_envelopes'],'leaf_ids':sorted(active),'witnesses':{str(s):{'mass':incumbent[s][1],'value':incumbent[s][0],
              'likelihood':likelihood(d,incumbent[s][1])} for s in (-1,1)},
            'lower_bracket':[-ub[-1],-incumbent[-1][0]],'upper_bracket':[incumbent[1][0],ub[1]],
            'tolerance':tol,'max_nodes':max_nodes,'precision_reached':max(ub[s]-incumbent[s][0] for s in (-1,1))<=tol,
            'nonlinear_proposals':calls,'phase_LP_proposals':phase_calls,'vertex_enumeration_called':False,
            'seconds':time.perf_counter()-tic}


def verify_document(doc,rebuild=True):
    """No optimization, no root search, no vertex enumeration. Exact checks."""
    data=doc['input'];assert digest(data)==doc['input_sha256'];d=prepare(data)
    d['record_envelopes']=verify_record_envelopes(d,doc['record_envelopes'])
    if rebuild:
        from integrated import compile_input
        rebuilt,_,_,_=compile_input(data['design'],data['target_index'])
        assert digest(rebuilt)==doc['input_sha256'],'Compiled model/predictor mismatch'
    root=tuple((F(1) if all(v) else 1-e,F(1)) for v,e in zip(d['V'],d['epsilon']))
    assert root==tuple(tuple(map(F,x)) for x in doc['root_box'])
    nodes=doc['nodes'];seen=set();leaves={};checks=0;certs=0;excluded=0
    def walk(idx,box):
        nonlocal checks,certs,excluded
        assert idx not in seen and 0<=idx<len(nodes);seen.add(idx);node=nodes[idx]
        assert idx==node['id'] and box==tuple(tuple(map(F,x)) for x in node['box'])
        if node['status']=='split':
            j=node['axis'];l,u=box[j];assert l<u;mid=(l+u)/2;x=list(box);y=list(box);x[j]=(l,mid);y[j]=(mid,u)
            assert len(node['children'])==2;walk(node['children'][0],tuple(x));walk(node['children'][1],tuple(y));return
        cc=node['certificates']
        if node['status']=='excluded':
            assert any(v['kind']=='empty' for v in cc.values())
            for sign,c0 in cc.items():checks+=check_dual(d,box,int(sign),c0);certs+=1
            excluded+=1;return
        assert node['status']=='leaf' and set(cc)=={'-1','1'}
        ub={}
        for sign,c0 in cc.items():
            assert c0['kind']=='upper';checks+=check_dual(d,box,int(sign),c0);certs+=1;ub[int(sign)]=F(c0['upper_bound'])
        leaves[idx]=ub
    walk(0,root);assert len(seen)==len(nodes) and sorted(leaves)==doc['leaf_ids']
    inc={}
    for sign,w in doc['witnesses'].items():
        sign=int(sign);m=tuple(map(F,w['mass']));L=likelihood(d,m)
        assert L is not None and L>=d['tau'] and L==F(w['likelihood'])
        inc[sign]=sign*dot(d['c'],m);assert inc[sign]==F(w['value'])
    ub={s:max([inc[s]]+[a[s] for a in leaves.values()]) for s in (-1,1)}
    assert list(map(F,doc['lower_bracket']))==[-ub[-1],-inc[-1]]
    assert list(map(F,doc['upper_bracket']))==[inc[1],ub[1]]
    assert doc['precision_reached']==(max(ub[s]-inc[s] for s in (-1,1))<=F(doc['tolerance']))
    return {'passed':True,'nodes':len(nodes),'column_inequalities':checks,'leaf_certificates':certs,
            'excluded_leaves':excluded,'parent_witnesses':2,'precision_reached':doc['precision_reached'],
            'optimization_called':False,'vertex_enumeration_called':False}
