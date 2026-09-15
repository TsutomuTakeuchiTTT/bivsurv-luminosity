"""Stage21: reparametrized convex relaxations and retained column-dual proofs.

Only proposal calculations use floating point. The model, exact likelihood
threshold, original mass simplex and tail restrictions are unchanged.
The 'relative' mode normalizes OBSERVATION expressions by union visibility;
completely invisible parent mass is NOT deleted from variables or target.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import sys,math,json,time
sys.set_int_max_str_digits(0)
PRIOR=Path(__file__).parent/'prior/bivariate_survival_reconstruction_stage20'
sys.path.insert(0,str(PRIOR));sys.path.insert(0,str(PRIOR/'vendor19'))
import numpy as np
from scipy.optimize import linprog,minimize
from exact import dot,logint,encode
from multilimit import prepare,likelihood,rational_proposal,digest
from column_certificates import verify_record_envelopes


def factorized_threshold(data):
    """Log(tau) without passing a many-thousand-digit product into logint."""
    a=F(data['alpha']);t=a;lo,hi=logint(a)
    for probs,ns in zip(data['prediction'],data['validation']):
        for p,n in zip(probs,ns):
            if n:
                p=F(p);assert p>0;t*=p**n;l,u=logint(p);lo+=n*l;hi+=n*u
    assert t==F(data['tau'])
    return lo,hi


class Refinement:
    def __init__(self,prior_doc,mode='relative'):
        self.old=prior_doc;self.data=prior_doc['input'];self.d=prepare(self.data)
        self.d['record_envelopes']=verify_record_envelopes(self.d,prior_doc['record_envelopes'])
        d=self.d;self.H=len(d['c']);self.J=len(d['V']);self.N=d['N'];self.mode=mode
        assert mode in ('raw','relative')
        self.star=tuple(max(v[h] for v in d['V']) for h in range(self.H))
        self.norm=tuple([1]*self.H) if mode=='raw' else self.star
        self.root=tuple((F(1),F(1)) if tuple(v)==self.norm else (1-e,F(1)) for v,e in zip(d['V'],d['epsilon']))
        self.tlo,self.thi=(x/self.N for x in factorized_threshold(self.data))
        self.threshold=float((self.tlo+self.thi)/2)
        self.A=np.array(d['A'],float);self.V=np.array(d['V'],float);self.c=np.array(d['c'],float);self.vnorm=np.array(self.norm,float)
        self.w=np.array(d['counts'],float)/self.N;self.nj=np.array(d['stratum_counts'],float)/self.N
        self.m0=tuple(d['initial_mass']);self.calls=0;self.lpcalls=0;self.refcache={};self.traces=[]
        self.inc={s:(F(prior_doc['witnesses'][str(s)]['value']),tuple(map(F,prior_doc['witnesses'][str(s)]['mass']))) for s in (-1,1)}
        assert all(likelihood(d,m)>=d['tau'] for val,m in self.inc.values())
        self.oldouter=(F(prior_doc['lower_bracket'][0]),F(prior_doc['upper_bracket'][1]))
        self.nodes=[];self.active={};self.globalrefs=[]

    def base(self,box):
        G=[];b=[]
        for v,e,(l,u) in zip(self.d['V'],self.d['epsilon'],box):
            G.append(tuple(-F(x) for x in v));b.append(e-1)
            G.append(tuple(F(vh)-u*vh0 for vh,vh0 in zip(v,self.norm)));b.append(F(0))
            G.append(tuple(l*vh0-F(vh) for vh,vh0 in zip(v,self.norm)));b.append(F(0))
        for a,j,(l,u) in zip(self.d['A'],self.d['row_stratum'],self.d['record_envelopes']):
            v=self.d['V'][j]
            G.extend([tuple(l*vh-ah for vh,ah in zip(v,a)),tuple(ah-u*vh for vh,ah in zip(v,a))]);b.extend([F(0),F(0)])
        # Previously certified GLOBAL target bounds, valid for every child box.
        G.extend([tuple(-x for x in self.d['c']),self.d['c']]);b.extend([-self.oldouter[0],self.oldouter[1]])
        return G,b

    def numerator(self,ref):
        ref=tuple(ref)
        if ref in self.refcache:return self.refcache[ref]
        q0=dot(self.norm,ref);assert q0>0
        coef=[F(0)]*self.H;beta=F(0)
        for a,n in zip(self.d['A'],self.d['counts']):
            p=dot(a,ref)/q0;assert p>0;w=F(n,self.N)
            beta+=w*(logint(p)[1]-1)
            for h,ah in enumerate(a):
                if ah:coef[h]+=w/p
        self.refcache[ref]=(beta,tuple(coef));return self.refcache[ref]

    def cut(self,box,ref):
        beta,t=self.numerator(ref);t=list(t)
        for v,n,(l,u) in zip(self.d['V'],self.d['stratum_counts'],box):
            w=F(n,self.N)
            if l==u:beta-=w*logint(l)[0]
            else:
                sl=-logint(l)[0];su=-logint(u)[0];slope=(su-sl)/(u-l)
                beta+=w*(sl-slope*l)
                for h,vh in enumerate(v):
                    if vh:t[h]+=w*slope
        # Multiply the affine necessary inequality in m/qstar by qstar.
        return tuple((self.tlo-beta)*vh-th for vh,th in zip(self.norm,t))

    def dual(self,G,b,sign,empty=False):
        gf=np.array(G,float);bf=np.array(b,float);obj=np.zeros(self.H) if empty else sign*self.c
        scale=np.maximum(1.,np.maximum(np.max(abs(gf),axis=1),abs(bf)))
        res=linprog(-obj,A_ub=gf/scale[:,None],b_ub=bf/scale,A_eq=np.ones((1,self.H)),b_eq=[1.],bounds=(0,None),method='highs')
        self.lpcalls+=1
        if res.success and not empty:yy=-res.ineqlin.marginals/scale;kind='upper'
        else:
            R=len(G)
            aux=linprog(np.r_[bf,1.],A_ub=np.vstack((np.c_[-gf.T,-np.ones(self.H)],np.r_[np.ones(R),0.])),b_ub=np.r_[np.zeros(self.H),1.],bounds=[(0,None)]*R+[(None,None)],method='highs');self.lpcalls+=1
            if aux.success:
                yy=aux.x[:-1];kind='empty';obj=np.zeros(self.H)
            else:yy=np.zeros(len(G));kind='upper';obj=sign*self.c
        y=tuple(F(format(max(0.,float(x)),'.13g')) if x>1e-10 else F(0) for x in yy)
        objective=tuple([F(0)]*self.H) if kind=='empty' else tuple(sign*x for x in self.d['c'])
        z=max(objective[h]-sum((v*g[h] for v,g in zip(y,G) if v),F(0)) for h in range(self.H))
        U=z+dot(y,b)
        if kind=='empty' and U>=0:
            kind='upper';y=tuple(F(0) for _ in G);z=max(sign*x for x in self.d['c']);U=z
        cert={'multipliers':y,'normalization_multiplier':z,'upper_bound':U,'kind':kind}
        return cert,res

    def original(self,m):
        a=self.A@m;q=self.V@m
        if min(a)<=0 or min(q)<=0:return -1e60,np.zeros(self.H)
        return float(self.w@np.log(a)-self.nj@np.log(q)),(self.w/a)@self.A-(self.nj/q)@self.V

    def secant_float(self,box):
        beta=0.;t=np.zeros(self.H)
        for v,n,(l,u) in zip(self.V,self.nj,box):
            l=float(l);u=float(u)
            if l==u:beta-=n*math.log(l)
            else:
                sl=-math.log(l);slope=(-math.log(u)-sl)/(u-l);beta+=n*(sl-slope*l);t+=n*slope*v
        return beta,t

    def perspective(self,m,box):
        a=self.A@m;q=float(self.vnorm@m);db,dt=self.secant_float(box)
        if min(a)<=0 or q<=0:return -1e60,np.zeros(self.H)
        v=float(self.w@np.log(a/q)+db-self.threshold)
        return q*v+float(dt@m),self.vnorm*(v-1)+q*(self.w/a)@self.A+dt

    def as_reference(self,z):
        m=rational_proposal(z)
        if m is None:return None
        if min(dot(a,m) for a in self.d['A'])<=0:
            e=F(1,10**6);m=tuple((1-e)*x+e*y for x,y in zip(m,self.m0))
        return m

    def attempt(self,z):
        p=rational_proposal(z)
        if p is None:return
        # Mixtures propose a new feasible point; only ORIGINAL exact L accepts.
        for e in [F(0),F(1,10**9),F(1,10**7),F(1,10**5),F(1,1000),F(1,100),F(1,2)]:
            m=tuple((1-e)*x+e*y for x,y in zip(p,self.m0))
            improvements=[s for s in (-1,1) if s*dot(self.d['c'],m)>self.inc[s][0]]
            if not improvements:return
            L=likelihood(self.d,m)
            if L is not None and L>=self.d['tau']:
                for s in improvements:self.inc[s]=(s*dot(self.d['c'],m),m)
                return

    def inner_search(self,z,s,maxiter=300):
        H=self.H;A=self.A;V=self.V;lo=np.array([float(1-e) for e in self.d['epsilon']])
        cc=[{'type':'eq','fun':lambda m:m.sum()-1,'jac':lambda m:np.ones(H)}, {'type':'ineq','fun':lambda m:V@m-lo,'jac':lambda m:V}, {'type':'ineq','fun':lambda m:A@m-1e-13,'jac':lambda m:A}, {'type':'ineq','fun':lambda m:self.original(m)[0]-self.threshold,'jac':lambda m:self.original(m)[1]}]
        r=minimize(lambda m:-s*self.c@m,np.array(z,float),jac=lambda m:-s*self.c,bounds=[(0.,1.)]*H,constraints=cc,method='SLSQP',options={'ftol':1e-12,'maxiter':maxiter});self.calls+=1;self.attempt(r.x)
        return r

    def process(self,box,tol,maxiter=250):
        idx=len(self.nodes);node={'id':idx,'box':box,'status':'leaf','certificates':{}};self.nodes.append(node)
        G,b=self.base(box);gf=np.array(G,float);bf=np.array(b,float);H=self.H
        # Positive phase point for convex relaxed optimization; not a proof.
        phase=linprog(np.r_[np.zeros(H),-1.],A_ub=np.vstack((np.c_[gf,np.zeros(len(G))],np.c_[-self.A,np.ones(len(self.A))])),b_ub=np.r_[bf,np.zeros(len(self.A))],A_eq=np.r_[np.ones(H),0.][None,:],b_eq=[1.],bounds=[(0,None)]*(H+1),method='highs');self.lpcalls+=1
        refs=[self.m0]+[m for _,m in self.inc.values()]
        z=phase.x[:H] if phase.success and phase.x[-1]>1e-13 else np.array(self.m0,float)
        rr=self.as_reference(z)
        if rr is not None:refs.append(rr)
        mask=np.max(abs(gf),axis=1)>0;gf2=gf[mask];bf2=bf[mask]
        cc=[{'type':'eq','fun':lambda m:m.sum()-1,'jac':lambda m:np.ones(H)}, {'type':'ineq','fun':lambda m:bf2-gf2@m,'jac':lambda m:-gf2}, {'type':'ineq','fun':lambda m:self.A@m-1e-13,'jac':lambda m:self.A}]
        if not phase.success or self.perspective(z,box)[0]<0:
            rr=minimize(lambda m:-self.perspective(m,box)[0],z,jac=lambda m:-self.perspective(m,box)[1],bounds=[(0,1)]*H,constraints=cc,method='SLSQP',options={'ftol':1e-12,'maxiter':maxiter});self.calls+=1;z=rr.x
            rp=self.as_reference(z)
            if rp is not None:refs.append(rp)
        upper={};proposal_notes=[]
        for s in (-1,1):
            w0=np.array(self.inc[s][1],float)
            start=w0 if min(bf-gf@w0)>=-1e-10 else z
            cx=cc+[{'type':'ineq','fun':lambda m:self.perspective(m,box)[0],'jac':lambda m:self.perspective(m,box)[1]}]
            res=minimize(lambda m:-s*self.c@m,start,jac=lambda m:-s*self.c,bounds=[(0,1)]*H,constraints=cx,method='SLSQP',options={'ftol':1e-12,'maxiter':maxiter});self.calls+=1
            rp=self.as_reference(res.x)
            if rp is not None:refs.append(rp)
            self.attempt(res.x)
            if res.success and s*self.c@res.x>float(self.inc[s][0])+float(tol)/5:
                self.inner_search(res.x,s,min(350,maxiter))
            refs=list(dict.fromkeys(refs))
            # Jointly retain valid tangents, not the best of single-cut bounds.
            GG=G+[self.cut(box,r) for r in refs];bb=b+[F(0)]*len(refs)
            cert,lp=self.dual(GG,bb,s)
            # A few cutting-plane refinements at LP optimizers, strictly proposal-only.
            for it in range(6):
                if cert['kind']=='empty' or cert['upper_bound']-self.inc[s][0]<=tol:break
                if not lp.success:break
                rp=self.as_reference(lp.x)
                if rp is None or rp in refs:break
                refs.append(rp);GG.append(self.cut(box,rp));bb.append(F(0))
                candidate,lp=self.dual(GG,bb,s)
                if candidate['kind']=='empty' or candidate['upper_bound']<=cert['upper_bound']:cert=candidate
                else:break
            # If a previous best used fewer rows, zero-pad its multipliers.
            yy=list(cert['multipliers']);yy += [F(0)]*(len(GG)-len(yy))
            keep=[i for i in range(len(refs)) if yy[len(G)+i]!=0]
            cert['references']=[refs[i] for i in keep]
            cert['multipliers']=tuple(yy[:len(G)]+[yy[len(G)+i] for i in keep])
            node['certificates'][str(s)]=cert
            proposal_notes.append({'sign':s,'slsqp_success':bool(res.success),'iterations':int(res.nit),'proposed_objective':float(s*self.c@res.x),'relaxed_residual':self.perspective(res.x,box)[0],'original_log_margin':self.original(res.x)[0]-self.threshold})
            if cert['kind']=='empty':node['status']='excluded';node['proposal_diagnostics']=proposal_notes;return idx
            upper[s]=cert['upper_bound']
        node['proposal_diagnostics']=proposal_notes;self.active[idx]=upper
        return idx

    def run(self,tolerance=F(1,1000),max_nodes=401,maxiter=250):
        tic=time.perf_counter();tol=F(tolerance)
        for s in (-1,1):self.inner_search(self.inc[s][1],s,600)
        self.process(self.root,tol,maxiter)
        while self.active:
            ub={s:max([self.inc[s][0]]+[v[s] for v in self.active.values()]) for s in (-1,1)}
            gap=max(ub[s]-self.inc[s][0] for s in (-1,1))
            if gap<=tol or len(self.nodes)+2>max_nodes:break
            idx=max(self.active,key=lambda i:max(self.active[i][s]-self.inc[s][0] for s in (-1,1)))
            box=self.nodes[idx]['box'];ids=[j for j,(l,u) in enumerate(box) if l<u]
            if not ids:break
            j=max(ids,key=lambda j:float(F(self.d['stratum_counts'][j],self.N)*(box[j][1]-box[j][0])**2/box[j][0]**2))
            l,u=box[j];mid=(l+u)/2;left=list(box);right=list(box);left[j]=(l,mid);right[j]=(mid,u)
            del self.active[idx];a=self.process(tuple(left),tol,maxiter);b=self.process(tuple(right),tol,maxiter)
            self.nodes[idx].update({'status':'split','axis':j,'children':[a,b]})
            if len(self.nodes)%20==1:
                print(self.data['name'],len(self.nodes),'gap',float(gap),'inc',[-float(self.inc[-1][0]),float(self.inc[1][0])],flush=True)
        ub={s:max([self.inc[s][0]]+[v[s] for v in self.active.values()]) for s in (-1,1)}
        return {'schema':'stage21_perspective_multicut_v1','input':self.data,'input_sha256':digest(self.data),'prior_input_sha256':self.old['input_sha256'],'prior_outer_bounds':self.oldouter,'mode':self.mode,'normalizer_row':self.norm,'root_box':self.root,'record_envelopes':self.d['record_envelopes'],'threshold_log_interval':[self.tlo,self.thi],'nodes':self.nodes,'leaf_ids':sorted(self.active),'witnesses':{str(s):{'mass':m,'value':val,'likelihood':likelihood(self.d,m)} for s,(val,m) in self.inc.items()},'lower_bracket':[-ub[-1],-self.inc[-1][0]],'upper_bracket':[self.inc[1][0],ub[1]],'precision_reached':max(ub[s]-self.inc[s][0] for s in (-1,1))<=tol,'tolerance':tol,'max_nodes':max_nodes,'maxiter':maxiter,'nonlinear_proposals':self.calls,'LP_proposals':self.lpcalls,'seconds':time.perf_counter()-tic}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('name');p.add_argument('--mode',choices=['raw','relative'],default='relative');p.add_argument('--nodes',type=int,default=401);p.add_argument('--maxiter',type=int,default=250);p.add_argument('--outdir',default=str(Path(__file__).parent/'outputs'));a=p.parse_args()
    old=json.loads((PRIOR/'outputs/certificates'/f'{a.name}.json').read_text());o=Refinement(old,a.mode);result=o.run(max_nodes=a.nodes,maxiter=a.maxiter)
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True);(out/f'{a.name}.json').write_text(json.dumps(encode(result),sort_keys=True,indent=2)+'\n')
    print(a.name,'DONE',result['precision_reached'],[float(x) for x in result['lower_bracket']],[float(x) for x in result['upper_bracket']],len(result['nodes']),result['seconds'],flush=True)
