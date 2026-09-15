"""Independent arithmetic verifier for Stage 20.

Does NOT call the search's affine/constraint/certificate routines, root search,
LP, nonlinear optimization, or vertex enumeration. It evaluates the log tangent
at each simplex basis column, shifts the cut to zero RHS using normalization,
and checks the resulting dual in that alternative representation.
The proven log enclosures and the previous full-space observation compiler are
shared arithmetic/model primitives, not numerical optimization oracles.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib,argparse
sys.path.insert(0,str(Path(__file__).resolve().parent/'vendor19'))
from exact import dot,logint,encode
from integrated import compile_input


def digest(x):return hashlib.sha256(json.dumps(encode(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()

def verify(doc):
    data=doc['input'];assert digest(data)==doc['input_sha256']
    rebuilt,*_=compile_input(data['design'],data['target_index']);assert digest(rebuilt)==doc['input_sha256']
    A=tuple(tuple(map(F,r)) for r in data['A']);V=tuple(tuple(map(F,r)) for r in data['V'])
    c=tuple(map(F,data['c']));ns=data['counts'];js=data['row_stratum'];N=sum(ns);H=len(c);J=len(V)
    nj=[sum(n for n,jj in zip(ns,js) if jj==j) for j in range(J)]
    tau=F(data['tau']);eps=tuple(map(F,data['epsilon']));sat=F(1)
    for n,j in zip(ns,js):sat*=F(n,nj[j])**n
    assert 0<tau<=sat
    rho=tau/sat;env=tuple(tuple(map(F,r)) for r in doc['record_envelopes']);assert len(env)==len(A)
    for (l,u),k,j in zip(env,ns,js):
        a=F(k,nj[j]);t=rho*a**k*(1-a)**(nj[j]-k)
        assert 0<=l<=a<=u<=1
        assert l==0 or l**k*(1-l)**(nj[j]-k)<=t
        assert u==1 or u**k*(1-u)**(nj[j]-k)<=t
    tlo=logint(tau)[0]/N
    root=tuple((F(1) if all(v) else 1-e,F(1)) for v,e in zip(V,eps))
    assert root==tuple(tuple(map(F,r)) for r in doc['root_box'])
    seen=set();leaves={};cols=certcount=emptycount=0

    def certificate(box,sign,cert):
        nonlocal cols,certcount
        G=[];b=[]
        for v,(l,u) in zip(V,box):G.extend((v,tuple(-x for x in v)));b.extend((u,-l))
        for a,j,(l,u) in zip(A,js,env):
            G.extend((tuple(l*x-y for x,y in zip(V[j],a)),tuple(y-u*x for x,y in zip(V[j],a))))
            b.extend((F(0),F(0)))
        ref=cert['reference'];shift=F(0)
        if ref is not None:
            ref=tuple(map(F,ref));assert len(ref)==H and min(ref)>=0 and sum(ref)==1
            av=[dot(a,ref) for a in A];assert min(av)>0
            lv=[logint(a)[1] for a in av]
            # Direct values T(e_h): distinct algebraic representation from builder.
            T=[]
            for h in range(H):
                value=sum((F(n,N)*(l+a[h]/p-1) for n,l,a,p in zip(ns,lv,A,av)),F(0))
                for v,n,(l,u) in zip(V,nj,box):
                    if l==u:term=-logint(l)[0]
                    else:term=-(u-v[h])/(u-l)*logint(l)[0]-(v[h]-l)/(u-l)*logint(u)[0]
                    value+=F(n,N)*term
                T.append(value)
            beta=F(cert['affine']['constant']);coef=tuple(map(F,cert['affine']['coefficients']))
            assert F(cert['affine']['threshold_lower'])==tlo
            assert tuple(beta+x for x in coef)==tuple(T)
            G.append(tuple(tlo-x for x in T));b.append(F(0))
            shift=beta-tlo
        else:assert cert['affine'] is None
        y=tuple(map(F,cert['multipliers']));assert len(y)==len(G) and min(y)>=0
        z=F(cert['normalization_multiplier'])+(y[-1]*shift if ref is not None else 0)
        objective=[0]*H if cert['kind']=='empty' else [sign*x for x in c]
        assert all(z+sum((yy*g[h] for yy,g in zip(y,G)),F(0))>=objective[h] for h in range(H))
        upper=z+dot(y,b);assert upper==F(cert['upper_bound'])
        if cert['kind']=='empty':assert upper<0
        else:assert cert['kind']=='upper'
        cols+=H;certcount+=1
        return upper

    def walk(idx,box):
        nonlocal emptycount
        assert idx not in seen and 0<=idx<len(doc['nodes']);seen.add(idx);node=doc['nodes'][idx]
        assert node['id']==idx and tuple(tuple(map(F,r)) for r in node['box'])==box
        if node['status']=='split':
            j=node['axis'];l,u=box[j];assert l<u;mid=(l+u)/2
            left=list(box);right=list(box);left[j]=(l,mid);right[j]=(mid,u)
            assert len(node['children'])==2;walk(node['children'][0],tuple(left));walk(node['children'][1],tuple(right));return
        cc=node['certificates']
        if node['status']=='excluded':
            assert any(z['kind']=='empty' for z in cc.values())
            for s,z in cc.items():certificate(box,int(s),z)
            emptycount+=1
        else:
            assert node['status']=='leaf' and set(cc)=={'-1','1'} and all(z['kind']=='upper' for z in cc.values())
            leaves[idx]={int(s):certificate(box,int(s),z) for s,z in cc.items()}
    walk(0,root);assert len(seen)==len(doc['nodes']) and sorted(leaves)==doc['leaf_ids']
    inc={}
    assert set(doc['witnesses'])=={'-1','1'}
    for s,w in doc['witnesses'].items():
        m=tuple(map(F,w['mass']));assert len(m)==H and min(m)>=0 and sum(m)==1
        q=[dot(v,m) for v in V];assert all(x>=1-e for x,e in zip(q,eps))
        L=F(1)
        for a,n,j in zip(A,ns,js):
            p=dot(a,m);assert p>0;L*=(p/q[j])**n
        assert L==F(w['likelihood']) and L>=tau
        inc[int(s)]=int(s)*dot(c,m);assert inc[int(s)]==F(w['value'])
    upper={s:max([inc[s]]+[w[s] for w in leaves.values()]) for s in (-1,1)}
    assert list(map(F,doc['lower_bracket']))==[-upper[-1],-inc[-1]]
    assert list(map(F,doc['upper_bracket']))==[inc[1],upper[1]]
    precision=max(upper[s]-inc[s] for s in (-1,1))<=F(doc['tolerance'])
    assert precision==doc['precision_reached']
    return {'passed':True,'nodes':len(seen),'column_inequalities':cols,'certificates':certcount,
            'excluded_leaves':emptycount,'parent_witnesses':2,'record_envelopes_checked':len(env),
            'precision_reached':precision,'optimization_or_root_search_called':False,'vertex_enumeration_called':False,
            'alternative_zero_rhs_affine_representation':True}


def main(outdir):
    outdir=Path(outdir)
    # Fail immediately if any solver is accidentally invoked by validation.
    import scipy.optimize,integrated,exact,multilimit
    def forbidden(*a,**kw):raise AssertionError('Verifier must not call optimization or enumerate vertices')
    scipy.optimize.linprog=scipy.optimize.minimize=forbidden
    integrated.minimize=multilimit.minimize=forbidden
    exact.vertices=multilimit.vertices=forbidden
    result={}
    for path in sorted((outdir/'certificates').glob('*.json')):
        tic=time.perf_counter();v=verify(json.loads(path.read_text()));v['verification_seconds']=time.perf_counter()-tic;result[path.stem]=v
        print(path.stem,v['passed'],v['nodes'],v['column_inequalities'],flush=True)
    answer={'all_passed':True,'cases':result,'solvers_disabled_during_verification':True}
    (outdir/'independent_column_verification.json').write_text(json.dumps(answer,sort_keys=True,indent=2)+'\n')
    return answer
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',default=str(Path(__file__).parent/'outputs'));a=p.parse_args();main(a.outdir)
