"""Independent Stage21 certificate check (no optimization or root search).

Reconstructs supporting-plane coefficients by direct homogeneous evaluation at
simplex basis vectors. Does not import refinement's cut, constraints or dual.
Stage20 ancestor certificates are rechecked when verify_prior=True.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib,argparse
sys.set_int_max_str_digits(0)
BASE=Path(__file__).parent
P=BASE/'prior/bivariate_survival_reconstruction_stage20'
sys.path.insert(0,str(P));sys.path.insert(0,str(P/'vendor19'))
from exact import dot,logint,encode
from integrated import compile_input
from verify_stage20 import verify as verify_old

def digest(x):return hashlib.sha256(json.dumps(encode(x),sort_keys=True,separators=(',',':')).encode()).hexdigest()

def verify(doc,old,verify_prior=True):
    if verify_prior:old_result=verify_old(old)
    else:old_result=None
    data=doc['input'];assert digest(data)==doc['input_sha256']==old['input_sha256']==doc['prior_input_sha256']
    reconstructed,*_=compile_input(data['design'],data['target_index']);assert digest(reconstructed)==digest(data)
    A=tuple(tuple(map(F,r)) for r in data['A']);V=tuple(tuple(map(F,r)) for r in data['V']);c=tuple(map(F,data['c']));eps=tuple(map(F,data['epsilon']))
    ns=data['counts'];js=data['row_stratum'];N=sum(ns);H=len(c);J=len(V);nj=[sum(n for n,jj in zip(ns,js) if j==jj) for j in range(J)]
    tau=F(data['alpha']);tlo,thi=logint(tau)
    for row,cnt in zip(data['prediction'],data['validation']):
        for p,n in zip(row,cnt):
            if n:
                p=F(p);assert p>0;ll,uu=logint(p);tlo+=n*ll;thi+=n*uu;tau*=p**n
    assert tau==F(data['tau']);tlo/=N;thi/=N
    assert tuple(map(F,doc['threshold_log_interval']))==(tlo,thi)
    # A universal union-visibility row, not an assumed deep observed stratum.
    star=tuple(max(v[h] for v in V) for h in range(H))
    norm=tuple([F(1)]*H) if doc['mode']=='raw' else star
    assert doc['mode'] in ('raw','relative') and tuple(map(F,doc['normalizer_row']))==norm
    assert all(a[h]<=V[j][h]<=star[h] for a,j in zip(A,js) for h in range(H))
    root=tuple((F(1),F(1)) if v==norm else (1-e,F(1)) for v,e in zip(V,eps))
    assert root==tuple(tuple(map(F,r)) for r in doc['root_box'])
    ancestor=(F(old['lower_bracket'][0]),F(old['upper_bracket'][1]));assert ancestor==tuple(map(F,doc['prior_outer_bounds']))
    env=tuple(tuple(map(F,r)) for r in doc['record_envelopes']);assert env==tuple(tuple(map(F,r)) for r in old['record_envelopes'])
    sat=F(1)
    for n,j in zip(ns,js):sat*=F(n,nj[j])**n
    R=tau/sat;assert 0<R<=1
    for (lo,hi),k,j in zip(env,ns,js):
        a=F(k,nj[j]);th=R*a**k*(1-a)**(nj[j]-k)
        assert 0<=lo<=a<=hi<=1
        assert lo==0 or lo**k*(1-lo)**(nj[j]-k)<=th
        assert hi==1 or hi**k*(1-hi)**(nj[j]-k)<=th
    cache={};counts={'nodes':0,'certificates':0,'column_inequalities':0,'supporting_planes':0,'excluded_leaves':0,'parent_witnesses':0}
    def plane(box,ref):
        ref=tuple(map(F,ref));key=(box,ref)
        if key in cache:return cache[key]
        assert len(ref)==H and min(ref)>=0 and sum(ref)==1
        scale=dot(norm,ref);assert scale>0
        pp=[dot(a,ref)/scale for a in A];assert min(pp)>0
        logs=[logint(p)[1] for p in pp]
        coeff=[]
        for h in range(H):
            T=sum((F(n,N)*(norm[h]*(l-1)+a[h]/p) for n,l,a,p in zip(ns,logs,A,pp)),F(0))
            for v,n,(lo,hi) in zip(V,nj,box):
                if lo==hi:term=-norm[h]*logint(lo)[0]
                else:term=-(hi*norm[h]-v[h])/(hi-lo)*logint(lo)[0]-(v[h]-lo*norm[h])/(hi-lo)*logint(hi)[0]
                T+=F(n,N)*term
            coeff.append(tlo*norm[h]-T)
        cache[key]=tuple(coeff);return cache[key]
    def check(box,s,cert):
        G=[];b=[]
        for v,e,(lo,hi) in zip(V,eps,box):
            G.extend([tuple(-x for x in v),tuple(x-hi*t for x,t in zip(v,norm)),tuple(lo*t-x for x,t in zip(v,norm))]);b.extend([e-1,F(0),F(0)])
        for a,j,(lo,hi) in zip(A,js,env):
            G.extend([tuple(lo*x-y for x,y in zip(V[j],a)),tuple(y-hi*x for x,y in zip(V[j],a))]);b.extend([F(0),F(0)])
        G.extend([tuple(-x for x in c),c]);b.extend([-ancestor[0],ancestor[1]])
        for ref in cert['references']:G.append(plane(box,ref));b.append(F(0));counts['supporting_planes']+=1
        yy=tuple(map(F,cert['multipliers']));z=F(cert['normalization_multiplier']);U=F(cert['upper_bound'])
        assert len(yy)==len(G) and min(yy)>=0
        obj=tuple([F(0)]*H) if cert['kind']=='empty' else tuple(s*x for x in c)
        assert all(z+sum((y*g[h] for y,g in zip(yy,G) if y),F(0))>=obj[h] for h in range(H))
        assert z+dot(yy,b)==U
        assert cert['kind'] in ('upper','empty')
        if cert['kind']=='empty':assert U<0
        counts['certificates']+=1;counts['column_inequalities']+=H
        return U
    nodes=doc['nodes'];seen=set();leaves={}
    def walk(i,box):
        assert i not in seen and 0<=i<len(nodes);seen.add(i);node=nodes[i]
        assert node['id']==i and tuple(tuple(map(F,r)) for r in node['box'])==box
        if node['status']=='split':
            j=node['axis'];lo,hi=box[j];assert lo<hi;mid=(lo+hi)/2;aa=list(box);bb=list(box);aa[j]=(lo,mid);bb[j]=(mid,hi)
            assert len(node['children'])==2;walk(node['children'][0],tuple(aa));walk(node['children'][1],tuple(bb))
        else:
            cc=node['certificates']
            if node['status']=='excluded':
                assert any(v['kind']=='empty' for v in cc.values())
                for s,cert in cc.items():check(box,int(s),cert)
                counts['excluded_leaves']+=1
            else:
                assert node['status']=='leaf' and set(cc)=={'-1','1'}
                assert all(v['kind']=='upper' for v in cc.values())
                leaves[i]={int(s):check(box,int(s),cert) for s,cert in cc.items()}
    walk(0,root);assert len(seen)==len(nodes) and sorted(leaves)==doc['leaf_ids'];counts['nodes']=len(nodes)
    inc={}
    assert set(doc['witnesses'])=={'-1','1'}
    for s,w in doc['witnesses'].items():
        m=tuple(map(F,w['mass']));assert len(m)==H and min(m)>=0 and sum(m)==1
        q=[dot(v,m) for v in V];assert all(v>=1-e for v,e in zip(q,eps));L=F(1)
        for a,n,j in zip(A,ns,js):
            x=dot(a,m);assert x>0;L*=(x/q[j])**n
        assert L>=tau and L==F(w['likelihood']);inc[int(s)]=int(s)*dot(c,m);assert inc[int(s)]==F(w['value']);counts['parent_witnesses']+=1
    upper={s:max([inc[s]]+[w[s] for w in leaves.values()]) for s in (-1,1)}
    assert tuple(map(F,doc['lower_bracket']))==(-upper[-1],-inc[-1]);assert tuple(map(F,doc['upper_bracket']))==(inc[1],upper[1])
    assert doc['precision_reached']==(max(upper[s]-inc[s] for s in (-1,1))<=F(doc['tolerance']))
    return dict(passed=True,**counts,precision_reached=doc['precision_reached'],model_rebuilt=True,original_input_unchanged=True,prior_certificate_rechecked=verify_prior,optimization_or_root_search_called=False,vertex_enumeration_called=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',default=str(BASE/'outputs'));p.add_argument('--skip-prior',action='store_true');args=p.parse_args()
    import scipy.optimize,integrated,multilimit,exact
    def forbidden(*a,**k):raise AssertionError('No optimization or vertex enumeration allowed')
    scipy.optimize.minimize=scipy.optimize.linprog=forbidden;integrated.minimize=multilimit.minimize=forbidden;exact.vertices=multilimit.vertices=forbidden
    out=Path(args.outdir);result={}
    for name in ['crossing2_large','crossing3_irregular','crossing5_consistent']:
        tic=time.perf_counter();doc=json.loads((out/f'{name}.json').read_text());old=json.loads((P/'outputs/certificates'/f'{name}.json').read_text());v=verify(doc,old,not args.skip_prior);v['seconds']=time.perf_counter()-tic;result[name]=v;print(name,v,flush=True)
    (out/'independent_verification.json').write_text(json.dumps({'all_passed':True,'cases':result,'solvers_disabled':True},sort_keys=True,indent=2)+'\n')
