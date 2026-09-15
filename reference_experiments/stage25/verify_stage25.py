"""Independent global-gap verifier. No optimizer, EM updates, or root search.
Rebuilds the full parent observation law from the public input, then evaluates
homogeneous affine bounds at simplex basis points and all-column dual tests.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,copy
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'reference_stage24'))
from verify_stage24 import data_from_public
from exact import logint


def dot(a,b):return sum((F(x)*F(y) for x,y in zip(a,b)),F(0))
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(y) for y in x]
    return x


def verify_case(path):
    cert=json.loads(Path(path).read_text());name=cert['case']
    snapshot=BASE/'inputs'/(name+'.json');raw=json.loads(snapshot.read_text())
    assert hashlib.sha256(snapshot.read_bytes()).hexdigest()==cert['source_input_sha256']
    inp=BASE/'reference_stage24/inputs/observed'/(name+'.json')
    source=BASE/'reference_stage24/outputs/certificates'/(name+'.json')
    assert hashlib.sha256(inp.read_bytes()).hexdigest()==raw['input_sha256']
    assert hashlib.sha256(source.read_bytes()).hexdigest()==raw['stage24_certificate_sha256']
    assert json.loads(inp.read_text())==raw['public_input']
    source=json.loads(source.read_text());assert source['independent_candidate']==raw['candidate']
    full=data_from_public(raw['public_input'])
    assert serial({k:full[k] for k in raw['compiled']})==raw['compiled']
    m=list(map(F,raw['candidate']['mass']));assert len(m)==len(full['initial']) and min(m)>=0 and sum(m)==1
    assert all(dot(v,m)>=F(1,2) for v in full['V'])
    eta=[max(v[h] for v in full['V']) for h in range(len(m))]
    assert eta==list(full['V'][2])
    visible=[h for h,e in enumerate(eta) if e]
    assert visible==cert['visible_indices']
    assert all(all(a[h]==0 for a in full['A']) for h,e in enumerate(eta) if not e)
    r=dot(eta,m);assert r>0
    mn=[m[h]/r for h in visible];assert mn==list(map(F,cert['normalized_candidate']))
    A=[[a[h] for h in visible] for a in full['A']]
    V=[[v[h] for h in visible] for v in full['V']]
    VR=[[v[h] for h in visible] for v in full['VR']]
    assert all(v==1 for v in V[2])
    assert min(dot(v,mn) for v in V)>=F(1,2)
    weights=[F(n,full['N']) for n in full['counts']]
    layer=[V.index(v) for v in VR]
    lam=[sum((w for w,jj in zip(weights,layer) if jj==j),F(0)) for j in range(3)]
    # Re-evaluate the lower bound at the ORIGINAL sixteen-mass candidate.
    lower=upper=F(0)
    for aa,vv,w,a,v in zip(full['A'],full['VR'],weights,A,VR):
        ratio=dot(aa,m)/dot(vv,m);assert ratio==dot(a,mn)/dot(v,mn)>0
        l,u=logint(ratio);lower+=w*l;upper+=w*u
    assert [lower,upper]==list(map(F,cert['mean_lower_interval']))
    nodes=cert['tree'];assert all(n['id']==k for k,n in enumerate(nodes))
    seen=set();leaves=[];columns=0
    def visit(k,box):
        nonlocal columns
        assert k not in seen and 0<=k<len(nodes);seen.add(k)
        node=nodes[k];bb=tuple(tuple(map(F,t)) for t in node['box']);assert bb==box
        assert len(box)==2 and all(F(1,2)<=l<=u<=1 for l,u in box)
        proof=node['proof'];v=list(map(F,proof['reference']))
        assert len(v)==len(visible) and min(v)>=0 and sum(v)==1
        ai=[dot(a,v) for a in A];assert min(ai)>0
        coeff=[]
        # Independent basis-value formula; no call to tangent_coeff/polytope.
        for h in range(len(v)):
            value=sum((w*(logint(x)[1]+F(a[h])/x-1) for a,x,w in zip(A,ai,weights)),F(0))
            for j,(l,u) in enumerate(box):
                aup=-logint(l)[0];bup=-logint(u)[0]
                if l==u:value+=lam[j]*aup
                else:value+=lam[j]*((u-V[j][h])*aup+(V[j][h]-l)*bup)/(u-l)
            coeff.append(value)
        y=list(map(F,proof['multipliers']));z=F(proof['simplex_multiplier'])
        assert len(y)==4 and min(y)>=0
        for h in range(len(v)):
            val=z+sum((V[j][h]*(y[2*j]-y[2*j+1]) for j in range(2)),F(0))
            assert val>=coeff[h];columns+=1
        bound=z+sum((y[2*j]*box[j][1]-y[2*j+1]*box[j][0] for j in range(2)),F(0))
        assert bound==F(proof['upper'])
        if 'children' in node:
            j=node['split_axis'];assert j in (0,1) and len(node['children'])==2
            l,u=box[j];assert l<u;mid=(l+u)/2
            b1=list(box);b2=list(box);b1[j]=(l,mid);b2[j]=(mid,u)
            visit(node['children'][0],tuple(b1));visit(node['children'][1],tuple(b2))
        else:leaves.append(bound)
    visit(0,((F(1,2),F(1)),(F(1,2),F(1))))
    assert len(seen)==len(nodes)
    global_upper=max(leaves);gap=global_upper-lower
    assert global_upper==F(cert['mean_upper']) and gap==F(cert['gap_upper']) and gap>=0
    assert (gap<=F(cert['tolerance']))==cert['precision_achieved']
    return {'case':name,'nodes':len(nodes),'leaves':len(leaves),'column_inequalities':columns,
       'gap_upper':str(gap),'gap_upper_float':float(gap),'precision_achieved':cert['precision_achieved'],
       'original_parent_mass_unchanged':True,'likelihood_ratios_preserved':len(A)}


def run(certdir,out):
    import scipy.optimize
    def forbidden(*a,**k):raise AssertionError('Optimizer prohibited during independent verification')
    for name in ('minimize','linprog','root','brentq','differential_evolution'):setattr(scipy.optimize,name,forbidden)
    paths=sorted(Path(certdir).glob('coarse_*.json'));results=[verify_case(p) for p in paths]
    # Harmful modifications, including changes that could fake full coverage.
    first=json.loads(paths[0].read_text());tmp=BASE/'outputs/tamper_test.json';rejected=[]
    for typ in ('global_upper','negative_dual','root_box','child_coverage','candidate','missing_column','tangent'):
        t=copy.deepcopy(first)
        if typ=='global_upper':t['mean_upper']=str(F(t['mean_lower_interval'][0])-1)
        if typ=='negative_dual':t['tree'][0]['proof']['multipliers'][0]='-1'
        if typ=='root_box':t['tree'][0]['box'][0][0]='3/4'
        if typ=='child_coverage':t['tree'][1]['box'][0][0]='7/8'
        if typ=='candidate':t['normalized_candidate'][0]=str(F(t['normalized_candidate'][0])+F(1,10))
        if typ=='missing_column':t['visible_indices']=t['visible_indices'][:-1]
        if typ=='tangent':t['tree'][0]['proof']['simplex_multiplier']=str(F(t['tree'][0]['proof']['simplex_multiplier'])-1)
        tmp.write_text(json.dumps(t))
        try:verify_case(tmp)
        except (AssertionError,ValueError,ZeroDivisionError,IndexError):rejected.append(typ)
        else:raise AssertionError('Harmful tamper accepted '+typ)
    tmp.unlink()
    result={'passed':True,'cases':results,'case_count':len(results),'nodes':sum(r['nodes'] for r in results),
       'leaves':sum(r['leaves'] for r in results),'column_inequalities':sum(r['column_inequalities'] for r in results),
       'maximum_gap_upper':max((r['gap_upper_float'] for r in results),default=0),
       'harmful_tampers_rejected':rejected,'optimizer_called':False,'GEM_refitted':False,'new_samples_generated':False}
    Path(out).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print('VERIFIED',len(results),result['nodes'],result['column_inequalities'],result['maximum_gap_upper'])
    return result
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--certdir',default=str(BASE/'outputs/main_certificates'));p.add_argument('--output',default=str(BASE/'outputs/independent_verification.json'));a=p.parse_args();run(a.certdir,a.output)
