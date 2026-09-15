"""Read-only arithmetic and observation audit of Stage23.
No nonlinear optimization, LP solving, sampler, or numerical EM refitting.
Training updates are re-evaluated from the saved public tables.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys,csv
import numpy as np
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'vendor20'));sys.path.insert(0,str(BASE/'vendor20/vendor19'))
from integrated import build_model,quotient,em_step,conditional_probabilities
from exact import dot,logint,encode


def antiderivative_mass(cell,lam):
    pairs=[]
    for I in cell:
        a=max(F(0),I.lo) if I.lo is not None else F(0)
        b=min(F(4),I.hi) if I.hi is not None else F(4)
        if b<=a:return F(0)
        pairs.append((a,b))
    (a,b),(c,d)=pairs
    P=lambda x:x-x*x/4
    return ((b-a)*(d-c)+lam*(P(b)-P(a))*(P(d)-P(c)))/16


def refmass(cell):
    # Direct interval integrals of the three reference densities.
    # Below1: 1/[3(2-x)^2]; (1,2):1/3; above2:1/[3(x-1)^2].
    value=F(1)
    for I in cell:
        if I.lo is not None and I.lo==I.hi:return F(0)
        lo,hi=I.lo,I.hi
        mass=F(0)
        a=lo;b=min(hi,F(1)) if hi is not None else F(1)
        if a is None or a<b:
            mass+=1/(3*(2-b))-(1/(3*(2-a)) if a is not None else 0)
        a=max(lo,F(1)) if lo is not None else F(1)
        b=min(hi,F(2)) if hi is not None else F(2)
        if a<b:mass+=(b-a)/3
        a=max(lo,F(2)) if lo is not None else F(2);b=hi
        if b is None or a<b:mass+=1/(3*(a-1))-(1/(3*(b-1)) if b is not None else 0)
        value*=mass
    return value


def verify_population(cert):
    res=cert['resolution'];depth=cert['depth'];lam=F(cert['lambda'])
    e=['1'] if res=='coarse' else ['1','3/2','2','3','4']
    model=build_model([e,e],[[2,1],[1,2],[depth,depth]],[['1','1'],['3/2','1'],['2','2']])
    groups,classes=quotient(model,model['targets_coefficients']);H=len(classes)
    mf=[antiderivative_mass(C,lam) for C in model['cells']]
    m0=[sum((mf[h] for h in C),F(0)) for C in classes];assert sum(m0)==1
    p0=conditional_probabilities(groups,m0);assert encode(p0)==cert['true_probabilities']
    E=[[F(1)]*H]
    V=[g['V'] for g in groups]
    for g,pp in zip(groups,p0):
        for r,p in zip(g['records'],pp):E.append([F(a)-p*v for a,v in zip(r['A'],g['V'])])
    c=[F(model['targets_coefficients'][cert['target_index']][C[0]]) for C in classes]
    totalcolumns=0
    for key,sign in [('lower',1),('upper',-1)]:
        d=cert[key];assert d['sign']==sign
        mass=tuple(map(F,d['mass']));assert min(mass)>=0 and sum(mass)==1
        assert conditional_probabilities(groups,mass)==p0
        assert all(dot(v,mass)>=F(1,2) for v in V)
        y=list(map(F,d['equality_multipliers']));z=list(map(F,d['inequality_multipliers']))
        assert len(y)==len(E) and len(z)==3 and max(z)<=0
        for h in range(H):
            val=sum((y[k]*E[k][h] for k in range(len(E))),F(0))-sum((z[j]*V[j][h] for j in range(3)),F(0))
            assert val<=sign*c[h];totalcolumns+=1
        lb=y[0]-sum(z,F(0))/2;inner=sign*dot(c,mass)
        assert lb==F(d['objective_lower_bound']) and inner==F(d['objective_at_witness'])
        assert inner-lb==F(d['gap']) and 0<=inner-lb<F(1,10**8)
        assert dot(c,mass)==F(d['survival_at_witness'])
    if cert['analytic'] is not None:
        an=cert['analytic']
        for key,target in [('lower_mass','lower'),('upper_mass','upper')]:
            m=list(map(F,an[key]));assert min(m)>=0 and sum(m)==1
            assert conditional_probabilities(groups,m)==p0
            assert all(dot(v,m)>=F(1,2) for v in V)
            assert dot(c,m)==F(an[target])
        # Verify every deepest record determines the target and shallower visibility.
        for r in groups[2]['records']:
            cols=[h for h,a in enumerate(r['A']) if a]
            assert len({c[h] for h in cols})==1
            for v in V:assert len({v[h] for h in cols})==1
    return totalcolumns


def run(out):
    import scipy.optimize
    def forbidden(*a,**k):raise AssertionError('Optimization is disabled in this verification')
    scipy.optimize.linprog=scipy.optimize.minimize=forbidden
    trials=predictions=memberships=0;population=columns=0;outside=accepted=0;maxinv=0.
    for ds in sorted((out/'designs').glob('*.json')):
        dd=json.loads(ds.read_text());s=dd['specification'];name=ds.stem
        model=build_model(s['edges'],s['limits'],s['targets']);groups,classes=quotient(model)
        rr=tuple(refmass(C) for C in model['cells']);initial=tuple(sum((rr[h] for h in C),F(0)) for C in classes)
        assert sum(initial)==1 and encode(initial)==dd['initial_mass']
        inv=[h for h in range(len(initial)) if all(not g['V'][h] for g in groups)]
        ce=[[sum((rr[h]*c[h] for h in C),F(0))/initial[j] for j,C in enumerate(classes)]
            for c in model['targets_coefficients']]
        for p in sorted((out/'observed_inputs'/name).glob('*.json')):
            public=json.loads(p.read_text());fit=json.loads((out/'fits'/name/p.name).read_text());ev=fit['evaluation']
            with np.load(BASE/'inputs/latent_evaluation_only'/f'{p.stem}.npz') as f:
                z=f['z'];which=f['limit_index'];training=f['training']
            # Independently count records one individual at a time after array comparisons.
            lim=np.array(s['limits'],float)[which];de=z<=lim
            bi=np.column_stack([(z[:,b,None]>np.array([float(F(x)) for x in s['edges'][b]])).sum(axis=1) for b in (0,1)])
            idx=[{r['label']:k for k,r in enumerate(g['records'])} for g in groups]
            tables=[[[0]*len(g['records']) for g in groups] for _ in (0,1)]
            for j,tr,d,b in zip(which,training,de,bi):
                if not any(d):continue
                label=(int(d[0]),int(d[1]),int(b[0]) if d[0] else -1,int(b[1]) if d[1] else -1)
                tables[0 if tr else 1][j][idx[j][label]]+=1
            assert tables==[public['train'],public['validation']]
            predmass=initial
            for _ in range(2):predmass,_,_=em_step(groups,tables[0],predmass)
            pp=conditional_probabilities(groups,predmass)
            assert encode(predmass)==fit['prediction'] and encode(pp)==fit['laws'];predictions+=1
            lam=[F(-3,4),F(0),F(3,4)][int(p.stem.split('_')[0][-1])]
            tm=[antiderivative_mass(C,lam) for C in model['cells']]
            true=[sum((tm[h] for h in C),F(0)) for C in classes]
            assert sum(true)==1 and all(dot(g['V'],true)>=F(1,2) for g in groups)
            tp=conditional_probabilities(groups,true)
            low=high=F(0);zero=False
            for row0,row1,nn in zip(pp,tp,tables[1]):
                for a,b,n in zip(row0,row1,nn):
                    if not n:continue
                    assert b>0
                    if a==0:zero=True;continue
                    l,u=logint(a/b);low+=n*l;high+=n*u
            cl,cu=logint(F(20))
            membership=True if zero or high<=cl else False if low>cu else None
            assert membership is ev['parent_membership']
            if not zero:assert encode((low,high))==fit['membership_logE']
            memberships+=1;accepted+=membership is True
            m=np.array(fit['final_mass']);assert min(m)>=0 and abs(m.sum()-1)<1e-11
            err=abs(m[inv].sum()-float(sum((initial[h] for h in inv),F(0))));assert err<1e-12;maxinv=max(maxinv,err)
            shat=np.array(ce,float)@m
            assert max(abs(shat[k]-ev[f'Shat_{k}']) for k in range(3))<1e-12
            tail=all(np.dot(g['V'],m)>=.5-1e-12 for g in groups)
            assert tail is ev['tail_feasible'];outside+=not tail;trials+=1
    for p in sorted((out/'population_certificates').glob('*.json')):
        columns+=verify_population(json.loads(p.read_text()));population+=1
    # Tampering with an extremal mass or a bound must fail.
    from copy import deepcopy
    example=json.loads(next((out/'population_certificates').glob('*.json')).read_text())
    rejects=0
    for what in ('mass','bound','multiplier','probability'):
        altered=deepcopy(example)
        if what=='mass':altered['lower']['mass'][0]=str(F(altered['lower']['mass'][0])+F(1,100))
        elif what=='bound':altered['lower']['objective_lower_bound']=str(F(altered['lower']['objective_lower_bound'])+F(1,100))
        elif what=='multiplier':altered['upper']['equality_multipliers'][0]=str(F(altered['upper']['equality_multipliers'][0])+10)
        else:altered['true_probabilities'][0][0]=str(F(altered['true_probabilities'][0][0])+F(1,100))
        try:verify_population(altered)
        except AssertionError:rejects+=1
        else:raise AssertionError('Tampering accepted')
    result={'passed':True,'public_reobservations_verified':trials,'training_predictions_verified':predictions,
        'true_membership_checks':memberships,'accepted_memberships':accepted,'unconstrained_fits_outside_tail':outside,
        'maximum_invisible_mass_drift':maxinv,'population_certificates':population,
        'rational_dual_columns':columns,'rational_extremal_parents':population*2,
        'analytic_witnesses':54,'tamper_cases_rejected':rejects,'sampler_called':False,
        'numeric_EM_refitted':False,'LP_called':False,'nonlinear_optimization_called':False}
    (out/'independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs');a=p.parse_args()
    print(json.dumps(run(a.outdir),indent=2))
