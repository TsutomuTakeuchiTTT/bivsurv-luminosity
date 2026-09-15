"""Read-only Stage22 verification: no sampling, optimization or EM refitting.

Reconstructs public counts from saved individuals, validates the generator
through its conditional CDF identity, calculates truth-cell integrals using
an independent density-antiderivative formula, and recomputes training-only
predictions and rational confidence-membership inequalities.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,argparse
import numpy as np
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'vendor20'));sys.path.insert(0,str(BASE/'vendor20/vendor19'))
from integrated import build_model,quotient,em_step,observation_initial,conditional_probabilities
from exact import dot,logint,encode

def write(p,x):p.write_text(json.dumps(encode(x),sort_keys=True,indent=2)+'\n')
def integral(cell,lam):
    bounds=[]
    for I in cell:
        a=max(F(0),I.lo) if I.lo is not None else F(0)
        b=min(F(4),I.hi) if I.hi is not None else F(4)
        if b<=a:return F(0)
        bounds.append((a,b))
    (a,b),(c,d)=bounds
    G=lambda x,y:y-x-(y*y-x*x)/4
    return ((b-a)*(d-c)+lam*G(a,b)*G(c,d))/16

def run(out):
    totals={'trials':0,'parents':0,'training_predictions':0,'true_membership_checks':0,
            'raw_count_reconstructions':0,'EM_last_state_checks':0}
    max_cdf_error=0.;max_invisible_error=0.;truth_cases={};covered=0;raw_min_delta=0.
    for p in sorted((out/'observed_inputs').glob('*.json')):
        public=json.loads(p.read_text());fit=json.loads((out/'fits'/p.name).read_text())
        seed=fit['seed'];lam=[F(-3,4),F(0),F(3,4)][seed[1]]
        model=build_model(public['edges'],public['limits'],public['targets']);groups,classes=quotient(model)
        truthfine=tuple(integral(c,lam) for c in model['cells']);assert sum(truthfine)==1
        truth=tuple(sum((truthfine[h] for h in C),F(0)) for C in classes)
        law=conditional_probabilities(groups,truth)
        truth_cases[str(lam)]={'cell_integral_sum':str(sum(truthfine)),
             'selection':[str(dot(g['V'],truth)) for g in groups]}
        with np.load(out/'latent_evaluation_only'/f'{p.stem}.npz') as aa:
            z=aa['z'];which=aa['limit_index'];trflag=aa['training'];stored_delta=aa['delta'];inc0=aa['included']
        lim=np.array(public['limits'],float)[which];delta=z<=lim;inc=delta.any(axis=1)
        assert np.array_equal(stored_delta,delta) and np.array_equal(inc0,inc)
        # Algebraic inverse-CDF check, not rerunning the sampling function.
        s1,s2,s3=np.random.SeedSequence(seed).spawn(3)
        rng=np.random.default_rng(s1);u0=rng.random(len(z));w0=rng.random(len(z));u=z[:,0]/4;v=z[:,1]/4
        assert np.array_equal(u,u0)
        err=np.max(abs(v+float(lam)*(1-2*u)*v*(1-v)-w0));max_cdf_error=max(max_cdf_error,float(err));assert err<1e-12
        assert np.array_equal(which,np.random.default_rng(s2).integers(0,3,len(z)))
        assert np.array_equal(trflag,np.random.default_rng(s3).random(len(z))<.5)
        labels=np.column_stack([delta.astype(int),np.where(delta,(z>1).astype(int),-1)])
        tables=[]
        for flag in (trflag,~trflag):
            table=[]
            for j,g in enumerate(groups):
                table.append([int(np.sum((which==j)&inc&flag&np.all(labels==r['label'],axis=1))) for r in g['records']])
            tables.append(table)
        assert tables==[public['train'],public['validation']]
        assert sum(map(sum,tables[0]))+sum(map(sum,tables[1]))==int(inc.sum())
        predictor=observation_initial(groups)
        for _ in range(2):predictor,_,_=em_step(groups,tables[0],predictor)
        assert tuple(map(F,fit['prediction_mass']))==predictor
        pred=conditional_probabilities(groups,predictor)
        assert encode(pred)==fit['prediction_laws']
        lo=hi=F(0);zero=False
        for pp,tt,ns in zip(pred,law,tables[1]):
            for g,t,n in zip(pp,tt,ns):
                if not n:continue
                assert t>0
                if g==0:zero=True;break
                ll,uu=logint(g/t);lo+=n*ll;hi+=n*uu
        cutoff=logint(F(1)/F(public['alpha']))
        assert encode(cutoff)==fit['log_rejection_threshold_interval']
        if zero:member=True
        else:
            assert encode((lo,hi))==fit['truth_membership_logE_interval']
            if hi<=cutoff[0]:member=True
            elif lo>cutoff[1]:member=False
            else:member=None
        assert member is fit['evaluation']['true_in_confidence_set']
        assert all(dot(g['V'],truth)>=1-F(e) for g,e in zip(groups,public['epsilon']))
        covered+=member is True
        mass=np.array(fit['EM_mass']);hist=fit['history'];assert mass.min()>=0 and abs(mass.sum()-1)<1e-12
        inv=[h for h in range(len(mass)) if all(not g['V'][h] for g in groups)]
        inv_error=abs(float(mass[inv].sum())-len(inv)/len(mass));max_invisible_error=max(max_invisible_error,inv_error);assert inv_error<1e-12
        last=hist[-1];assert np.array_equal(mass,np.array(last['mass']))
        raw_min_delta=min(raw_min_delta,*np.diff([h['ell'] for h in hist]))
        totals['trials']+=1;totals['parents']+=len(z)
        for k in ('training_predictions','true_membership_checks','raw_count_reconstructions','EM_last_state_checks'):totals[k]+=1
    # Reverify both generations of endpoint proof, disabling numeric optimizers.
    import scipy.optimize,integrated,multilimit,exact
    def forbidden(*a,**k):raise AssertionError('Numerical optimization is forbidden during verification')
    scipy.optimize.minimize=scipy.optimize.linprog=forbidden
    integrated.minimize=multilimit.minimize=forbidden;exact.vertices=multilimit.vertices=forbidden
    sys.path.insert(0,str(BASE/'vendor21'))
    from verify_stage21 import verify
    proofs={}
    for p in sorted((out/'refined_projection_certificates').glob('*.json')):
        new=json.loads(p.read_text());old=json.loads((out/'projection_certificates'/p.name).read_text())
        proofs[p.stem]=verify(new,old,True)
    res={'all_scoped_checks_passed':True,'counts':totals,'covered_parent_sets':covered,
         'inverse_conditional_cdf_max_abs_error':max_cdf_error,
         'invisible_initial_mass_max_abs_error':max_invisible_error,
         'raw_min_ell_increment':float(raw_min_delta),'truth_cell_integrals':truth_cases,
         'projection_certificates':proofs,'optimization_called':False,
         'point_EM_refitted':False,'sampler_called':False,
         'scope':'Saved input and proof verification, not an independent proof of sampling coverage.'}
    write(out/'independent_verification.json',res);return res
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs');a=p.parse_args()
    print(json.dumps(run(a.outdir),indent=2))
