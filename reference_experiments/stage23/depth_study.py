"""Stage23: paired re-observation of the saved Stage22 parent samples.

No new latent draws. Fixed common reference measure across six designs;
original D-completion EM. Population fibers are separate from sample inference.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import sys, json, csv, hashlib, time, platform
from fractions import Fraction as F
import numpy as np
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'vendor20'));sys.path.insert(0,str(BASE/'vendor20/vendor19'))
from integrated import build_model, quotient, aggregate_masses, em_step, numeric_em, conditional_probabilities
from coarsened_model import FGMLaw, measure_masses, record, Interval
from exact import dot, logint, encode

TARGETS=[['1','1'],['3/2','1'],['2','2']]
LAMBDAS=(F(-3,4),F(0),F(3,4))
ALPHA=F(1,20)

def save(p:Path,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(encode(x),ensure_ascii=False,sort_keys=True,indent=2)+'\n')

def save_csv(p:Path,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

class ReferenceLaw:
    """A fixed continuous probability law, NOT fitted and NOT the true law.
    Independent margins put 1/3 in each original interval (-inf,1],(1,2],(2,inf).
    All subsequent partitions are pushforwards of this one measure.
    """
    @staticmethod
    def cdf(x):
        if x is None:return F(1)
        x=F(x)
        if x<=1:return 1/(3*(2-x))
        if x<=2:return x/3
        return 1-1/(3*(x-1))
    def rectangle(self,cell):
        value=F(1)
        for I in cell:
            value*=self.cdf(I.hi)-(self.cdf(I.lo) if I.lo is not None else 0)
        return value


def setup(resolution,depth):
    edges=['1'] if resolution=='coarse' else ['1','3/2','2','3','4']
    spec={'resolution':resolution,'depth':depth,'edges':[edges,edges],
          'limits':[[2,1],[1,2],[depth,depth]],'targets':TARGETS,
          'epsilon':['1/2']*3,'alpha':'1/20',
          'initial_measure':'fixed_reference_R; not uniform after refining'}
    model=build_model(spec['edges'],spec['limits'],spec['targets'])
    obs,classes=quotient(model)
    allg,jclasses=quotient(model,model['targets_coefficients'])
    ref_fine=measure_masses(model,ReferenceLaw())
    init=aggregate_masses(ref_fine,classes)
    assert sum(init)==1 and min(init)>0
    # conditional reference target expectation inside an observation class.
    coeff=tuple(tuple(sum((ref_fine[h]*c[h] for h in C),F(0))/init[j]
                       for j,C in enumerate(classes)) for c in model['targets_coefficients'])
    inv=tuple(h for h in range(len(init)) if all(not g['V'][h] for g in obs))
    assert len(inv)==1
    assert init[inv[0]]==F(1,9*(depth-1)**2)
    v0,v1=obs[0]['V'],obs[1]['V']
    assert any(a>b for a,b in zip(v0,v1)) and any(b>a for a,b in zip(v0,v1))
    assert all((1-obs[2]['V'][h]) <= (1-v0[h]) and
               (1-obs[2]['V'][h]) <= (1-v1[h]) for h in range(len(init)))
    return spec,model,obs,classes,init,coeff,inv


def public_counts(z,which,training,model,groups):
    lim=np.asarray([[float(x) for x in c] for c in model['limits']])[which]
    delta=z<=lim; included=delta.any(axis=1)
    labels=np.column_stack([delta.astype(int)]+[
        np.where(delta[:,b],np.searchsorted(np.array(model['edges'][b],float),z[:,b],side='left'),-1)
        for b in (0,1)])
    tables=[]
    for flag in (training,~training):
        tables.append([[int(np.sum((which==j)&included&flag&np.all(labels==r['label'],axis=1)))
                         for r in g['records']] for j,g in enumerate(groups)])
    assert sum(map(sum,tables[0]))+sum(map(sum,tables[1]))==int(included.sum())
    return *tables,delta,included,labels


def evaluate_membership(pred,truth,counts):
    lo=hi=F(0)
    for pp,tt,ns in zip(pred,truth,counts):
        for p,t,n in zip(pp,tt,ns):
            if not n:continue
            if t<=0:raise AssertionError('impossible true record')
            if p==0:return True,None,'zero_prediction_noninformative'
            l,u=logint(p/t);lo+=n*l;hi+=n*u
    cl,cu=logint(1/ALPHA)
    if hi<=cl:return True,(lo,hi),'included'
    if lo>cu:return False,(lo,hi),'excluded'
    return None,(lo,hi),'arithmetic_undecided'


def run(out):
    start=time.perf_counter();out.mkdir(parents=True,exist_ok=True)
    paths=sorted((BASE/'inputs/latent_evaluation_only').glob('*.npz'))
    if len(paths)!=120:raise ValueError('Exactly the 120 saved Stage22 samples are required')
    latent={}
    for p in paths:
        with np.load(p) as d:latent[p.stem]={k:d[k].copy() for k in ('z','limit_index','training')}
    rows=[];checks={'baseline_public_exact':0,'baseline_prediction_exact':0,
       'baseline_fit_checks':0,'exact_record_integrals':0,'geometry_comparisons':0,
       'paired_detection_monotonicity':0,'fixed_prior_pushforward':0}
    originals={p.stem:json.loads(p.read_text()) for p in (BASE/'inputs/stage22_fits').glob('*.json')}
    actual_difference=0.;all_states={};designs=[]
    for res in ('coarse','aligned'):
        for depth in (2,3,4):
            name=f'{res}_d{depth}'
            spec,model,obs,classes,init,coeff,inv=setup(res,depth)
            specs={'specification':spec,'full_cells':len(model['cells']),
                   'observation_classes':len(classes),'target_classes':len(quotient(model,model['targets_coefficients'])[1]),
                   'records':[len(g['records']) for g in obs], 'initial_mass':init,
                   'initial_invisible_mass':init[inv[0]],'reference_target_coefficients':coeff}
            save(out/'designs'/f'{name}.json',specs);designs.append(specs)
            truths=[]
            for lam in LAMBDAS:
                law=FGMLaw(lam);mass=aggregate_masses(measure_masses(model,law),classes)
                probs=conditional_probabilities(obs,mass)
                for j,g in enumerate(obs):
                    direct,q=law.direct_law(g['limit'],model['bins'])
                    assert q==dot(g['V'],mass) and q>=F(1,2)
                    for r,p in zip(g['records'],probs[j]):
                        assert direct[r['label']]/q==p;checks['exact_record_integrals']+=1
                truths.append((mass,probs,[law.survival(t) for t in model['targets']]))
            for case,arr in latent.items():
                li=int(case.split('_')[0][-1]);nparent=int(case.split('_')[1][1:])
                lam=LAMBDAS[li];truth,truthlaw,trueS=truths[li]
                z=arr['z'];which=arr['limit_index'];trflag=arr['training']
                tr,va,delta,inc,labels=public_counts(z,which,trflag,model,obs)
                oldinc=all_states.get((res,depth-1,case))
                if oldinc is not None:
                    assert np.all(~oldinc|inc);checks['paired_detection_monotonicity']+=1
                all_states[(res,depth,case)]=inc
                # Compare a bounded set of exact-coordinate records with vectorization.
                for k in range(8):
                    zz=tuple(F.from_float(float(x)) for x in z[k])
                    lab=record(zz,model['limits'][which[k]],model['bins'])
                    assert lab==(tuple(map(int,labels[k])) if inc[k] else None)
                    checks['geometry_comparisons']+=1
                public={**spec,'case':case,'train':tr,'validation':va}
                # Inference functions take only public counts and predeclared initial measure.
                save(out/'observed_inputs'/name/f'{case}.json',public)
                public=json.loads((out/'observed_inputs'/name/f'{case}.json').read_text())
                tr,va=public['train'],public['validation']
                if min(map(sum,tr))==0 or min(map(sum,va))==0:
                    raise RuntimeError('Empty layer must have an explicitly reported inference policy')
                predictor=init
                for _ in range(2):predictor,_,_=em_step(obs,tr,predictor)
                pred=conditional_probabilities(obs,predictor)
                member,loge,decision=evaluate_membership(pred,truthlaw,va)
                counts=[[a+b for a,b in zip(x,y)] for x,y in zip(tr,va)]
                m,hist=numeric_em(obs,counts,init,max_iter=4000,tol=1e-8)
                shat=np.array(coeff,float)@m
                nt=sum(map(sum,counts));tv=0.
                for g,nn,p0 in zip(obs,counts,truthlaw):
                    q=np.dot(g['V'],m)
                    pp=np.array([np.dot(r['A'],m)/q for r in g['records']])
                    tv+=sum(nn)/nt*.5*np.abs(pp-np.array(p0,float)).sum()
                selected=np.array([np.mean(np.all(z[inc]>np.array(t,float),axis=1)) for t in model['targets']])
                evalrow={'design':name,'resolution':res,'depth':depth,'case':case,'lambda':str(lam),
                    'parent_n':nparent,'observed_n':int(inc.sum()),'train_n':sum(map(sum,tr)),
                    'validation_n':sum(map(sum,va)),'region_D':int((~inc).sum()),
                    'parent_membership':member,'membership_status':decision,'iterations':len(hist)-1,
                    'score':hist[-1]['positive_score'],'stopped':bool(hist[-1]['positive_score']<=1e-8),
                    'tail_feasible':all(np.dot(g['V'],m)>=.5-1e-12 for g in obs),
                    'invisible_initial':float(init[inv[0]]),'invisible_final':float(m[inv[0]]),
                    'true_invisible':float(truth[inv[0]]),'observed_law_TV':float(tv),
                    'reference_RMSE':float(np.sqrt(np.mean((shat-np.array(trueS,float))**2))),
                    'selected_oracle_RMSE':float(np.sqrt(np.mean((selected-np.array(trueS,float))**2))),
                    'min_raw_ell_gain':float(min(np.diff([h['ell'] for h in hist]),default=0.))}
                for k in range(3):evalrow[f'Shat_{k}']=float(shat[k]);evalrow[f'Strue_{k}']=float(trueS[k])
                if res=='coarse' and depth==2:
                    original=json.loads((BASE/'inputs/stage22_observed_inputs'/f'{case}.json').read_text())
                    assert tr==original['train'] and va==original['validation'];checks['baseline_public_exact']+=1
                    prev=originals[case]
                    assert encode(predictor)==prev['prediction_mass'];checks['baseline_prediction_exact']+=1
                    err=max(float(np.max(np.abs(m-np.array(prev['EM_mass'])))),
                            *(abs(evalrow[f'Shat_{k}']-prev['evaluation'][f'Shat_{k}']) for k in range(3)))
                    actual_difference=max(actual_difference,err)
                    assert err<1e-12;checks['baseline_fit_checks']+=1
                assert abs(m[inv[0]]-float(init[inv[0]]))<1e-12
                rows.append(evalrow)
                save(out/'fits'/name/f'{case}.json',{'prediction':predictor,'laws':pred,'membership_logE':loge,
                     'initial_mass':init,'final_mass':m.tolist(),'evaluation':evalrow,
                     'history_ell':[h['ell'] for h in hist]})
            print('completed',name,len(obs[0]['V']),'masses',flush=True)
    # The same reference measure gives identical coarse class totals under actual refinements.
    for depth in (2,3,4):
        _,cm,cg,cc,ci,_,_=setup('coarse',depth)
        _,fm,fg,fc,fi,_,_=setup('aligned',depth)
        # Map each fine class using its exact representative into the coarse observation signature.
        for j,C in enumerate(cc):
            labels=[record(cm['representatives'][C[0]],g['limit'],cm['bins']) for g in cg]
            inds=[k for k,Fc in enumerate(fc) if
                  [record(fm['representatives'][Fc[0]],g['limit'],cm['bins']) for g in cg]==labels]
            assert sum((fi[k] for k in inds),F(0))==ci[j]
            checks['fixed_prior_pushforward']+=1
    summary=[]
    for res in ('coarse','aligned'):
        for depth in (2,3,4):
            for lam in LAMBDAS:
                for n in (400,1600):
                    rr=[r for r in rows if r['resolution']==res and r['depth']==depth and r['lambda']==str(lam) and r['parent_n']==n]
                    s={k:rr[0][k] for k in ('resolution','depth','lambda','parent_n')}
                    s.update(replicates=len(rr),membership_accepted=sum(r['parent_membership'] is True for r in rr),
                             membership_undecided=sum(r['parent_membership'] is None for r in rr),
                             em_stopped=sum(r['stopped'] for r in rr),outside_tail=sum(not r['tail_feasible'] for r in rr))
                    for k in ('observed_n','observed_law_TV','reference_RMSE','invisible_initial','true_invisible','iterations','Shat_0','Shat_1','Shat_2'):
                        s['mean_'+k]=float(np.mean([r[k] for r in rr]))
                    summary.append(s)
    save_csv(out/'trials.csv',rows);save_csv(out/'summary.csv',summary)
    result={'stage':'23','paired_distinct_parent_samples':120,'unique_parent_objects':120000,
            'design_evaluations':len(rows),'survival_reference_values':3*len(rows),'new_random_draws':0,
            'recording_designs':6,'all_sharp_sample_projections_performed':False,
            'checks':checks,'baseline_fit_max_abs_difference':actual_difference,
            'all_invisible_mass_invariants_checked':True,
            'membership_accepted':sum(r['parent_membership'] is True for r in rows),
            'membership_undecided':sum(r['parent_membership'] is None for r in rows),
            'all_em_stopped':all(r['stopped'] for r in rows),'outside_tail':sum(not r['tail_feasible'] for r in rows),
            'max_iterations':max(r['iterations'] for r in rows),'minimum_iterations':min(r['iterations'] for r in rows),
            'software':{'python':platform.python_version(),'numpy':np.__version__},'seconds':time.perf_counter()-start}
    save(out/'sampling_checks.json',result)
    return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs')
    args=p.parse_args();print(json.dumps(run(args.outdir),indent=2))
