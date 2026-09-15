"""Stage 22A: continuous-parent, nonnested-limit sampling pilot.

Actual parent pairs are simulated, then independently assigned limits and a
training/validation flag BEFORE union selection. Only public coarsened counts
enter the unchanged learner/EM. Truth and discarded D objects are evaluation-only.
This pilot does not claim global optimality of numerical EM or asymptotic coverage.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import sys,json,hashlib,argparse,csv,time,platform
from fractions import Fraction as F
import numpy as np
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'vendor20'));sys.path.insert(0,str(BASE/'vendor20/vendor19'))
from integrated import (build_model,quotient,aggregate_masses,em_step,
                        observation_initial,numeric_em,conditional_probabilities)
from coarsened_model import FGMLaw,measure_masses,record
from exact import dot,logint,encode

EDGES=[['1'],['1']]
LIMITS=[[2,1],[1,2],[2,2]]
TARGETS=[['1','1'],['3/2','1'],['2','2']]
EPS=[F(1,2)]*3
ALPHA=F(1,20)

def write_json(p:Path,x:object)->None:
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(encode(x),ensure_ascii=False,sort_keys=True,indent=2)+'\n')

def sample_parent(lam:F,n:int,rng:np.random.Generator)->np.ndarray:
    """Inverse conditional FGM sampling; model has density on (0,4)^2."""
    u=rng.random(n);w=rng.random(n);a=float(lam)*(1-2*u)
    v=2*w/(1+a+np.sqrt((1+a)**2-4*a*w))
    return 4*np.column_stack([u,v])

def public_counts(z:np.ndarray,which:np.ndarray,training:np.ndarray,groups):
    """Public records only. No true masses, support box, or parent n in output."""
    lim=np.asarray(LIMITS,float)[which];delta=z<=lim
    labels=np.column_stack([delta.astype(int),np.where(delta,(z>1).astype(int),-1)])
    included=np.any(delta,axis=1)
    tr=[];va=[]
    for j,g in enumerate(groups):
        tr.append([int(np.sum((which==j)&included&training&np.all(labels==r['label'],axis=1))) for r in g['records']])
        va.append([int(np.sum((which==j)&included&~training&np.all(labels==r['label'],axis=1))) for r in g['records']])
    return tr,va,delta,included,labels

def log_e_interval(pred,truth,counts):
    lo=hi=F(0)
    for g,p,nn in zip(pred,truth,counts):
        for v,t,n in zip(g,p,nn):
            if not n:continue
            if t<=0:raise AssertionError('Impossible outcome under generator')
            if v==0:return None
            l,u=logint(v/t);lo+=n*l;hi+=n*u
    return lo,hi

def run(outdir:Path,reps:int=20)->dict:
    start=time.perf_counter();outdir.mkdir(parents=True,exist_ok=True)
    model=build_model(EDGES,LIMITS,TARGETS)
    obs,classes=quotient(model)
    joint,jclasses=quotient(model,model['targets_coefficients'])
    parents={h:j for j,C in enumerate(classes) for h in C}
    jparent=[parents[C[0]] for C in jclasses]
    multiplicity=[jparent.count(j) for j in range(len(classes))]
    init=observation_initial(obs)
    invisible=[j for j in range(len(classes)) if all(not g['V'][j] for g in obs)]
    assert len(invisible)==1
    # Verify nonnested inclusion sets with actual positive incidence columns.
    v0,v1=obs[0]['V'],obs[1]['V']
    assert any(a>b for a,b in zip(v0,v1)) and any(b>a for a,b in zip(v0,v1))
    rows=[];checks={'record_comparisons':0,'hidden_value_invariance':0,'exact_probability_checks':0,'split_disjointness':0}
    models=[]
    for li,lam in enumerate([F(-3,4),F(0),F(3,4)]):
        law=FGMLaw(lam);truefine=measure_masses(model,law)
        truth=aggregate_masses(truefine,classes);trueprob=conditional_probabilities(obs,truth)
        trueS=[law.survival(t) for t in model['targets']]
        qtrue=[dot(g['V'],truth) for g in obs]
        assert all(q>=1-e for q,e in zip(qtrue,EPS))
        # Direct rectangles are computed independently of incidence aggregation.
        for j,g in enumerate(obs):
            direct,q=law.direct_law(g['limit'],model['bins']);assert q==qtrue[j]
            for r,p in zip(g['records'],trueprob[j]):
                assert direct[r['label']]/q==p;checks['exact_probability_checks']+=1
        models.append({'lambda':lam,'true_q':qtrue,'true_D':[1-q for q in qtrue],
                       'true_S':trueS,'truth_mass':truth,'full_cell_mass':truefine})
        for nparent in (400,1600):
            for rep in range(reps):
                seed=[20260913,li,nparent,rep]
                s1,s2,s3=np.random.SeedSequence(seed).spawn(3)
                z=sample_parent(lam,nparent,np.random.default_rng(s1))
                which=np.random.default_rng(s2).integers(0,3,size=nparent)
                training=np.random.default_rng(s3).random(nparent)<0.5
                tr,va,delta,inc,labels=public_counts(z,which,training,obs)
                case=f'lambda{li}_N{nparent}_r{rep:03d}'
                (outdir/'latent_evaluation_only').mkdir(exist_ok=True)
                np.savez_compressed(outdir/'latent_evaluation_only'/f'{case}.npz',
                    z=z,limit_index=which,training=training,delta=delta,included=inc)
                spec={'name':case,'edges':EDGES,'limits':LIMITS,'targets':TARGETS,
                      'train':tr,'validation':va,'epsilon':['1/2']*3,'alpha':'1/20','training_steps':2}
                write_json(outdir/'observed_inputs'/f'{case}.json',spec)
                assert sum(map(sum,tr))+sum(map(sum,va))==int(inc.sum())
                checks['split_disjointness']+=1
                # Hidden values may be changed without changing observable input.
                changed=z.copy();changed[~delta]+=100
                tr2,va2,*_=public_counts(changed,which,training,obs)
                assert (tr2,va2)==(tr,va);checks['hidden_value_invariance']+=1
                # Check vectorized recorder against original exact-coordinate code.
                for k in range(min(24,nparent)):
                    lab=record(tuple(F.from_float(float(x)) for x in z[k]),
                               tuple(map(F,LIMITS[which[k]])),model['bins'])
                    assert lab==(tuple(map(int,labels[k])) if inc[k] else None)
                    checks['record_comparisons']+=1
                # Re-read public input; fitting doesn't receive truth or latent arrays.
                public=json.loads((outdir/'observed_inputs'/f'{case}.json').read_text())
                tr,va=public['train'],public['validation']
                if min(map(sum,tr))==0 or min(map(sum,va))==0:
                    rows.append({'case':case,'lambda':str(lam),'parent_n':nparent,'rep':rep,'status':'empty_stratum'})
                    continue
                predictor=init
                for _ in range(2):predictor,_,_=em_step(obs,tr,predictor)
                pred=conditional_probabilities(obs,predictor)
                logE=log_e_interval(pred,trueprob,va);crit=logint(1/ALPHA)
                if logE is None:member=True;decision='zero_prediction_noninformative'
                elif logE[1]<=crit[0]:member=True;decision='included'
                elif logE[0]>crit[1]:member=False;decision='excluded'
                else:member=None;decision='arithmetic_undecided'
                counts=[[a+b for a,b in zip(x,y)] for x,y in zip(tr,va)]
                m,hist=numeric_em(obs,counts,init,max_iter=4000,tol=1e-8)
                assert min(m)>=0
                # Reference curve: preserve observation-class totals through target splitting.
                lifted=np.array([m[k]/multiplicity[k] for k in jparent])
                cmat=np.array([[coeff[C[0]] for C in jclasses] for coeff in model['targets_coefficients']])
                Shat=cmat@lifted
                estimateprob=[]
                for g in obs:
                    q=float(np.dot(g['V'],m));estimateprob.append(np.array([np.dot(r['A'],m)/q for r in g['records']]))
                ntotal=sum(map(sum,counts))
                tv=sum(sum(cc)/ntotal*np.abs(p-np.array(pt,float)).sum()/2 for cc,p,pt in zip(counts,estimateprob,trueprob))
                selected=np.mean((z[inc,None,:]>np.array([[float(F(x)) for x in t] for t in TARGETS])[None,:,:]).all(axis=2),axis=0)
                full_empirical=np.mean((z[:,None,:]>np.array([[float(F(x)) for x in t] for t in TARGETS])[None,:,:]).all(axis=2),axis=0)
                tailok=all(float(np.dot(g['V'],m))>=float(1-e)-1e-12 for g,e in zip(obs,EPS))
                row={'case':case,'lambda':str(lam),'parent_n':nparent,'rep':rep,'status':'computed',
                     'A':int(np.sum(delta.all(axis=1))),'B':int(np.sum(delta[:,0]&~delta[:,1])),
                     'C':int(np.sum(~delta[:,0]&delta[:,1])),'D':int(np.sum(~inc)),
                     'train_n':sum(map(sum,tr)),'validation_n':sum(map(sum,va)),
                     'true_in_confidence_set':member,'coverage_decision':decision,
                     'EM_iterations':len(hist)-1,'EM_stop_tolerance_met':bool(hist[-1]['positive_score']<=1e-8),
                     'EM_positive_score':hist[-1]['positive_score'],'EM_tail_class_final':tailok,
                     'EM_global_optimality':'not_certified','EM_invisible_mass':float(m[invisible].sum()),
                     'true_invisible_mass':float(truth[invisible[0]]),
                     'observed_law_TV':float(tv),'reference_survival_RMSE':float(np.sqrt(np.mean((Shat-np.array(trueS,float))**2))),
                     'selected_oracle_RMSE':float(np.sqrt(np.mean((selected-np.array(trueS,float))**2))),
                     'full_parent_empirical_RMSE':float(np.sqrt(np.mean((full_empirical-np.array(trueS,float))**2))),
                     'raw_min_ell_increment':float(np.min(np.diff([h['ell'] for h in hist]))) if len(hist)>1 else 0.,
                     'max_mass_sum_error':float(max(abs(h['mass_sum']-1) for h in hist))}
                for k in range(len(TARGETS)):row[f'Shat_{k}']=float(Shat[k]);row[f'Strue_{k}']=float(trueS[k])
                rows.append(row)
                write_json(outdir/'fits'/f'{case}.json',{'input_file':f'observed_inputs/{case}.json','seed':seed,
                    'prediction_mass':predictor,'prediction_laws':pred,
                    'truth_membership_logE_interval':logE,'log_rejection_threshold_interval':crit,
                    'EM_mass':m.tolist(),'reference_target_mass':lifted.tolist(),
                    'history':hist,'evaluation':row})
            print('completed lambda',lam,'parent_n',nparent,'reps',reps,flush=True)
    keys=sorted(set().union(*(r.keys() for r in rows)))
    with (outdir/'trials.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    summary=[]
    for lam in ('-3/4','0','3/4'):
        for n in (400,1600):
            rr=[r for r in rows if r['lambda']==lam and r['parent_n']==n]
            valid=[r for r in rr if r['status']=='computed']
            s={'lambda':lam,'parent_n':n,'trials':len(rr),'computed':len(valid),
               'parent_covered':sum(r['true_in_confidence_set'] is True for r in valid),
               'parent_rejected':sum(r['true_in_confidence_set'] is False for r in valid),
               'membership_undecided':sum(r['true_in_confidence_set'] is None for r in valid),
               'EM_stopped':sum(r['EM_stop_tolerance_met'] for r in valid),
               'EM_outside_tail':sum(not r['EM_tail_class_final'] for r in valid)}
            for k in ['A','B','C','D','train_n','validation_n','EM_iterations','observed_law_TV',
                      'reference_survival_RMSE','selected_oracle_RMSE','full_parent_empirical_RMSE','EM_invisible_mass','true_invisible_mass']:
                s['mean_'+k]=float(np.mean([r[k] for r in valid])) if valid else None
            summary.append(s)
    with (outdir/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    result={'stage':'22A','scope':'continuous-parent sampling, observed-law check, confidence membership; projection separate',
            'software':{'python':platform.python_version(),'numpy':np.__version__},
            'parent_n_options':[400,1600],'replicates':reps,'total_trials':len(rows),
            'full_cells':len(model['cells']),'observation_classes':len(classes),'target_classes':len(jclasses),
            'records_per_layer':[len(g['records']) for g in obs],
            'checks':checks,'truths':models,'summary':summary,
            'no_hidden_truth_given_to_fit':True,'all_record_categories_retained':True,
            'new_point_global_certification':False,'new_sharp_projection_for_every_trial':False,
            'seconds':time.perf_counter()-start}
    write_json(outdir/'stage22_sampling_results.json',result)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs');p.add_argument('--reps',type=int,default=20)
    args=p.parse_args();print(json.dumps(encode(run(args.outdir,args.reps)),ensure_ascii=False,indent=2))
