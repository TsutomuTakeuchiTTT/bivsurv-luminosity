"""New Gaussian luminosity benchmark; simulation truth is evaluator-only."""
from __future__ import annotations
from pathlib import Path
import json, csv, time, platform, hashlib
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr, expit
from .model import GridModel
from .em import fit_em


def save_json(path, obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)+'\n')

def write_csv(path, rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def build_design(config, name):
    spec=config['standardized_edges']
    w_edges=np.linspace(spec['lower'],spec['upper'],int(round((spec['upper']-spec['lower'])/spec['step']))+1)
    mean=np.asarray(config['mean_log10L']);sd=np.asarray(config['sigma_dex'])
    x_edges=tuple(mean[b]+sd[b]*w_edges for b in (0,1))
    z_edges=tuple(-e[::-1] for e in x_edges)
    wlimits=np.asarray(config['limits_standardized_X'][name])
    # Use declared edge indices rather than recalculate nearly equal floating values.
    li=[[int(np.flatnonzero(w_edges==t)[0]) for t in row] for row in wlimits]
    limits_x=np.array([[x_edges[b][row[b]] for b in (0,1)] for row in li])
    model=GridModel(z_edges,-limits_x)
    # Independent standard logistic reference in standardized X; it is not Gaussian.
    probabilities=np.diff(np.r_[0.,expit(w_edges),1.])
    reference=np.outer(probabilities,probabilities)[::-1,::-1].ravel().copy()
    return model, reference, w_edges, x_edges


def normal_interval(a,b):
    if a>0:
        return ndtr(-a)-ndtr(-b)
    return ndtr(b)-ndtr(a)


def gaussian_cell_masses(w_edges,rho,epsabs=1e-13):
    """Deterministic 1-D conditional Gaussian integrals, including infinite tails.
    QUADPACK error estimates are numerical diagnostics, not rigorous intervals.
    Returns masses ordered in increasing X, not Z.
    """
    if not -1<rho<1:raise ValueError('Nondegenerate correlation required')
    edges=np.r_[-np.inf,w_edges,np.inf]
    mass=np.empty((len(edges)-1,)*2);error=np.zeros_like(mass)
    scale=np.sqrt(1-rho*rho)
    for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        for j,(c,d) in enumerate(zip(edges[:-1],edges[1:])):
            def integrand(t):
                return np.exp(-t*t/2)/np.sqrt(2*np.pi)*normal_interval((c-rho*t)/scale,(d-rho*t)/scale)
            mass[i,j],error[i,j]=quad(integrand,a,b,epsabs=epsabs,epsrel=1e-10,limit=200)
    if np.any(mass<0) or abs(mass.sum()-1)>1e-10:
        raise ArithmeticError('Truth quadrature failed')
    return mass,{'total_mass':float(mass.sum()),'sum_error_estimates':float(error.sum()),'max_error_estimate':float(error.max()),'epsabs':epsabs}


def bright_surface(mx):
    """Probability of both X components exceeding a bin lower boundary."""
    return mx[::-1,::-1].cumsum(0).cumsum(1)[::-1,::-1]


def metrics(model, mass, true_mass, counts, w_edges):
    vis=~model.invisible
    qhat=mass[vis].sum();qt=true_mass[vis].sum()
    qm=mass/qhat;qm[~vis]=0
    q0=true_mass/qt;q0[~vis]=0
    # Every finite bin in the deepest both-detected square has an identifiable Q mass.
    resolved=model.common_detected
    l1resolved=float(np.abs(qm[resolved]-q0[resolved]).sum())
    lower=np.r_[-np.inf,w_edges]
    metric_indices=np.flatnonzero((lower>=-0.5)&(lower<=2.5))
    mx=mass.reshape(model.shape)[::-1,::-1]
    truex=true_mass.reshape(model.shape)[::-1,::-1]
    cond=qm.reshape(model.shape)[::-1,::-1]
    cond0=q0.reshape(model.shape)[::-1,::-1]
    err=(bright_surface(cond)-bright_surface(cond0))[np.ix_(metric_indices,metric_indices)]
    errp=(bright_surface(mx)-bright_surface(truex))[np.ix_(metric_indices,metric_indices)]
    a,qq=model.probabilities(mass);a0,qq0=model.probabilities(true_mass)
    nt=sum(n.sum() for n in counts)
    tv=sum(n.sum()/nt*0.5*np.abs(aj/qj-aj0/qj0).sum() for n,aj,qj,aj0,qj0 in zip(counts,a,qq,a0,qq0))
    return {'resolved_bin_L1_Q':l1resolved,'bright_probability_RMSE_Q':float(np.sqrt(np.mean(err*err))),
            'bright_probability_RMSE_P_reference':float(np.sqrt(np.mean(errp*errp))),
            'record_law_TV':float(tv),'invisible_mass':float(mass[~vis].sum()),
            'true_invisible_mass':float(true_mass[~vis].sum()),'Q_mass_resolved':float(qm[resolved].sum())}


def evaluate_split(model,true_mass,train,val,alpha):
    # Add-one is a record predictor only; never modifies the EM mass update.
    pred=[(n+1)/(n.sum()+len(n)) for n in train]
    a,q=model.probabilities(true_mass)
    loge=0.
    for n,g,aj,qj in zip(val,pred,a,q):
        pos=n>0
        if np.any(aj[pos]<=0):raise ArithmeticError('Observed zero-probability truth bin')
        loge+=float(n[pos]@(np.log(g[pos])-np.log(aj[pos]/qj)))
    return {'true_membership':bool(loge<=np.log(1/alpha)), 'logE_true':loge,
            'distance_to_threshold':float(np.log(1/alpha)-loge),'prediction':'stratumwise_add_one',
            'arithmetic':'float64 evaluation; not a rational certificate'}


def run(config_path,output,replicates=None,sizes=None):
    config=json.loads(Path(config_path).read_text());out=Path(output);out.mkdir(parents=True,exist_ok=True)
    if replicates is not None:config['replicates']=int(replicates)
    if sizes is not None:config['parent_sizes']=list(map(int,sizes))
    save_json(out/'run_config.json',config)
    designs={name:build_design(config,name) for name in config['limits_standardized_X']}
    w_edges=next(iter(designs.values()))[2]
    true={};quadrature={}
    for ri,rho in enumerate(config['rho']):
        mx,diag=gaussian_cell_masses(w_edges,rho)
        true[ri]=mx[::-1,::-1].ravel().copy();quadrature[str(rho)]=diag
    save_json(out/'truth_quadrature.json',quadrature)
    for name,(model,initial,_,xe) in designs.items():
        save_json(out/'models'/f'{name}.json',model.as_dict())
        np.savez_compressed(out/'models'/f'{name}_arrays.npz',initial=initial,x_edges1=xe[0],x_edges2=xe[1],w_edges=w_edges)
    rows=[];start=time.perf_counter();tot=0
    mean=np.asarray(config['mean_log10L']);sd=np.asarray(config['sigma_dex'])
    for ri,rho in enumerate(config['rho']):
        for nparent in config['parent_sizes']:
            for rep in range(config['replicates']):
                rng_z,rng_f,rng_s=[np.random.default_rng(s) for s in np.random.SeedSequence([config['seed'],ri,nparent,rep]).spawn(3)]
                u=rng_z.standard_normal((nparent,2))
                wx=np.column_stack((u[:,0],rho*u[:,0]+np.sqrt(1-rho*rho)*u[:,1]))
                x=mean+sd*wx;z=-x
                field=rng_f.choice(3,size=nparent,p=config['field_probabilities']);training=rng_s.random(nparent)<0.5
                case=f'rho{rho:.2f}_N{nparent}_r{rep:03d}'
                (out/'latent_evaluation_only').mkdir(exist_ok=True)
                np.savez_compressed(out/'latent_evaluation_only'/f'{case}.npz',x=x,field=field,training=training)
                previous_inc=None
                for name,(model,initial,_,xe) in designs.items():
                    inc,delta,labels,idx=model.record_catalogue(z,field)
                    if previous_inc is not None and np.any(previous_inc & ~inc):
                        raise AssertionError('Deep selection must contain shallow selection')
                    previous_inc=inc
                    counts=model.counts_from_indices(field,idx)
                    tr=model.counts_from_indices(field,idx,training);va=model.counts_from_indices(field,idx,~training)
                    public={'case':case,'design':name,'counts':[n.tolist() for n in counts],
                            'train':[n.tolist() for n in tr],'validation':[n.tolist() for n in va]}
                    save_json(out/'public_counts'/name/f'{case}.json',public)
                    # Read back ONLY the public record counts. Neither latent values nor truth passed.
                    public=json.loads((out/'public_counts'/name/f'{case}.json').read_text())
                    counts=[np.asarray(n) for n in public['counts']]
                    t=time.perf_counter()
                    full=fit_em(model,counts,initial,**config['em'],correction=True,tail_upper=config['tail_upper_bound_diagnostic'])
                    elapsed=time.perf_counter()-t
                    t=time.perf_counter()
                    noD=fit_em(model,counts,initial,**config['em'],correction=False)
                    elapsed_noD=time.perf_counter()-t
                    # Deliberately naive both-detected histogram, from PUBLIC A records only.
                    naive=np.zeros(model.H)
                    for j,(lab,nn,code) in enumerate(zip(model.labels,counts,model.record_index)):
                        for r in np.flatnonzero((lab[:,0]==1)&(lab[:,1]==1)):
                            hs=np.flatnonzero(code==r)
                            if len(hs)!=1:raise AssertionError('A record must identify one cell in aligned design')
                            naive[hs[0]]+=nn[r]
                    if naive.sum()==0:raise RuntimeError('No A records: cannot evaluate A-only histogram')
                    naive/=naive.sum()
                    membership=evaluate_split(model,true[ri],tr,va,config['alpha_split_membership'])
                    region={s:int(np.sum((delta[:,0]==a)&(delta[:,1]==b))) for s,a,b in [('A',1,1),('B',1,0),('C',0,1),('D',0,0)]}
                    base={'case':case,'rho':rho,'design':name,'parent_n':nparent,'replicate':rep,'observed_n':int(inc.sum()),**region,
                          'train_n':int(sum(n.sum() for n in tr)),'validation_n':int(sum(n.sum() for n in va)),**membership}
                    for method,fit,m,secs in [('D_completion',full,full.mass,elapsed),('no_D_ablation',noD,noD.mass,elapsed_noD),('A_only_histogram',None,naive,0.)]:
                        row={**base,'method':method,**metrics(model,m,true[ri],counts,w_edges),'seconds':secs,
                             'iterations':fit.iterations if fit else 0,'score':fit.score if fit else 0.,
                             'status':fit.status if fit else 'deliberately_uncorrected_histogram',
                             'tail_feasible':fit.tail_feasible if fit else None,
                             'min_mean_loglik_gain':fit.min_loglik_increment if fit else None,
                             'invisible_drift':fit.invisible_mass_drift if fit else None}
                        rows.append(row)
                    folder=out/'fits'/name;folder.mkdir(parents=True,exist_ok=True)
                    np.savez_compressed(folder/f'{case}.npz',mass=full.mass,history=full.history,noD=noD.mass,noD_history=noD.history,naive=naive,truth=true[ri])
                    save_json(folder/f'{case}.json',{'evaluation':rows[-3:], 'source_public_counts_sha256':hashlib.sha256((out/'public_counts'/name/f'{case}.json').read_bytes()).hexdigest()})
                    tot+=1
                if rep%5==0:print(f'rho={rho} N={nparent} rep={rep} completed={tot}',flush=True)
    write_csv(out/'all_metrics.csv',rows)
    summary=[]
    for rho in config['rho']:
        for name in designs:
            for n in config['parent_sizes']:
                for method in ('D_completion','no_D_ablation','A_only_histogram'):
                    chosen=[r for r in rows if (r['rho'],r['design'],r['parent_n'],r['method'])==(rho,name,n,method)]
                    item={'rho':rho,'design':name,'parent_n':n,'method':method,'replicates':len(chosen)}
                    for metric in ['observed_n','A','B','C','D','resolved_bin_L1_Q','bright_probability_RMSE_Q','bright_probability_RMSE_P_reference','record_law_TV','seconds','iterations']:
                        v=np.array([r[metric] for r in chosen]);item[metric+'_mean']=float(v.mean());item[metric+'_sd']=float(v.std(ddof=1)) if len(v)>1 else 0.
                    item['budget_exhausted']=sum(r['status']=='iteration_budget_exhausted' for r in chosen)
                    item['tail_violations']=sum(r['tail_feasible'] is False for r in chosen)
                    item['true_membership_count']=sum(r['true_membership'] for r in chosen)
                    summary.append(item)
    write_csv(out/'summary.csv',summary)
    import scipy,matplotlib
    manifest={'new_parent_samples':len(config['rho'])*len(config['parent_sizes'])*config['replicates'],
              'new_parent_objects':len(config['rho'])*sum(config['parent_sizes'])*config['replicates'],
              'paired_design_evaluations':tot,'EM_fits':2*tot,'elapsed_seconds':time.perf_counter()-start,
              'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,
              'estimator':'indexed D-completion EM; matched against recovered Fraction implementation',
              'global_optimality_certificates':'Not computed for these new high-resolution fits',
              'random_model':'Gaussian parent; deterministic quadrature truth; selected latent values not used as parent truth'}
    save_json(out/'execution.json',manifest)
    return manifest
