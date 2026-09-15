#!/usr/bin/env python3
"""Recompute public records and diagnostics; optionally compare a fresh rerun.
This is a float64 validation of the new benchmark, NOT a global certificate.
"""
from pathlib import Path
import sys,json,csv,argparse,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bivsurv.model import GridModel
from bivsurv.em import quantities,fit_em
from bivsurv.benchmark import save_json,gaussian_cell_masses,metrics,build_design,write_csv

def verify(root,rerun=None):
    root=Path(root);config=json.loads((root/'run_config.json').read_text())
    tested={'catalogues':0,'public_tables':0,'fit_invariants':0,'metric_comparisons':0,'dense_step_comparisons':0}
    maxgrad=0.;maxdense=0.;quaddiff=0.;margdiff=0.
    for name in config['limits_standardized_X']:
        model,initial,we,xe=build_design(config,name)
        A=model.dense_operators();V=model.visible.astype(float)
        if name==next(iter(config['limits_standardized_X'])):
            from scipy.special import ndtr
            for ri,rho in enumerate(config['rho']):
                mx,diag=gaussian_cell_masses(we,rho,epsabs=1e-14)
                t,_=gaussian_cell_masses(we,rho,epsabs=1e-13)
                quaddiff=max(quaddiff,float(abs(mx-t).max()))
                margdiff=max(margdiff,float(abs(mx.sum(axis=1)-np.diff(np.r_[0,ndtr(we),1])).max()))
        for p in sorted((root/'public_counts'/name).glob('*.json')):
            public=json.loads(p.read_text());case=p.stem
            with np.load(root/'latent_evaluation_only'/f'{case}.npz') as d:
                x,field,training=[d[k].copy() for k in ('x','field','training')]
            inc,delta,labels,idx=model.record_catalogue(-x,field)
            for key,subset in [('counts',None),('train',training),('validation',~training)]:
                found=model.counts_from_indices(field,idx,subset)
                for a,b in zip(found,public[key]):np.testing.assert_array_equal(a,b)
                tested['public_tables']+=1
            # Nondetections and unrecorded objects may be moved arbitrarily further below detection.
            changed=-x.copy();changed[~delta]+=123
            after=model.record_catalogue(changed,field)
            for a,b in zip((inc,delta,labels,idx),after):np.testing.assert_array_equal(a,b)
            tested['catalogues']+=1
            counts=[np.asarray(a) for a in public['counts']]
            with np.load(root/'fits'/name/f'{case}.npz') as d:
                mass,history,truth=[d[k].copy() for k in ('mass','history','truth')]
                np.testing.assert_allclose(mass.sum(),1,atol=2e-12)
                assert mass.min()>=0 and np.min(np.diff(history[:,1]))>=-1e-11
                assert np.max(abs(history[:,4]-initial[model.invisible].sum()))<1e-12
                _,ell,g,q,K=quantities(model,counts,mass)
                assert abs(ell-history[-1,1])<1e-12
                assert np.max(g)<=config['em']['score_tolerance']*1.001
                maxgrad=max(maxgrad,float(g.max()))
                assert np.min(q)>=.5-1e-10
                tested['fit_invariants']+=1
                result=metrics(model,mass,truth,counts,we)
                saved=json.loads((root/'fits'/name/f'{case}.json').read_text())['evaluation'][0]
                for key,value in result.items():
                    assert abs(value-saved[key])<1e-12;tested['metric_comparisons']+=1
            # Compare the optimized indexed map against the original dense update, at initialization.
            aa=[a@initial for a in A];qq=V@initial;njs=np.array([n.sum() for n in counts])
            num=np.zeros(model.H)
            for a,y,n in zip(A,aa,counts):
                r=np.zeros_like(y);pos=n>0;r[pos]=n[pos]/y[pos];num+=r@a
            complete=initial*(num+(njs/qq)@(1-V))/(njs/qq).sum()
            indexed,*_=quantities(model,counts,initial)
            diff=float(np.max(abs(indexed-complete)));maxdense=max(maxdense,diff);assert diff<1e-13
            tested['dense_step_comparisons']+=1
    comparison={}
    if rerun:
        rerun=Path(rerun);json_count=array_count=0
        for p in (root/'public_counts').rglob('*.json'):
            assert p.read_bytes()==(rerun/p.relative_to(root)).read_bytes();json_count+=1
        for directory in ['latent_evaluation_only','fits','models']:
            for p in (root/directory).rglob('*.npz'):
                with np.load(p) as a,np.load(rerun/p.relative_to(root)) as b:
                    assert a.files==b.files
                    for key in a.files:np.testing.assert_array_equal(a[key],b[key]);array_count+=1
        def deterministic_rows(path):
            return [{k:v for k,v in r.items() if k not in ('seconds',)} for r in csv.DictReader(path.open())]
        assert deterministic_rows(root/'all_metrics.csv')==deterministic_rows(rerun/'all_metrics.csv')
        comparison={'public_JSON_byte_matches':json_count,'array_exact_matches':array_count,'all_metric_rows_equal_except_seconds':True}
    result={'passed':True,**tested,'max_score':maxgrad,'max_dense_indexed_update_difference':maxdense,
            'quadrature_tightening_max_difference':quaddiff,'true_marginal_max_difference':margdiff,
            'fresh_rerun_comparison':comparison,'global_certificate_for_new_fits':False}
    save_json(root/'verification.json',result);print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='outputs/visual');p.add_argument('--compare-rerun');a=p.parse_args();verify(a.output,a.compare_rerun)
