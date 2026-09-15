"""Reproducibility and deterministic arithmetic controls, separate from search."""
from pathlib import Path
import json,sys,csv,hashlib
from fractions import Fraction as F
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE))
from global_bound import prepare,loglik_interval
from verify_stage25 import verify_case
from exact import dot,logint,encode

def main():
    results=[];all_same=True;ratio_count=0;candidate_count=0;scale_checks=0
    for path in sorted((BASE/'outputs/main_certificates').glob('coarse_*.json')):
        c=json.loads(path.read_text());other=json.loads((BASE/'outputs/rerun_certificates'/path.name).read_text())
        aa={k:v for k,v in c.items() if k!='seconds'};bb={k:v for k,v in other.items() if k!='seconds'}
        assert aa==bb
        source=json.loads((BASE/'inputs'/path.name).read_text());d=prepare(source)
        results.append({'case':c['case'],'N':d['N'],'source_first_order_upper':source['candidate']['stationarity']['directional_gap_upper'],
          'nodes':len(c['tree']),'leaves':sum('children' not in n for n in c['tree']),
          'mean_loglik_lower':c['mean_lower_interval'][0],'global_mean_loglik_upper':c['mean_upper'],
          'global_mean_gap_upper':c['gap_upper'],'global_mean_gap_upper_float':float(F(c['gap_upper'])),
          'global_total_gap_upper_float':float(F(c['gap_upper'])*d['N']),
          'GEM80_global_mean_gap_lower':source['candidate']['ascent_from_GEM'][0],
          'GEM80_global_mean_gap_upper':str(F(source['candidate']['ascent_from_GEM'][1])+F(c['gap_upper'])),
          'GEM80_global_mean_gap_upper_float':float(F(source['candidate']['ascent_from_GEM'][1])+F(c['gap_upper'])),
          'precision_achieved':c['precision_achieved'],'seconds':c['seconds']})
        full=source['compiled'];H=len(full['initial']);initial=list(map(F,full['initial']));ids=d['visible_ids']
        for trial in range(24):
            z=[F(1+((trial+1)*(h+3)+7*h*h)%23) for h in range(H)];z=[v/sum(z) for v in z]
            m=[F(3,4)*x+F(1,4)*y for x,y in zip(initial,z)]
            assert min(dot(v,m) for v in full['V'])>=F(1,2)
            r=dot(full['V'][2],m);s=[m[h]/r for h in ids]
            assert sum(s)==1 and min(s)>=0
            assert all(dot(v,s)>=dot(vv,m) for v,vv in zip(d['V'],full['V']));scale_checks+=3
            for a,v,af,vf in zip(d['A'],d['VR'],full['A'],full['VR']):
                assert dot(a,s)/dot(v,s)==dot(af,m)/dot(vf,m);ratio_count+=1
            li=loglik_interval(d,s);assert li[1]<=F(c['mean_upper']);candidate_count+=1
    out=BASE/'outputs'
    with (out/'global_gap_comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    # The tight and budget-limited probes are valid certificates, but NOT precise.
    controls=[]
    for kind in ('budget_control','strict_control'):
        p=next((out/kind).glob('coarse_*.json'));r=verify_case(p)
        assert r['precision_achieved'] is False;controls.append({'kind':kind,**r})
    # Exact likelihood fibre of the saved interior candidate, without changing it.
    name='coarse_d3__lambda0_N400_r010'
    rec=json.loads((BASE/'inputs'/(name+'.json')).read_text());d=prepare(rec);full=rec['compiled']
    vmax=d['ref'];pmin=min(dot(v,vmax) for v in d['V']);theta_max=1-F(1,2)/pmin
    assert theta_max>0
    invisible=[h for h in range(len(full['initial'])) if h not in d['visible_ids']];assert len(invisible)==1
    parents=[]
    for th in (F(0),theta_max/2,theta_max):
        mass=[F(0)]*len(full['initial'])
        for h,x in zip(d['visible_ids'],vmax):mass[h]=(1-th)*x
        mass[invisible[0]]=th
        assert sum(mass)==1 and min(mass)>=0 and min(dot(v,mass) for v in full['V'])>=F(1,2)
        assert all(dot(a,mass)/dot(v,mass)==dot(aa,vmax)/dot(vv,vmax) for a,v,aa,vv in zip(full['A'],full['VR'],d['A'],d['VR']))
        parents.append({'theta':th,'mass':mass})
    (out/'same_likelihood_family.json').write_text(json.dumps(encode({'case':name,'theta_max':theta_max,'theta_max_float':float(theta_max),'target':'P(D_deep)=S(d,d)','parents':parents}),indent=2)+'\n')
    summary={'passed':True,'main_cases':13,'all_precision_achieved':all(r['precision_achieved'] for r in results),
        'nodes':sum(r['nodes'] for r in results),'leaves':sum(r['leaves'] for r in results),
        'maximum_global_mean_gap_upper':max(r['global_mean_gap_upper_float'] for r in results),
        'maximum_total_gap_upper':max(r['global_total_gap_upper_float'] for r in results),
        'GEM80_minimum_global_gap_lower':min(float(F(r['GEM80_global_mean_gap_lower'])) for r in results),
        'GEM80_maximum_global_gap_upper':max(r['GEM80_global_mean_gap_upper_float'] for r in results),
        'deterministic_feasible_candidates_checked':candidate_count,'exact_ratio_invariance_checks':ratio_count,
        'selection_probability_monotonicity_checks':scale_checks,'controls':controls,
        'all_13_searches_reexecuted':True,'certificate_fields_identical_except_seconds':True,
        'unchanged_stage24_candidate_masses':True,'new_random_samples':False,'GEM_refitting':False,
        'confidence_sets_recomputed':False,'same_likelihood_theta_max':float(theta_max)}
    (out/'audit_results.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='controls'},indent=2))
if __name__=='__main__':main()
