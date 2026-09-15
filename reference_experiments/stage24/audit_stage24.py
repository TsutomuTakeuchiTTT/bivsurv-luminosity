"""Run the Stage24 predeclared 13-case constrained-GEM audit.
No new parent draws, no training refits, no replacement of Stage23 records.
Independent boundary optimization is separately labeled, not a GEM step.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from core import *
import time,hashlib,csv,itertools,platform

def save(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(encode(x),ensure_ascii=False,sort_keys=True,indent=2)+'\n')
def csvsave(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run(out,steps=80):
 out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();rows=[];active_counts={};cert_count=0
 names=json.loads((BASE/'inputs/selection_manifest.json').read_text())['names']
 for name in names:
  inp_path=BASE/'inputs/observed'/f'{name}.json';inp=json.loads(inp_path.read_text())
  old=json.loads((BASE/'inputs/prior_fits'/f'{name}.json').read_text());d=compile_model(inp)
  if encode(d['initial'])!=old['initial_mass']:raise AssertionError('Different initial measure')
  m=d['initial'];hist=[];max_gap=F(0);min_ll=None;min_q=None
  for k in range(steps):
   w,K=expected(d,m);raw,cert=closed_mstep(w,d['sig']);u,mix=rounded_candidate(raw,d,d['initial'])
   sg=surrogate_interval(w,m,u);lg=gain_interval(d,m,u)
   if sg[0]<0 or lg[0]<0:raise ArithmeticError(f'No certified ascent: {name} step{k}')
   if lg[0]<K*sg[1]/d['N']:raise ArithmeticError('Minorization check inconclusive')
   # A concave tangent at the accepted proposal bounds its surrogate loss vs exact M-step.
   mgap=sum((a*(b/c-1) for a,b,c in zip(w,raw,u) if a),F(0))
   assert mgap>=0
   max_gap=max(max_gap,mgap);min_ll=lg[0] if min_ll is None else min(min_ll,lg[0]);min_q=min(qvals(d,u)) if min_q is None else min(min_q,min(qvals(d,u)))
   hist.append({'iteration':k+1,'mass':u,'active':cert['active'],'kkt_eta':cert['eta'],'kkt_tail_multipliers':cert['multipliers'],
      'normalized_surrogate_gain':sg,'mean_likelihood_gain':lg,'exact_Mstep_surrogate_loss_upper':mgap,'rounding_mix':mix})
   active_counts[cert['active']]=active_counts.get(cert['active'],0)+1;m=u;cert_count+=1
  gemstat=stationary_certificate(d,m)
  candidate,stat,search=independent_candidates(d,m)
  extra=gain_interval(d,m,candidate)
  assert extra[0]>=0
  assert stat['directional_gap_upper']<F(1,10**8), (name,float(stat['directional_gap_upper']))
  # Earlier raw floating fit is rationalized ONLY for a diagnostic, not overwritten.
  om=tuple(F.from_float(float(x)) for x in old['final_mass']);om=tuple(x/sum(om) for x in om)
  iq=[h for h in range(len(om)) if all(not v[h] for v in d['V'])];assert len(iq)==1
  h=iq[0];theta=om[h];noI=tuple(F(0) if k==h else x/(1-theta) for k,x in enumerate(om))
  oldq=qvals(d,om);newq=qvals(d,candidate);gemq=qvals(d,m)
  rawdiff=gain_interval(d,candidate,om)
  trueS=[old['evaluation'][f'Strue_{k}'] for k in range(3)]
  refnew=[float(dot(c,candidate)) for c in d['coeff']];refgem=[float(dot(c,m)) for c in d['coeff']]
  rows.append({'name':name,'observed_n':d['N'],'mass_count':len(m),'old_min_q':float(min(oldq)),
   'old_theta':float(theta),'old_no_invisible_min_q':float(min(qvals(d,noI))),
   'invisible_removal_alone_feasible':feasible(d,noI),
   'certified_GEM_steps':steps,'GEM_min_q':float(min(gemq)),'GEM_directional_gap_upper':float(gemstat['directional_gap_upper']),
   'candidate_min_q':float(min(newq)),'candidate_theta':float(candidate[h]),
   'candidate_raw_positive_score':float(max(F(0),stat['maximum_raw_score'])),
   'candidate_directional_gap_upper':float(stat['directional_gap_upper']),
   'mean_likelihood_gain_after_GEM_lower':float(extra[0]),
   'old_minus_constrained_mean_loglik_lower':float(rawdiff[0]),
   'old_minus_constrained_mean_loglik_upper':float(rawdiff[1]),
   'old_reference_RMSE':old['evaluation']['reference_RMSE'],
   'GEM_reference_RMSE':float(np.sqrt(np.mean((np.array(refgem)-trueS)**2))),
   'candidate_reference_RMSE':float(np.sqrt(np.mean((np.array(refnew)-trueS)**2)))})
  output={'case':name,'input_sha256':hashlib.sha256(inp_path.read_bytes()).hexdigest(),
   'old_fit_sha256':hashlib.sha256((BASE/'inputs/prior_fits'/f'{name}.json').read_bytes()).hexdigest(),
   'scope':'80 certified feasible generalized EM updates plus separately labeled independent stationary likelihood candidate',
   'compiled':{'A':d['A'],'VR':d['VR'],'counts':d['counts'],'V':d['V'],'initial':d['initial'],'sig':d['sig'],'N':d['N']},
   'GEM_steps':hist,'GEM_final_stationarity':gemstat,'minimum_mean_likelihood_gain':min_ll,
   'maximum_normalized_Mstep_loss_upper':max_gap,
   'independent_candidate':{'mass':candidate,'stationarity':stat,'ascent_from_GEM':extra,'numerical_searches':search,
                            'global_maximum_certified':False,'is_GEM_step':False},'comparison':rows[-1]}
  save(out/'certificates'/f'{name}.json',output)
  print(name,'GEM',steps,'q>=',float(min_q),'stationary',float(stat['directional_gap_upper']), 'elapsed',time.perf_counter()-start,flush=True)
 # Exact four-group KKT coverage, not random astrophysical data.
 cats={};comparison=0;max_numeric=0.
 for vals in itertools.product(range(1,5),repeat=4):
  w=tuple(F(v,sum(vals)) for v in vals);p,c=closed_mstep(w,(0,1,2,3));cats[c['active']]=cats.get(c['active'],0)+1
  if comparison<16 and len(cats)==4:
   WW=np.array(w,float);D=np.array([[0,0,1,1],[0,1,0,1]],float)
   rr=minimize(lambda z:-float(WW@np.log(z)),np.ones(4)/4,jac=lambda z:-WW/z,method='SLSQP',bounds=[(1e-14,1)]*4,
    constraints=[{'type':'eq','fun':lambda z:z.sum()-1,'jac':lambda z:np.ones(4)},
      {'type':'ineq','fun':lambda z:.5-D@z,'jac':lambda z:-D}],options={'ftol':1e-14,'maxiter':300})
   err=abs(float(WW@np.log(np.array(p,float)))-float(WW@np.log(rr.x)));max_numeric=max(max_numeric,err);comparison+=1
 result={'stage':24,'selected_violating_cases':len(names),'certified_GEM_updates':cert_count,'GEM_active_sets':active_counts,
  'all_steps_feasible_and_certified_increasing':True,
  'independent_stationary_candidates':len(rows),'max_stationarity_upper':max(r['candidate_directional_gap_upper'] for r in rows),
  'global_likelihood_maxima_certified':False,'fixed_step_GEM_claimed_converged':False,
  'invisible_removal_alone_feasible_count':sum(r['invisible_removal_alone_feasible'] for r in rows),
  'exact_four_group_tests':256,'four_group_active_sets':cats,'numeric_Mstep_comparisons':comparison,'max_numeric_Mstep_objective_difference':max_numeric,
  'new_parent_draws':0,'training_predictors_changed':False,'confidence_sets_changed':False,'new_coverage_evaluations':0,
  'software':{'python':platform.python_version(),'numpy':np.__version__},'seconds':time.perf_counter()-start}
 csvsave(out/'comparison.csv',rows);save(out/'results.json',result);print(result,flush=True)
 return result

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs');p.add_argument('--steps',type=int,default=80)
 a=p.parse_args();run(a.outdir,a.steps)
