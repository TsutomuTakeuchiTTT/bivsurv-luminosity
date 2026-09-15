#!/usr/bin/env python3
"""Additional deterministic reference-initialization sensitivity checks."""
from pathlib import Path
import sys,argparse,json
import numpy as np
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bivsurv.benchmark import build_design,write_csv
from bivsurv.em import fit_em

def run(root):
 root=Path(root);config=json.loads((root/'run_config.json').read_text());rows=[]
 for rho in config['rho']:
  case=f"rho{rho:.2f}_N{config['representative']['parent_size']}_r000"
  for name in config['limits_standardized_X']:
   model,initial,we,xe=build_design(config,name)
   public=json.loads((root/'public_counts'/name/f'{case}.json').read_text());counts=[np.asarray(n) for n in public['counts']]
   with np.load(root/'fits'/name/f'{case}.npz') as a:ref=a['mass'].copy()
   vr=~model.invisible;qr=ref.copy();qr[~vr]=0;qr/=qr.sum()
   for shift in [-1.,0.,1.]:
    f=np.diff(np.r_[0,expit(we-shift),1]);init=np.outer(f,f)[::-1,::-1].ravel().copy()
    fit=fit_em(model,counts,init,score_tolerance=1e-9,tail_upper=.5)
    q=fit.mass.copy();q[~vr]=0;q/=q.sum()
    def laws(m):
     a,b=model.probabilities(m);return np.concatenate([x/y for x,y in zip(a,b)])
    rows.append({'rho':rho,'design':name,'shift_of_reference_in_sigma_units':shift,'initial_invisible':float(init[~vr].sum()),'final_invisible':float(fit.mass[~vr].sum()),'max_resolved_Q_mass_difference_from_default':float(abs(q[model.common_detected]-qr[model.common_detected]).max()),'max_record_law_difference':float(abs(laws(fit.mass)-laws(ref)).max()),'mean_loglik':float(fit.history[-1,1]),'score':fit.score,'iterations':fit.iterations,'tail_feasible':fit.tail_feasible})
 write_csv(root/'initialization_sensitivity.csv',rows)
 print('Completed',len(rows),'initialization checks.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default='outputs/visual');a=p.parse_args();run(a.output)
