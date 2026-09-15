#!/usr/bin/env python3
"""Fit from public counts, geometry and an initial mass ONLY."""
from pathlib import Path
import sys,json,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bivsurv import GridModel,fit_em
from bivsurv.benchmark import save_json
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--model',required=True);p.add_argument('--input',required=True);p.add_argument('--initial',required=True);p.add_argument('--output',required=True)
p.add_argument('--score-tolerance',type=float,default=1e-7);p.add_argument('--max-iterations',type=int,default=20000)
a=p.parse_args()
model=GridModel.from_dict(json.loads(Path(a.model).read_text()))
public=json.loads(Path(a.input).read_text())
with np.load(a.initial) as d:initial=d['initial'].copy()
f=fit_em(model,[np.asarray(n) for n in public['counts']],initial,score_tolerance=a.score_tolerance,max_iterations=a.max_iterations,tail_upper=.5)
save_json(a.output,{'mass':f.mass.tolist(),'status':f.status,'iterations':f.iterations,'score':f.score,'q':f.q.tolist(),'tail_feasible_diagnostic':f.tail_feasible,'global_optimality_certified':False,'invisible_mass_drift':f.invisible_mass_drift})
print(a.output)
