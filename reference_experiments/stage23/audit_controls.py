"""Additional structural controls for depth, recording, and unobservability."""
from depth_study import *
from itertools import product


def run_controls(out):
    results=[];checks=0
    # A complete geometry grid containing every limit, recording, and target cut.
    allcuts=['1','3/2','2','3','4']
    fine=build_model([allcuts,allcuts],[[2,1],[1,2],[2,2],[3,3],[4,4]],TARGETS)
    for res in ('coarse','aligned'):
        e=['1'] if res=='coarse' else allcuts
        md=build_model([e,e],[[2,1],[1,2],[2,2],[3,3],[4,4]],TARGETS)
        for old,new in ((2,3),(3,4)):
            mappings={};examples={}
            for zz in fine['representatives']:
                a=record(zz,(F(old),F(old)),md['bins'])
                b=record(zz,(F(new),F(new)),md['bins'])
                mappings.setdefault(b,set()).add(a);examples.setdefault((b,a),zz);checks+=1
            conflicts=[k for k,v in mappings.items() if len(v)>1]
            if res=='aligned':assert not conflicts
            else:assert conflicts
            ex=[]
            if conflicts:
                k=conflicts[0]
                # deterministic JSON order
                ex=[{'point':examples[(k,a)],'old_record':a,'new_record':k}
                    for a in sorted(mappings[k],key=str)[:2]]
            results.append({'resolution':res,'old_depth':old,'new_depth':new,
                            'old_record_recoverable_from_new':not bool(conflicts),
                            'conflicting_new_records':len(conflicts),'example':ex})
    # Independent absolute ratio integral using only marginal CDF values.
    analytic=[]
    for lam in LAMBDAS:
        for d in (2,3,4):
            tail=(1-F(d,4))**2*(1+lam*F(d,4)**2)
            qmin=F(5,8)-F(3,64)*lam
            theta_max=1-(1-tail)/(2*qmin)
            analytic.append({'lambda':lam,'depth':d,'true_invisible':tail,
                             'maximum_compatible_invisible_mass':theta_max})
    # Confirm the 13 unconstrained estimates outside the common parent class are not roundoff.
    vals=[]
    for file in (out/'fits').glob('*/*.json'):
        fit=json.loads(file.read_text())
        ev=fit['evaluation']
        if not ev['tail_feasible']:
            _,md,gg,*_=setup(ev['resolution'],ev['depth'])
            qs=[float(np.dot(g['V'],fit['final_mass'])) for g in gg]
            vals.append({'design':ev['design'],'case':ev['case'],'selection_probabilities':qs,
                         'max_violation':max(.5-q for q in qs)})
    result={'depth_record_maps':results,'map_cell_comparisons':checks,
            'analytic_tail_controls':analytic,'unconstrained_outside_class':vals,
            'largest_tail_violation':max(x['max_violation'] for x in vals)}
    save(out/'structural_controls.json',result)
    return result
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs');a=p.parse_args()
    print(json.dumps(encode(run_controls(a.outdir)),indent=2))
