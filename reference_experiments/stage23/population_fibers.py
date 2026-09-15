"""Population identification ranges, NOT confidence intervals.

Exact true record probabilities are evaluator inputs. LP proposes extrema;
primal masses are reconstructed rationally and every dual column is checked.
Analytic endpoint formulas for the aligned grid are independently compared.
"""
from __future__ import annotations
from depth_study import *
import scipy
from scipy.optimize import linprog
from sympy import Matrix, Rational


def rows_for(model,classes,groups,trueprobs):
    H=len(classes)
    E=[tuple(F(1) for _ in range(H))]
    for g,pp in zip(groups,trueprobs):
        for rr,p in zip(g['records'],pp):
            E.append(tuple(F(a)-p*v for a,v in zip(rr['A'],g['V'])))
    return E,[F(1)]+[F(0)]*(len(E)-1),[tuple(map(F,g['V'])) for g in groups]


def exact_primal(res,E,b,V):
    supp=[i for i,x in enumerate(res.x) if x>1e-9]
    active=[j for j,v in enumerate(V) if abs(np.dot(np.array(v,float),res.x)-.5)<1e-7]
    rows=E+[V[j] for j in active];rhs=b+[F(1,2)]*len(active)
    def sym(x):return Rational(x.numerator,x.denominator)
    M=Matrix([[sym(r[h]) for h in supp] for r in rows]);d=Matrix(list(map(sym,rhs)))
    sol,parameters=M.gauss_jordan_solve(d)
    if len(parameters):
        # A numerical LP extreme point should have a uniquely determined support.
        raise ArithmeticError('Degenerate rational reconstruction requires explicit handling')
    mass=[F(0)]*len(V[0])
    for h,x in zip(supp,sol):mass[h]=F(int(x.p),int(x.q))
    assert min(mass)>=0 and all(dot(r,mass)==x for r,x in zip(E,b))
    assert all(dot(v,mass)>=F(1,2) for v in V)
    return tuple(mass)


def solve_bound(c,E,b,V,sign):
    # minimize sign*c, with -V*m <= -1/2.
    obj=np.array(c,float)*sign
    res=linprog(obj,A_eq=np.array(E,float),b_eq=np.array(b,float),
                A_ub=-np.array(V,float),b_ub=-.5*np.ones(len(V)),bounds=(0,None),method='highs')
    if not res.success:raise ArithmeticError(res.message)
    mass=exact_primal(res,E,b,V)
    eq=[F(float(x)).limit_denominator(10**12) for x in res.eqlin.marginals]
    iq=[min(F(0),F(float(x)).limit_denominator(10**12)) for x in res.ineqlin.marginals]
    # Exact correction along the normalization equality, not a constraint relaxation.
    residual=[F(sign)*c[h]-sum((eq[k]*E[k][h] for k in range(len(E))),F(0))
              +sum((iq[j]*V[j][h] for j in range(len(V))),F(0)) for h in range(len(c))]
    eq[0]+=min(residual)
    columns=[sum((eq[k]*E[k][h] for k in range(len(E))),F(0))
             -sum((iq[j]*V[j][h] for j in range(len(V))),F(0)) for h in range(len(c))]
    assert all(x<=sign*y for x,y in zip(columns,c)) and max(iq)<=0
    lower=dot(eq,b)-sum(iq,F(0))/2
    objinner=sign*dot(c,mass)
    assert lower<=objinner and objinner-lower<F(1,10**8)
    return {'sign':sign,'mass':mass,'equality_multipliers':eq,'inequality_multipliers':iq,
            'objective_lower_bound':lower,'objective_at_witness':objinner,
            'survival_at_witness':dot(c,mass),'gap':objinner-lower}


def run_population(out):
    rows=[];maximum_gap=F(0);analytic_count=0;witness_count=0
    for resolution in ('coarse','aligned'):
        for depth in (2,3,4):
            spec,model,obs,oc,initial,coeff,inv=setup(resolution,depth)
            groups,classes=quotient(model,model['targets_coefficients'])
            for li,lam in enumerate(LAMBDAS):
                law=FGMLaw(lam);truth=aggregate_masses(measure_masses(model,law),classes)
                probs=conditional_probabilities(groups,truth)
                E,b,V=rows_for(model,classes,groups,probs)
                d=law.survival((F(depth),F(depth)));qmin=min(dot(v,truth) for v in V)
                theta_max=1-(1-d)/(2*qmin)
                assert 0<=theta_max<F(1)
                invj=[h for h in range(len(truth)) if all(not g['V'][h] for g in groups)]
                assert len(invj)==1
                Q=tuple(F(0) if h==invj[0] else x/(1-d) for h,x in enumerate(truth))
                Pupper=tuple((1-theta_max)*x+(theta_max if h==invj[0] else 0) for h,x in enumerate(Q))
                for ti,t in enumerate(model['targets']):
                    c=tuple(F(model['targets_coefficients'][ti][C[0]]) for C in classes)
                    low=solve_bound(c,E,b,V,1);up=solve_bound(c,E,b,V,-1)
                    maximum_gap=max(maximum_gap,low['gap'],up['gap'])
                    analytic=None
                    if resolution=='aligned':
                        # The deepest record determines the target indicator and all V_j.
                        for rec in groups[2]['records']:
                            inside=[h for h,a in enumerate(rec['A']) if a]
                            assert len({c[h] for h in inside})==1
                            for v in V:assert len({v[h] for h in inside})==1
                        lower=(law.survival(t)-d)/(1-d)
                        upper=1-(1-law.survival(t))/(2*qmin)
                        assert low['objective_lower_bound']<=lower<=low['survival_at_witness']
                        assert up['survival_at_witness']<=upper<=-up['objective_lower_bound']
                        assert dot(c,Q)==lower and dot(c,Pupper)==upper
                        for mm in (Q,Pupper):
                            assert conditional_probabilities(groups,mm)==probs
                            assert all(dot(v,mm)>=F(1,2) for v in V)
                            witness_count+=1
                        analytic={'lower':lower,'upper':upper,'theta_max':theta_max,
                                  'lower_mass':Q,'upper_mass':Pupper}
                        analytic_count+=1
                    cert={'resolution':resolution,'depth':depth,'lambda':lam,'target_index':ti,
                          'target':t,'true_probabilities':probs,'lower':low,'upper':up,'analytic':analytic}
                    save(out/'population_certificates'/f'{resolution}_d{depth}_lambda{li}_t{ti}.json',cert)
                    rows.append({'resolution':resolution,'depth':depth,'lambda':str(lam),'target':str(tuple(t)),
                       'target_index':ti,'true_survival':str(law.survival(t)),
                       'lower_out':str(low['objective_lower_bound']),'lower_in':str(low['survival_at_witness']),
                       'upper_in':str(up['survival_at_witness']),'upper_out':str(-up['objective_lower_bound']),
                       'lower_display':float(low['survival_at_witness']),
                       'upper_display':float(up['survival_at_witness']),
                       'gap_max':float(max(low['gap'],up['gap'])),
                       'analytic_lower':str(analytic['lower']) if analytic else '',
                       'analytic_upper':str(analytic['upper']) if analytic else ''})
            print('population',resolution,depth,flush=True)
    save_csv(out/'population_ranges.csv',rows)
    result={'population_ranges':len(rows),'rational_lp_witnesses':2*len(rows),
            'analytic_aligned_ranges':analytic_count,'analytic_endpoint_witnesses':witness_count,
            'maximum_endpoint_gap':maximum_gap,'all_models_include_full_real_plane':True,
            'sample_zeros_used_as_population_zeros':False,
            'scipy_version':scipy.__version__}
    save(out/'population_checks.json',result)
    return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs');args=p.parse_args()
    print(json.dumps(encode(run_population(args.outdir)),indent=2))
