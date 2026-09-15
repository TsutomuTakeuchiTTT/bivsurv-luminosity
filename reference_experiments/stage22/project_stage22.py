"""Three prespecified diagnostic projections, not a coverage simulation.
Uses the Stage20 common-parent engine unchanged; every certificate is checked
by the separately written Stage20 verifier. Resource-limited brackets remain.
"""
from pathlib import Path
import sys,json,csv,time,argparse
from fractions import Fraction as F
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'vendor20'));sys.path.insert(0,str(BASE/'vendor20/vendor19'))
from integrated import compile_input
from exact import encode
from column_certificates import project_columns
from verify_stage20 import verify

def run(outdir):
    summaries=[]
    for li in range(3):
        case=f'lambda{li}_N400_r000'
        spec=json.loads((outdir/'observed_inputs'/f'{case}.json').read_text())
        # The representative is chosen in advance: first replicate, smaller n,
        # one per parent law, same target (1,1). No selection for easy solutions.
        data,*_=compile_input(spec,0)
        doc=project_columns(data,tolerance=F(1,1000),max_nodes=15,proposal_maxiter=100)
        directory=outdir/'projection_certificates';directory.mkdir(exist_ok=True)
        p=directory/f'{case}.json';p.write_text(json.dumps(encode(doc),sort_keys=True,indent=2)+'\n')
        checked=verify(json.loads(p.read_text()))
        ev=json.loads((outdir/'fits'/f'{case}.json').read_text())['evaluation']
        truth=F(str(ev['Strue_0']))
        # Truth values used only AFTER the certificate is completed.
        included=F(doc['lower_bracket'][0])<=truth<=F(doc['upper_bracket'][1])
        row={'case':case,'nodes':len(doc['nodes']),'precision_reached':doc['precision_reached'],
             'lower_outer':str(doc['lower_bracket'][0]),'lower_inner':str(doc['lower_bracket'][1]),
             'upper_inner':str(doc['upper_bracket'][0]),'upper_outer':str(doc['upper_bracket'][1]),
             'true_survival':str(truth),'outer_covers_truth':included,'seconds':doc['seconds'],
             'independent_verification':checked}
        summaries.append(row);print(case,'nodes',len(doc['nodes']),'precision',doc['precision_reached'],
            'lower',[float(x) for x in doc['lower_bracket']],
            'upper',[float(x) for x in doc['upper_bracket']],'truth covered',included,flush=True)
    (outdir/'projection_results.json').write_text(json.dumps(encode(summaries),sort_keys=True,indent=2)+'\n')
    return summaries
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs');a=p.parse_args();run(a.outdir)
