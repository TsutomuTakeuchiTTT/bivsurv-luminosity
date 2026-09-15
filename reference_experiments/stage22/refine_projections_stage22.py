"""Stage21 multiple-plane refinement of three prespecified Stage22 projections."""
from pathlib import Path
import sys,json
from fractions import Fraction as F
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'vendor21'))
from refine import Refinement
from verify_stage21 import verify
from exact import encode

def run(out):
    results=[]
    for li in range(3):
        case=f'lambda{li}_N400_r000'
        old=json.loads((out/'projection_certificates'/f'{case}.json').read_text())
        doc=Refinement(old,mode='relative').run(tolerance=F(1,1000),max_nodes=65,maxiter=150)
        p=out/'refined_projection_certificates'/f'{case}.json';p.parent.mkdir(exist_ok=True)
        p.write_text(json.dumps(encode(doc),sort_keys=True,indent=2)+'\n')
        vv=verify(json.loads(p.read_text()),old,True)
        results.append({'case':case,'verification':vv,'lower_bracket':doc['lower_bracket'],
                        'upper_bracket':doc['upper_bracket'],'seconds':doc['seconds']})
        print(case,'NEW',doc['precision_reached'],len(doc['nodes']),
            [float(x) for x in doc['lower_bracket']],[float(x) for x in doc['upper_bracket']],flush=True)
    (out/'refined_projection_results.json').write_text(json.dumps(encode(results),sort_keys=True,indent=2)+'\n')
if __name__=='__main__':run(BASE/'outputs')
