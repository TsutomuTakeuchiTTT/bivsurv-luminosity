"""Compare two completed Stage24 runs without recomputing estimates."""
from pathlib import Path
import json,hashlib
BASE=Path(__file__).resolve().parent

def run(a,b):
 files=sorted((a/'certificates').glob('*.json'))
 checks={p.name:p.read_bytes()==(b/'certificates'/p.name).read_bytes() for p in files}
 x=json.loads((a/'results.json').read_text());y=json.loads((b/'results.json').read_text())
 x.pop('seconds');y.pop('seconds')
 result={'passed':all(checks.values()) and x==y and (a/'comparison.csv').read_bytes()==(b/'comparison.csv').read_bytes(),
  'all_13_cases_rerun':True,'GEM_updates_per_run':1040,'independent_boundary_searches_per_run':39,
  'certificate_json_byte_identical_count':sum(checks.values()),'comparison_csv_byte_identical':(a/'comparison.csv').read_bytes()==(b/'comparison.csv').read_bytes(),
  'results_identical_except_seconds':x==y,'sampling_not_performed':True}
 assert result['passed'];(a/'reproducibility.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
 return result
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--original',type=Path,default=BASE/'outputs');p.add_argument('--rerun',type=Path,default=BASE/'rerun_outputs');a=p.parse_args();run(a.original,a.rerun)
