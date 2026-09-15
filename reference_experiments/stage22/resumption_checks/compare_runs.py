"""Compare preserved Stage22 outputs with the complete resumption rerun."""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import json, hashlib, re, csv, zipfile
import numpy as np
BASE=Path(__file__).resolve().parents[1]
a=BASE/'outputs'; b=BASE/'resumption_checks/sampling_rerun'
def without_times(x):
    if isinstance(x,dict): return {k:without_times(v) for k,v in x.items() if k!='seconds'}
    if isinstance(x,list): return [without_times(v) for v in x]
    return x
checks={}
for folder in ('observed_inputs','fits'):
    aa=sorted((a/folder).glob('*.json'))
    assert len(aa)==120
    checks[folder+'_byte_identical']=all(p.read_bytes()==(b/folder/p.name).read_bytes() for p in aa)
for name in ('trials.csv','summary.csv'):
    checks[name+'_byte_identical']=(a/name).read_bytes()==(b/name).read_bytes()
for folder in ('projection_certificates','refined_projection_certificates'):
    pp=sorted((a/folder).glob('*.json'))
    assert len(pp)==3
    checks[folder+'_identical_excluding_seconds']=all(without_times(json.loads(p.read_text()))==without_times(json.loads((b/folder/p.name).read_text())) for p in pp)
for name in ('stage22_sampling_results.json','projection_results.json','refined_projection_results.json','independent_verification.json'):
    checks[name+'_identical_excluding_seconds']=without_times(json.loads((a/name).read_text()))==without_times(json.loads((b/name).read_text()))
array_count=0
for p in sorted((a/'latent_evaluation_only').glob('*.npz')):
    with np.load(p) as old,np.load(b/'latent_evaluation_only'/p.name) as new:
        assert sorted(old.files)==sorted(new.files)
        for k in old.files:
            assert np.array_equal(old[k],new[k]),(p,k)
            array_count+=1
checks['600_latent_arrays_identical']=array_count==600
# Original Stage20/21 routines really are byte-identical to earlier archives.
for stage,files in [(20,['column_certificates.py','verify_stage20.py']),(21,['refine.py','verify_stage21.py'])]:
    with zipfile.ZipFile(f'/mnt/data/bivariate_survival_reconstruction_stage{stage}.zip') as zz:
        for name in files:
            candidates=[x for x in zz.namelist() if x.endswith('/'+name)]
            direct=[x for x in candidates if x==f'bivariate_survival_reconstruction_stage{stage}/'+name]
            if len(direct)!=1:raise RuntimeError((stage,name,candidates))
            checks[f'unchanged_stage{stage}_{name}']=zz.read(direct[0])==(BASE/f'vendor{stage}'/name).read_bytes()
assert all(checks.values()),checks
report={'all_comparisons_passed':True,
  'scope':'Resumption rerun: 120 parent samples and fits; 3 Stage20 initial projections; 3 Stage21 refined projections; independent verification of all saved and rerun results.',
  'observed_input_json_count':120,'fit_json_count':120,'latent_array_count':array_count,
  'initial_projection_rerun_count':3,'refined_projection_rerun_count':3,
  'all_120_sharp_projections_computed':False,
  'checks':checks}
(BASE/'resumption_checks/resumption_reproducibility.json').write_text(json.dumps(report,indent=2)+'\n')
# LaTeX structural lint only (not a compilation claim).
s=Path('/mnt/data/bivariate_survival_stage21_note.tex').read_text()
sn='\n'.join(line.split('%')[0] for line in s.splitlines())
labels=re.findall(r'\\label\{([^}]+)\}',sn)
refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',sn)
stack=[]
for m in re.finditer(r'\\(begin|end)\{([^}]+)\}',sn):
    what,env=m.groups()
    if what=='begin':stack.append(env)
    else:
        assert stack and stack.pop()==env,(m.start(),env)
assert not stack
lint={'duplicate_labels':[k for k,v in Counter(labels).items() if v>1],
      'undefined_internal_references':sorted(set(refs)-set(labels)),
      'environment_nesting_ok':True,'compiled':False,
      'new_citation_keys':re.findall(r'\\cite\w*\{([^}]+)\}',sn)}
assert not lint['duplicate_labels'] and not lint['undefined_internal_references']
(BASE/'resumption_checks/stage21_latex_lint.json').write_text(json.dumps(lint,indent=2)+'\n')
# Portable projection table: keep exact rational endpoints and explicitly outward decimals.
def floor8(q):q=Fraction(q);return f'{(q.numerator*10**8)//q.denominator/10**8:.8f}'
def ceil8(q):q=Fraction(q);return f'{-((-q.numerator*10**8)//q.denominator)/10**8:.8f}'
runs=json.loads((b/'refined_projection_results.json').read_text())
rows=[]
for lam,d in zip(('-3/4','0','3/4'),runs):
    lo=list(map(Fraction,d['lower_bracket']));hi=list(map(Fraction,d['upper_bracket']))
    row={'lambda':lam,'case':d['case'],'minimum_lower':str(lo[0]),'minimum_upper':str(lo[1]),
         'maximum_lower':str(hi[0]),'maximum_upper':str(hi[1]),
         'minimum_outward_8dp':f'[{floor8(lo[0])}, {ceil8(lo[1])}]',
         'maximum_outward_8dp':f'[{floor8(hi[0])}, {ceil8(hi[1])}]',
         'minimum_gap':str(lo[1]-lo[0]),'maximum_gap':str(hi[1]-hi[0]),
         'nodes':d['verification']['nodes'],'precision_reached':d['verification']['precision_reached']}
    rows.append(row)
with (BASE/'resumption_checks/projection_endpoint_table.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print(json.dumps(report,indent=2))
for r in rows:print(r['lambda'],r['minimum_outward_8dp'],r['maximum_outward_8dp'])
print('latex',lint)
