"""Independent arithmetic verification; no LP, optimizer or data generation.
Reconstructs observation matrices from public interval records. The M-step is
verified via its KKT equations, NOT by calling the proposal/active-set routine.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,copy
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'vendor'))
from coarsened_model import build_model,quotient
from exact import logint

def dot(a,b):return sum((F(x)*F(y) for x,y in zip(a,b)),F(0))
def serial(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [serial(a) for a in x]
 return x

def reference_mass(cell):
    # Independent integration of the three reference-density pieces.
    ans=F(1)
    for I in cell:
        if I.lo is not None and I.lo==I.hi:return F(0)
        mass=F(0);a=I.lo;b=min(I.hi,F(1)) if I.hi is not None else F(1)
        if a is None or a<b:mass+=1/(3*(2-b))-(1/(3*(2-a)) if a is not None else 0)
        a=max(I.lo,F(1)) if I.lo is not None else F(1);b=min(I.hi,F(2)) if I.hi is not None else F(2)
        if a<b:mass+=(b-a)/3
        a=max(I.lo,F(2)) if I.lo is not None else F(2);b=I.hi
        if b is None or a<b:mass+=1/(3*(a-1))-(1/(3*(b-1)) if b is not None else 0)
        ans*=mass
    return ans

def data_from_public(inp):
 mod=build_model(inp['edges'],inp['limits'],inp['targets']);gs,C=quotient(mod)
 ini=tuple(sum((reference_mass(mod['cells'][h]) for h in cc),F(0)) for cc in C)
 A=[];VR=[];ns=[];counts=[[x+y for x,y in zip(a,b)] for a,b in zip(inp['train'],inp['validation'])]
 for g,row in zip(gs,counts):
  for r,n in zip(g['records'],row):
   if n:A.append(r['A']);VR.append(g['V']);ns.append(n)
 V=[g['V'] for g in gs];D=[[1-x for x in v] for v in V];sig=[2*D[0][h]+D[1][h] for h in range(len(ini))]
 return {'A':A,'VR':VR,'counts':ns,'V':V,'D':D,'sig':sig,'initial':ini,'N':sum(ns),'groups':gs}

def admissible(d,m):
 assert min(m)>=0 and sum(m)==1
 assert all(dot(v,m)>=F(1,2) for v in d['V'])
 assert all(dot(a,m)>0 for a in d['A'])

def estep_independent(d,m,counts=None):
 H=len(m);nu=[F(0)]*H;K=F(0)
 if counts is None:
  rows=list(zip(d['A'],d['VR'],d['counts']))
 else:
  rows=[]
  for g,cc in zip(d['groups'],counts):
   for rr,n in zip(g['records'],cc):
    if n:rows.append((rr['A'],g['V'],n))
 for A,V,n in rows:
  a=dot(A,m);q=dot(V,m);assert a>0 and q>0;K+=n/q
  for h in range(H):
   if A[h]:nu[h]+=n*m[h]/a
   if not V[h]:nu[h]+=n*m[h]/q
 assert sum(nu)==K
 return tuple(v/K for v in nu),K

def likelihood_gain(d,m,u):
 lo=hi=F(0)
 for A,V,n in zip(d['A'],d['VR'],d['counts']):
  a=dot(A,m);ap=dot(A,u);q=dot(V,m);qp=dot(V,u)
  l,h=logint(ap*q/(a*qp));lo+=n*l;hi+=n*h
 return lo/d['N'],hi/d['N']

def check_step(d,m,step):
 u=tuple(map(F,step['mass']));admissible(d,u)
 w,K=estep_independent(d,m)
 eta=F(step['kkt_eta']);beta=tuple(map(F,step['kkt_tail_multipliers']));assert min(beta)>=0
 den=[eta+beta[0]*d['D'][0][h]+beta[1]*d['D'][1][h] for h in range(len(m))]
 assert min(den)>0
 star=tuple(v/t for v,t in zip(w,den));admissible(d,star)
 assert all(beta[j]*(dot(d['D'][j],star)-F(1,2))==0 for j in range(2))
 # All 16 coordinate KKT identities certify a GLOBAL surrogate maximum.
 assert all(w[h]/star[h]==den[h] for h in range(len(m)) if w[h])
 sl=su=F(0)
 for a,b,c in zip(w,m,u):
  if a:
   l,h=logint(c/b);sl+=a*l;su+=a*h
 ll,lh=likelihood_gain(d,m,u)
 assert [sl,su]==list(map(F,step['normalized_surrogate_gain']))
 assert [ll,lh]==list(map(F,step['mean_likelihood_gain']))
 assert sl>=0 and ll>=K*su/d['N']>=0
 gap=sum((v*(a/b-1) for v,a,b in zip(w,star,u) if v),F(0))
 assert gap==F(step['exact_Mstep_surrogate_loss_upper']) and gap>=0
 return u,gap,ll

def gradient(d,m):
 H=len(m);out=[F(0)]*H
 for A,V,n in zip(d['A'],d['VR'],d['counts']):
  a=dot(A,m);q=dot(V,m)
  for h in range(H):out[h]+=F(n,d['N'])*(F(A[h])/a-F(V[h])/q)
 assert dot(out,m)==0
 return out

def check_stationarity(d,m,c):
 admissible(d,m);g=gradient(d,m);beta=list(map(F,c['multipliers']));z=F(c['simplex_multiplier'])
 assert len(beta)==2 and min(beta)>=0
 for h in range(len(m)):assert beta[0]*d['D'][0][h]+beta[1]*d['D'][1][h]+z>=g[h]
 gap=z+sum(beta,F(0))/2-dot(g,m)
 assert gap==F(c['directional_gap_upper']) and gap>=0
 assert F(c['maximum_raw_score'])==max(g)
 return gap

def run(out):
 # Hard failure if a forbidden optimizer happens to be invoked by imported code.
 import scipy.optimize
 def forbidden(*a,**k):raise AssertionError('Optimizer is disabled during certificate verification')
 scipy.optimize.minimize=scipy.optimize.linprog=scipy.optimize.root=scipy.optimize.brentq=forbidden
 names=json.loads((BASE/'inputs/selection_manifest.json').read_text())['names'];updates=0;col=0;maxgap=F(0);maxstat=F(0);minmean=None
 first=None;prediction_checks=0;cases=[]
 for name in names:
  p=BASE/'inputs/observed'/f'{name}.json';inp=json.loads(p.read_text());oldp=BASE/'inputs/prior_fits'/f'{name}.json';old=json.loads(oldp.read_text())
  cert=json.loads((out/'certificates'/f'{name}.json').read_text());assert cert['case']==name
  assert cert['input_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
  assert cert['old_fit_sha256']==hashlib.sha256(oldp.read_bytes()).hexdigest()
  assert old['evaluation']['tail_feasible'] is False
  d=data_from_public(inp);assert serial({k:d[k] for k in cert['compiled']})==cert['compiled']
  assert [str(x) for x in d['initial']]==old['initial_mass']
  pred=d['initial']
  for _ in range(2):pred,K=estep_independent(d,pred,inp['train'])
  assert [str(x) for x in pred]==old['prediction'];prediction_checks+=1
  m=d['initial'];admissible(d,m)
  for k,step in enumerate(cert['GEM_steps']):
   assert step['iteration']==k+1
   if first is None:first=(copy.deepcopy(d),m,copy.deepcopy(step))
   m,gap,ll=check_step(d,m,step);updates+=1;col+=len(m)
   maxgap=max(maxgap,gap);minmean=ll if minmean is None else min(minmean,ll)
  gemgap=check_stationarity(d,m,cert['GEM_final_stationarity'])
  q=cert['independent_candidate'];u=tuple(map(F,q['mass']));g=check_stationarity(d,u,q['stationarity'])
  assert g<F(1,10**8) and q['global_maximum_certified'] is False and q['is_GEM_step'] is False
  ga=likelihood_gain(d,m,u);assert ga[0]>=0 and list(ga)==list(map(F,q['ascent_from_GEM']))
  maxstat=max(maxstat,g);cases.append({'case':name,'candidate_stationarity_upper':str(g),'GEM_stationarity_upper':str(gemgap)})
 # Cheap focused tamper checks on a first-step certificate.
 d,m,st=first;rejected=0
 for typ in ('mass','eta','negative_multiplier','gain','loss'):
  t=copy.deepcopy(st)
  if typ=='mass':t['mass'][0]=str(F(t['mass'][0])+F(1,10))
  if typ=='eta':t['kkt_eta']=str(F(t['kkt_eta'])+1)
  if typ=='negative_multiplier':t['kkt_tail_multipliers'][0]='-1'
  if typ=='gain':t['mean_likelihood_gain'][0]=str(F(t['mean_likelihood_gain'][0])+1)
  if typ=='loss':t['exact_Mstep_surrogate_loss_upper']='-1'
  try:check_step(d,m,t)
  except (AssertionError,ValueError,ZeroDivisionError):rejected+=1
  else:raise AssertionError('Harmful tamper accepted: '+typ)
 result={'passed':True,'case_count':len(names),'public_model_reconstructions':len(names),'unchanged_training_predictors':prediction_checks,
   'certified_GEM_updates':updates,'surrogate_KKT_coordinate_identities':col,
   'maximum_normalized_surrogate_loss_upper':str(maxgap),'maximum_normalized_surrogate_loss_upper_float':float(maxgap),
   'minimum_mean_likelihood_gain_lower':str(minmean),'minimum_mean_likelihood_gain_lower_float':float(minmean),
   'independent_stationary_candidates':len(cases),'maximum_directional_gap_upper_float':float(maxstat),
   'first_order_only_not_global_likelihood_gap':True,'harmful_tampers_rejected':rejected,
   'LP_called':False,'nonlinear_optimizer_called':False,'sampler_called':False,'Mstep_active_set_solver_called':False}
 (out/'independent_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)
 return result

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=BASE/'outputs');a=p.parse_args();run(a.outdir)
