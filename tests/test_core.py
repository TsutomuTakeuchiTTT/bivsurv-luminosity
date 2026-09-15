from fractions import Fraction as F
from pathlib import Path
import sys
import numpy as np
import pytest
from bivsurv.model import GridModel
from bivsurv.em import quantities,fit_em

@pytest.fixture
def model():
    return GridModel((np.arange(-2,3,dtype=float),)*2,np.array([[1.,-1.],[-1.,1.],[1.,1.]]))

def test_all_records_partition(model):
    for A,V in zip(model.dense_operators(),model.visible):
        np.testing.assert_array_equal(A.sum(axis=0),V)

def test_threshold_equality_is_detection(model):
    z=np.array([[1.,-1.],[np.nextafter(1.,np.inf),-1.],[1.,np.nextafter(-1.,np.inf)],[3.,3.]])
    inc,delta,lab,idx=model.record_catalogue(z,np.zeros(4,int))
    assert delta.tolist()==[[True,True],[False,True],[True,False],[False,False]]
    assert inc.tolist()==[True,True,True,False]

def test_public_record_does_not_leak_hidden_value(model):
    z=np.array([[3.,-1.],[-1.,3.],[3.,3.],[0.,0.]])
    fields=np.array([0,1,2,2])
    base=model.record_catalogue(z,fields)
    changed=z.copy(); changed[~base[1]]+=100
    after=model.record_catalogue(changed,fields)
    for a,b in zip(base,after):np.testing.assert_array_equal(a,b)

def test_reject_unaligned_threshold():
    with pytest.raises(ValueError):GridModel((np.array([0.,1.]),)*2,np.array([[.5,1.]]))

def test_same_likelihood_invisible_family(model):
    m=np.arange(1,model.H+1,dtype=float);m/=m.sum()
    Q=m.copy();Q[model.invisible]=0;Q/=Q.sum()
    R=np.zeros(model.H);R[model.invisible]=1/model.invisible.sum()
    for t in [0.,.1,.4]:
        a,q=model.probabilities((1-t)*Q+t*R)
        a0,q0=model.probabilities(m)
        for aa,qq,bb,rr in zip(a,q,a0,q0):np.testing.assert_allclose(aa/qq,bb/rr,atol=1e-14)

def test_old_fraction_update_and_atomic_geometry():
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'validation_reference'))
    import integrated as old
    for limits in [[[1,-1],[-1,1],[1,1]], [[0,0]], [[2,1],[1,2]]]:
        edges=[[-2,-1,0,1,2]]*2
        legacy=old.build_model(edges,limits,[])
        new=GridModel(tuple(np.asarray(e,float) for e in edges),np.asarray(limits,float))
        reps=legacy['representatives'];fields=np.repeat(np.arange(len(limits)),len(reps))
        zz=np.tile(np.array(reps,float),(len(limits),1))
        inc,delta,labels,idx=new.record_catalogue(zz,fields)
        for j in range(len(limits)):
            for h,z in enumerate(reps):
                label=old.record(z,legacy['limits'][j],legacy['bins'])
                assert label==(tuple(labels[j*len(reps)+h]) if inc[j*len(reps)+h] else None)
        # Positive atomic/open-cell masses test all boundary cases, not only continuous laws.
        fine=tuple(F(h+1) for h in range(len(reps)));total=sum(fine);fine=tuple(x/total for x in fine)
        mapping=np.ravel_multi_index(tuple(np.searchsorted(new.z_edges[b],np.array(reps,float)[:,b],side='left') for b in (0,1)),new.shape)
        mass=np.bincount(mapping,weights=np.array(fine,float),minlength=new.H)
        counts=[np.array([1+(k%3) for k in range(len(g['records']))]) for g in legacy['groups']]
        for _ in range(3):
            fine,nu,K=old.em_step(legacy['groups'],[n.tolist() for n in counts],fine)
            mass,_,_,_,Knew=quantities(new,counts,mass)
            expected=np.bincount(mapping,weights=np.array(fine,float),minlength=new.H)
            np.testing.assert_allclose(mass,expected,atol=3e-15,rtol=3e-14)
            assert abs(Knew-float(K))<1e-9

def test_monotonicity_and_invisible_mass(model):
    init=np.ones(model.H)/model.H
    counts=[np.arange(1,len(r)+1) for r in model.labels]
    fit=fit_em(model,counts,init,max_iterations=2000,score_tolerance=1e-9)
    assert fit.min_loglik_increment>=-1e-11
    assert fit.invisible_mass_drift<1e-12
    assert abs(fit.mass.sum()-1)<1e-12

def test_gradient_finite_difference(model):
    rng=np.random.default_rng(23);m=rng.dirichlet(3*np.ones(model.H))
    counts=[np.ones(len(r)) for r in model.labels]
    _,ell,g,_,_=quantities(model,counts,m)
    v=rng.normal(size=model.H);v-=v.mean();h=1e-7
    _,ep,_,_,_=quantities(model,counts,m+h*v)
    _,em,_,_,_=quantities(model,counts,m-h*v)
    assert abs((ep-em)/(2*h)-g@v)<2e-6
    assert abs(g@m)<1e-12

def test_row_and_coordinate_invariance(model):
    init=np.arange(1,model.H+1,dtype=float);init/=init.sum()
    counts=[np.arange(1,len(r)+1) for r in model.labels]
    m1,*_=quantities(model,counts,init)
    swapped=GridModel(model.z_edges[::-1],model.limits_z[:,::-1])
    cc=[]
    for j in range(model.J):
        cc.append(np.array([counts[j][model._lookup[j][(r[1],r[0],r[3],r[2])]] for r in swapped.labels[j]]))
    m2,*_=quantities(swapped,cc,init.reshape(model.shape).T.ravel())
    np.testing.assert_allclose(m1.reshape(model.shape).T.ravel(),m2,atol=1e-14)

def test_readback_model(model):
    again=GridModel.from_dict(model.as_dict())
    np.testing.assert_array_equal(model.visible,again.visible)

def test_bin_partition_noninformative_density(model):
    # Two completely invisible cells can hold arbitrarily different masses.
    ids=np.flatnonzero(model.invisible)
    assert len(ids)>1
    m=np.ones(model.H)/model.H;n=m.copy();s=n[ids].sum();n[ids]=0;n[ids[0]]=s
    a,q=model.probabilities(m);b,r=model.probabilities(n)
    for aa,bb in zip(a,b):np.testing.assert_allclose(aa,bb,atol=1e-14)
    np.testing.assert_allclose(q,r,atol=1e-14)
