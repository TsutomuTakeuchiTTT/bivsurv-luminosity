"""The manuscript's D-completion update, evaluated with indexed sums.

This routine is UNCONSTRAINED on the probability simplex. A supplied tail
bound is diagnosed, never silently imposed or repaired. See archived
reference code for the tail-constrained M-step and rigorous certificates.
No smoothing, pseudo-counts, floors, latent truth, or rejected counts enter.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .model import GridModel

@dataclass
class FitResult:
    mass: np.ndarray
    history: np.ndarray
    status: str
    iterations: int
    score: float
    q: np.ndarray
    min_loglik_increment: float
    invisible_mass_drift: float
    tail_feasible: bool | None

def quantities(model: GridModel, counts, mass, correction=True):
    a, q = model.probabilities(mass)
    N = sum(float(n.sum()) for n in counts)
    if N <= 0:
        raise ValueError('No included records')
    numerator = np.zeros(model.H)
    missing = np.zeros(model.H)
    penalty = np.zeros(model.H)
    ell, K = 0., 0.
    for j, (n, aj, idx, vis) in enumerate(zip(counts, a, model.record_index, model.visible)):
        nj = float(n.sum())
        if not nj:
            continue
        pos = n > 0
        if np.any(aj[pos] <= 0) or q[j] <= 0:
            raise ArithmeticError('A positive observed record has nonpositive model probability')
        ratio = np.zeros_like(aj)
        ratio[pos] = n[pos]/aj[pos]
        numerator[vis] += ratio[idx[vis]]
        ell += float(n[pos] @ np.log(aj[pos]))
        if correction:
            ell -= nj*np.log(q[j])
            missing[~vis] += nj/q[j]
            penalty[vis] += nj/q[j]
            K += nj/q[j]
        else:
            K += nj
    grad = (numerator-penalty)/N if correction else numerator/N-1
    new_mass = mass*(numerator+missing)/K
    return new_mass, ell/N, grad, q, K

def fit_em(model: GridModel, counts, initial, *, score_tolerance=1e-7, max_iterations=20000, correction=True, tail_upper=None):
    counts = [np.asarray(n, float) for n in counts]
    if len(counts) != model.J or any(n.shape != (len(r),) or np.any(n < 0) for n, r in zip(counts, model.labels)):
        raise ValueError('Invalid counts, including missing possible-record columns')
    m = model.validate_mass(initial).copy()
    inv0 = float(m[model.invisible].sum())
    history = []
    for k in range(max_iterations+1):
        new, ell, grad, q, K = quantities(model, counts, m, correction)
        residual = max(0., float(grad.max()))
        history.append((k, ell, residual, float(m.sum()), float(m[model.invisible].sum())))
        if residual <= score_tolerance:
            status = 'score_tolerance_reached_not_global_certificate'
            break
        if k == max_iterations:
            status = 'iteration_budget_exhausted'
            break
        if np.any(new < 0) or not np.isfinite(new).all() or abs(new.sum()-1) > 2e-10:
            raise ArithmeticError('EM mass invariant failed')
        # Do not renormalize, clip, or add a mass floor.
        m = new
    history = np.asarray(history, float)
    mingain = float(np.min(np.diff(history[:, 1]))) if len(history)>1 else 0.
    if mingain < -1e-11:
        raise ArithmeticError(f'Unexpected observed log-likelihood decrease {mingain}')
    drift = float(np.max(np.abs(history[:,4]-inv0)))
    if correction and drift>1e-9:
        raise ArithmeticError('Completely invisible mass was not preserved')
    feasible = None if tail_upper is None else bool(np.all(q >= 1-np.asarray(tail_upper)-1e-10))
    return FitResult(m,history,status,k,residual,q,mingain,drift,feasible)
