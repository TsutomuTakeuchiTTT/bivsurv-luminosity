# Stage 24 results

> English translation of [`reference_experiments/stage24/SUMMARY_ja.md`](../../../reference_experiments/stage24/SUMMARY_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision. For the later global near-optimality certificates of these same 13 candidates, see the separate Stage 25 documentation; this Stage 24 account retains its original scope.

All 13 design/sample cases outside the common Stage 23 tail class were reanalysed. They represent eight distinct parent samples. The same 16-mass observation models, count tables, initial measure R, and tail bounds were retained. Learning predictions and confidence sets were unchanged.

## Updates and stopping are separate questions

1. Eighty constrained GEM updates were executed for each case, giving 1040 updates. Rational feasibility, surrogate increase, and observed-likelihood increase were checked at every update.
2. The closed form for the four groups determined by the two D-membership indicators was used. KKT equalities established global maximization of each M-step. Because the implementation accepts finite-rational approximations, the executed iteration is called GEM.
3. After 80 updates, the certified upper bounds on the largest feasible-direction first-order residual ranged from 5.68e-5 to 5.33e-4. These iterates are not described as converged solutions.
4. Independently, 13 candidates were constructed by maximizing the same constrained observed likelihood from three starting points. Every candidate belonged to the parent class and had higher likelihood than the 80-step iterate. The largest constrained first-order residual bound was 5.85e-14. These candidates are not GEM updates and were not certified at this stage as global MLEs of the original nonconcave likelihood.

## Principal numerical results

| Quantity | Result |
|---|---:|
| Original minimum inclusion probabilities | 0.426358–0.490940 |
| Minimum inclusion probability of constrained candidates | At least 0.5 |
| Rationally checked GEM updates | 1040 |
| Smallest lower bound on mean observed-log-likelihood increase per update | 1.023829e-7 |
| Largest certified upper bound on normalized M-step loss | 2.452849e-24 |
| Boundary candidates constructed separately from GEM | 13 |
| Largest certified first-order residual bound of boundary candidates | 5.847021e-14 |
| Largest unconstrained mean positive score at those same candidates | 0.0328461 |
| Old candidates entering the tail class by removing invisible mass alone | 1/13 |

The 13 original results were valid nonnegative survival distributions; their failure was membership in the additionally specified tail class. Twelve retained a constraint violation even after all completely invisible mass was removed, so changing only the global scale did not resolve the problem. An unconstrained score can be large at a constrained boundary; stopping should be diagnosed through feasible directions.

Reference survival RMSE decreased in all 13 cases, but the analysis selected only prior class-violation cases. It does not establish a general gain in statistical efficiency or coverage.

## Independent verification

All 13 public models were reconstructed, and the original two-step learning EM predictions were confirmed unchanged. The verifier rechecked 1040 updates, 16640 M-step KKT coordinate equalities, and 13 boundary candidates. It invoked neither optimization nor random generation, and did not call the M-step proposal function. Five types of harmful alteration were rejected.

The 256 rational four-group controls included 90 inactive cases, 57 with only the first constraint active, 57 with only the second active, and 52 with both active. Sixteen comparisons with a separate SLSQP calculation gave a maximum objective difference of 6.67e-16.

All 13 cases were rerun in a separate output directory. The 13 certificate JSON files and comparison CSV were byte-identical; aggregate JSON agreed apart from time.

## What was not executed at this stage

There was no new parent-sample generation, refitting of all 720 evaluations, change to the parent confidence sets, new coverage estimation, or proof of global optimality for the 13 original-likelihood candidates. While preparing the Stage 23 record, the earlier 720 inputs and 54 identification certificates were rechecked with an independent verifier; that was not counted as rerunning all earlier searches.
