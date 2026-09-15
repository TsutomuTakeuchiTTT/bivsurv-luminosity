# Stage 22 initial sample validation: continuous parents, non-nested limits, and parent confidence sets

> English translation of [`reference_experiments/stage22/THEORY_ja.md`](../../../reference_experiments/stage22/THEORY_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

## 1. Scope

Through Stage 21, the most recent computational improvements had been tested on deterministic saved count tables. Stage 22 generated individual latent objects from continuous parent distributions, independently assigned limits and training/validation membership, and passed only interval records after A–D selection to the predictor and point estimator. The update formula was unchanged. True-parent membership was examined for 120 parent confidence sets, and certified survival projections were computed for three fixed representative samples. This did not complete sharp projections for all 120 samples or establish general EM consistency.

## 2. Declared observations and evaluation truth

Both coordinates use intervals `(-infinity, 1]` and `(1, infinity)`. The limits are (2,1), (1,2), and (2,2). The first two are non-nested; the complete record alphabets have 5, 5, and 8 categories. The targets are (1,1), (3/2,1), and (2,2). The representation uses 35 cells over the full real plane, nine equivalence classes preserving observations, and 11 classes preserving all targets as well. Boundary mass, negative coordinates, unbounded regions, and completely invisible regions are not removed from the candidates. The external parent-class constraints are `P(D_j) <= 1/2`, with `alpha_sig = 1/20`.

The evaluation-only density is

```text
f_lambda(z1,z2) = 1/16 [1 + lambda(1-z1/2)(1-z2/2)] on (0,4)^2,
lambda = -3/4, 0, 3/4,
```

and is zero elsewhere. Neither the true support nor lambda is passed to the estimation model. This density and the true cell masses are reserved for evaluation. The three true D probabilities are `[3/8+3lambda/64, 3/8+3lambda/64, 1/4+lambda/16]`; every truth belongs to the declared tail class.

Generate independent uniform variables U and W, set `a = lambda(1-2U)`, and use

```text
V = 2W / [1+a+sqrt((1+a)^2-4aW)],
Z = (4U,4V).
```

This solves `F(V|U) = V+a V(1-V) = W`. Random generation uses double precision. The maximum conditional-CDF identity residual for the stored Z values was 3.33067e-16. Probability evaluation uses rational values of rectangular density integrals.

Parent sample sizes are 400 and 1600, with 20 replicates per lambda and sample size: 120 samples and 120000 parent objects in total. Three generators branch from `SeedSequence([20260913,lambda_index,N_parent,rep])`. The limit is chosen uniformly from the three possibilities, and the training flag has probability 1/2; both assignments are independent of Z. Splitting precedes selection, and observed counts are collected after removing D. Catalogue size, per-limit counts, and training/validation counts are therefore random. This is a fixed-parent-size experiment, not a multinomial experiment with fixed catalogue size.

## 3. Separating learning prediction, point estimation, and truth evaluation

Learning predictions start with uniform mass over the nine observation classes and apply two rational steps of the existing EM to training records only. All-data point estimates apply the same original EM to the sum of training and validation counts, using a mean positive-score threshold of 1e-8 and at most 4000 updates. The all-data fit is not reused as the validation predictor. Neither procedure receives truth values or D objects.

The original EM is not a tail-constrained optimization. Tail-class membership of final candidates is diagnosed separately. All 120 final candidates satisfied the tail condition in this experiment, but general preservation of that condition is not claimed. A KKT-type stopping condition is not a proof of global optimality.

Reference survival values use a predeclared allocation that preserves each observation-class total and divides it equally among the target-preserving subclasses. The internal allocation at (3/2,1) is not identified by the estimation data. RMSE at the three fixed targets is saved separately from record-law TV weighted by per-layer catalogue counts. An oracle empirical distribution using the selected objects' hidden latent values and an empirical distribution of the full parent sample, including D, are stored separately for evaluation only.

## 4. True-parent membership in the confidence set

Use the set `C = {P: L_val(P) >= alpha_sig L_val(Pred)}` within the declared parent class, with the original validation likelihood. For the true P0, enclose `log E0 = sum n_jr log[g_jr/pi0_jr]` in a rational interval and compare it with `log(1/alpha_sig)`. Calibration retains independent training and validation and normalization over all possible records. Membership of the true distribution in the candidate model is checked separately through continuous-density cell integrals and tail constraints.

The true parent belonged to the confidence set in all 120 samples. This is an observed outcome for six settings with 20 replicates each, not a theoretical coverage probability of 1. Calibration uses the previously adopted split-likelihood theory; it is not reinvented in this experiment. The result is also not an exact enumeration of coverage over all possible count tables.

## 5. Observational fit is different from parent-distribution recovery

The region `I = (2, infinity)^2` is completely invisible at every limit. The original EM preserves its initial class mass 1/9, so `S_hat(2,2)` stays at 1/9 in all 120 candidates up to numerical error. The corresponding truths are 13/64, 1/4, and 19/64. Increasing sample size cannot remove this discrepancy through the update alone.

The record-law TV error was approximately halved when parent size increased from 400 to 1600, but reference parent-survival RMSE retained contributions from invisible mass and within-record allocation. This residual error is consistent with the confidence set containing the truth. A single point-estimate representative must be distinguished from the set of distributions allowed by the class.

## 6. Certified projections for three representative samples

For each lambda, fix `N_parent=400, rep=0` and project target (1,1). Keep that sample's learning prediction, tail conditions, and threshold unchanged. Running the unchanged Stage 20 engine with 15 nodes did not achieve the requested accuracy in any of the three cases. Retaining those certificates as ancestors, use the unchanged Stage 21 engine in relative-scale coordinates. Every case uses endpoint-width tolerance 1e-3, node limit 65, and proposal-iteration limit 150.

All three met the tolerance, using 23, 13, and 15 nodes. The new certificates comprise 51 nodes, 54 dual certificates, 594 full-column inequalities, 82 tangents, and six inner endpoint parent distributions.

Exact rational endpoints are stored in `outputs/refined_projection_results.json` and in the individual certificates. In all three cases the outer range includes the true survival probability. These three results are not extended into a projection-coverage rate for all 120 samples. The initial 15-node failures to attain accuracy remain in the archive rather than retaining successful computations only.

## 7. Independent checks and reproducibility

Public counts were reconstructed from saved latent arrays and agreed with all 120 inputs. Increasing undetected components and latent values of D objects left all 120 public inputs unchanged. An independent exact-coordinate record map was compared at 2880 points. Direct density integration and the finite operators agreed exactly for 54 record probabilities. A separate verifier recomputed true masses from antiderivatives of the cell density.

The 120 learning predictions and intervals for log E0 were recomputed, and final original-EM masses and survival quantities were checked. The maximum invisible-mass difference from initialization was 6.245e-16. The raw minimum likelihood increment, -4.5475e-13, was retained as a roundoff diagnostic rather than repairing the history. During projection verification, LP, nonlinear optimization, and vertex enumeration were prohibited; ancestor certificates and the Stage 21 independent verifier rechecked all three cases.

All 120 sample generations, predictions, point estimates, and membership checks were rerun into a separate output directory. Aggregate results excluding time, 240 JSON files, two CSV files, and 600 latent arrays agreed. Byte identity including timestamps inside compressed NPZ files is not claimed. At the time this theory note was written, the three projection searches had not been rerun; independent rechecking of their certificates was a separate operation.

## 8. Validation remaining at this historical stage

This was an initial integrated pilot for non-nested limits and continuous parents. With 20 replicates per setting, it was not a precise Monte Carlo estimate of coverage. At this stage, further work included sharp confidence projections for all 120 samples, various depths and dependence structures, small inclusion probabilities, physical distance distributions, and photometric errors. Where the full parent probability cannot be recovered as a point, initialization dependence should remain visible and the distributions and target ranges compatible with the same observation law should be reported.

See [RESUMPTION.md](RESUMPTION.md) for the later additional rerun of the three projection searches. The historical scope of this theory note has been preserved rather than silently updated.
