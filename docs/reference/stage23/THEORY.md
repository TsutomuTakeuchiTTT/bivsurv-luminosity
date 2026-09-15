# Stage 23: paired comparisons of observation depth and recording resolution

> English translation of [`reference_experiments/stage23/THEORY_ja.md`](../../../reference_experiments/stage23/THEORY_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

## 1. Scope and relationship to Stage 22

Reuse the 120 Stage 22 parent samples containing 120000 objects. Keep latent values, limit-layer assignments, and training/validation membership unchanged. Generate no new random numbers. Applying six observation designs to the same 120 parent samples gives 720 design evaluations, not 720 independent parent samples. D counts and parent sample size are not used in the likelihood.

Each design has limits (2,1), (1,2), and (d,d), with d = 2, 3, 4. The two recording schemes use the original boundary {1} or the boundaries {1,3/2,2,3,4}. The latter preserves all target coordinates, (1,1), (3/2,1), and (2,2), and the earlier depth thresholds. Candidates are Borel distributions on the entire real plane; they are not truncated to the true support (0,4)^2.

All designs impose the common conditions `P(D_(2,1)) <= 1/2` and `P(D_(1,2)) <= 1/2`. The additional region D_(d,d) is a subset of both fixed D regions. Consequently `P(D_(d,d)) <= 1/2` is redundant and the parent class is identical across designs. The bound for the deep layer is not replaced by its true D probability.

## 2. Fixing the initial measure

Do not reset to uniform masses for every representation. Use the same continuous measure R with independent coordinates and marginal CDF

```text
F_R(x) = 1/[3(2-x)]      (x <= 1),
         x/3             (1 < x <= 2),
         1-1/[3(x-1)]    (x > 2).
```

It is continuous and nondecreasing on the entire real line, allocating 1/3 to each of `(-inf,1]`, `(1,2]`, and `(2,inf)`. Thus every original observation class has mass 1/9, reproducing Stage 22's initialization. The split at 1.5 within (1,2] also gives equal masses, matching the earlier reference survival values. New observation-class masses are formed by summing exact cell integrals under R. Internal allocation for the reference survival values also retains R's conditional distributions. R specifies initialization; it is neither a Bayesian prior nor a shape restriction on candidate parents.

For every d, the completely invisible region `I_d = (d, infinity)^2` forms one class. The original EM preserves `R(I_d) = 1/[9(d-1)^2]`, namely 1/9, 1/36, and 1/81. The reduction reflects the shrinking invisible region in each design, not estimation of its total mass from observations.

## 3. Point estimation and true-parent membership

The all-data point estimate uses the unchanged original EM with mean positive-score tolerance 1e-8 and at most 4000 updates. Learning prediction applies two rational EM updates to training records only, starting from the same initial measure. For validation-likelihood membership, log E0 is enclosed by rational intervals. The significance level is 1/20. The true distribution's tail-class membership is checked separately by integrating the continuous density. The all-data fit is not fed back into learning prediction.

Record separately the candidate's TV error against the true record law, reference RMSE at three targets, numerical stopping, final tail-class membership, and preservation of initialized invisible mass. Do not remove a sample from coverage summaries when its original-EM candidate falls outside the tail class. Enforcing tail conditions in point estimation requires the separate constrained M-step developed in Stage 19. This stage diagnoses the unconstrained EM without post-hoc repair.

## 4. Compute population identified sets separately

Evaluate identification with the entire true record law pi0 known exactly, separately from the sample confidence set. Using complete cell-equivalence classes that also preserve survival coefficients c, impose

```text
sum(m) = 1, m >= 0,
(A_jr - pi0_jr V_j)m = 0   for every record,
V_j m >= 1/2,
```

and minimize or maximize `c^T m`. This is a linear program. A condition pi0 = 0 is an input from the true population law; it does not identify an unobserved sample category with zero probability.

The LP proposes numerical candidates. Solve the active-face equalities with rational linear algebra to reconstruct endpoint-attaining parent masses. Convert dual multipliers to rationals as well, and correct the normalization multiplier by the minimum full-column residual so that exact weak duality holds. Do not repair or relax the original parent-mass constraints.

Specifically, for a minimization problem with `E m = b`, `G m <= h`, and `m >= 0`, verify `E^T y + G^T z <= c` for unrestricted y and nonpositive z. Then

```text
b^T y + h^T z <= c^T m
```

holds for every candidate. Using objective -c gives an upper bound for maximization. Here `G = -V`, `h = -1/2`, and the first row of E is normalization. A numerical solver's success status is not a substitute for the guarantee.

## 5. A closed-form identified range when recording preserves the boundaries

Let `Q_d = P(. | I_d^c)`. Assume each interval record in the deep layer (d,d) determines the target survival indicator and the inclusion indicators of the two shallow layers. The fine-record design satisfies this condition on every cell.

The deep layer's population record law fixes

```text
s_d = Q_d(G_t),     v_j = Q_d(D_j^c).
```

Decompose the parent as `P = (1-theta) Q_d + theta R_I`. For targets `t <= (d,d)`, `I_d` is a subset of `G_t`, giving

```text
S_P(t) = theta + (1-theta)s_d,
(1-theta)v_j >= 1/2.
```

Therefore

```text
0 <= theta <= theta_max = 1-max_j[1/(2 v_j)],
```

and the sharp identified range is

```text
[s_d, theta_max+(1-theta_max)s_d].
```

The lower endpoint is attained at theta = 0 and the upper endpoint at theta = theta_max. To realize Q_d with the same deep-record law, one may use the true visible conditional distribution; R_I may be any probability measure supported within I_d. This construction is for proving and evaluating the identification result. It does not supply the unknown truth to the estimator.

In terms of the evaluation truth, write `d_0 = P0(I_d)`, `qmin = min_j P0(D_j^c)`, and `S0 = S_P0(t)`. Then

```text
s_d = (S0-d_0)/(1-d_0),
upper = 1-(1-S0)/(2 qmin).
```

In this design qmin is fixed by the shallow layers and does not depend on d. Greater depth brings the lower endpoint closer to the truth, while the upper endpoint can remain unchanged. In particular, a true I_d mass of zero at d = 4 does not determine the candidate parents' I_d mass to be zero, because candidates are not restricted to the true compact support.

For lambda = 0 and target (2,2), the exact ranges at d = 2, 3, 4 are

```text
[0,2/5], [1/5,2/5], [1/4,2/5].
```

The final width 3/20 remains even with zero sampling error. All true deep-layer objects being A does not prove candidate inclusion probability 1 under the conditional observation model, which does not use the number of unrecorded D objects.

## 6. Greater depth does not necessarily preserve old information

With coarse recording, for example, (5/4,5/4) and (5/2,5/2) produce the same high–high A record at the new limit (3,3). At the old limit (2,2), the first is A and the second is D. The old record cannot be reconstructed from the new one. Even when detection improves, discarding the details into coarse interval labels can discard information about the old thresholds.

For lambda = 0 and target (2,2), the coarse-record population identified ranges are

```text
d = 2: [0,2/5],     d = 3: [2/15,1/2],     d = 4: [3/16,1/2].
```

The upper endpoint increases from 2/5 to 1/2 despite the common parent class. This is a concrete example in which increasing depth does not imply nested identified sets. Each endpoint was checked using a rational realizing parent mass and a matching dual certificate.

By contrast, the fine recording scheme retains the old thresholds 2 and 3. The new record determines the old record, including exclusion as old D. Population identified sets are then nested. This does not assert that split-likelihood confidence sets from individual samples are nested.

## 7. Interpreting the execution results

The true parent was accepted in 719 of 720 design evaluations and rejected once. The rejected case is `aligned_d2`, lambda = 3/4, parent size 400, replicate 7, with `log E0` approximately 3.68197195 > log20. The 720 evaluations are not independent Bernoulli trials and must not be combined as an independent binomial experiment. Each parent-distribution/sample-size/observation-design setting is a 20-replicate pilot, not a precise coverage estimate.

The original EM stopped numerically in all 720 evaluations, but 13 candidates were outside the common tail class: five in `coarse_d3` and eight in `coarse_d4`, all with 400 parent objects. The maximum violation of the inclusion-probability lower bound was approximately 0.07364, not merely least-significant-digit roundoff. These candidates were retained as diagnostics and not called class-constrained MLEs.

Both endpoints of all 54 population projections were enclosed by rational primal/dual certificates. The maximum bracket width was `1/2703999999937808`, approximately 3.70e-16. The 27 fine-record projections were also checked with the closed form above and 54 explicit endpoint parent distributions. Sharp confidence projections for all 720 samples were not computed. The newly added optimization concerns population identification, not coverage.

## 8. Synthesis at this stage

The current model combines D-completion updates, an exact cell representation of the observation law, preservation of unidentified degrees of freedom, parent confidence sets, and projections with computational guarantees. The comparison makes clear that recording rules, not only observation depth, must be part of the specification. General unconstrained-EM consistency, ungrouped values and photometric errors, and actual distance/selection dependence are not inherited without conditions. These belong to distinct scope statements and extensions when formulating the main theorems, reported quantities, and returned states on failure.
