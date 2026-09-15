# Stage 24: verification of the constrained M-step for 13 out-of-class candidates

> English translation of [`reference_experiments/stage24/THEORY_ja.md`](../../../reference_experiments/stage24/THEORY_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

## 1. Targets and quantities left unchanged

The analysis includes all 13 Stage 23 tail-class failures: five `coarse_d3` and eight `coarse_d4` cases. Because these are reobservations of the same parent samples at different depths, there are eight distinct parent samples, not 13 newly generated independent samples.

Observation intervals, limits, training/validation counts, observation-equivalence classes over the full real plane, the common initial measure R, and the tail class remain fixed. All-data counts are used for point estimation, but the result is not fed back into learning predictions. Stage 23's two-step learning EM predictions and confidence sets are unchanged. True densities and latent values are not estimator inputs. Reference RMSE uses the three previously stored true target values only after fitting.

The 13 original-EM results remain nonnegative, normalized survival distributions. The issue is that they are outside the separately declared tail class, not that they are invalid survival distributions. Old results are not overwritten.

The parent class is `K = {m >= 0, 1'm = 1, D_j m <= 1/2, j=1,2,3}`. In these designs, D_3 is a subset of D_1 intersect D_2, so the third constraint is redundant. Each D_j is a membership-indicator row, `q_j = 1-D_j m`, and every case has 16 observation-equivalence classes.

## 2. The unchanged E-step and a concave M-step

With `a_i = A_i m`, `q_i = V_i m`, and positive counts n_i,

```text
nu_h = m_h sum_i n_i[A_ih/a_i + (1-V_ih)/q_i],
K = sum_i n_i/q_i = sum_h nu_h,     w_h = nu_h/K.
```

The constrained M-step maximizes `Qbar(u|m) = sum_h w_h log u_h` over the parent class. Although the parent likelihood is nonconcave, this surrogate is concave. The existing completion-EM lower bound gives, for feasible m and u,

```text
ell(u)-ell(m) >= K[Qbar(u|m)-Qbar(m|m)].
```

No new missing-data mechanism is introduced at this stage. As in the source note, K in the displayed normalization denotes expected total mass; the feasible parent class is specified separately above.

## 3. Aggregation into four D-membership groups

Group cells by `s = (D_1h,D_2h) = (a,b)` in {0,1}^2. Let `W_ab = sum_{h:s=(a,b)} w_h` and `p_ab = sum_{h:s=(a,b)} u_h`. When every W_ab is positive, the optimal within-group allocation is

```text
u_h = p_ab w_h/W_ab.
```

Thus the M-step reduces exactly to a four-variable problem:

```text
maximize sum_ab W_ab log p_ab,
p >= 0, sum p = 1,
p_10+p_11 <= 1/2,     p_01+p_11 <= 1/2.
```

This does not assume independence of the D flags. All A, B, C records and the expected completion quantities are retained in W_ab.

There are four candidate active sets:

- Inactive: `p_ab = W_ab`.
- Only the first constraint active: with `t = W_10+W_11`, `p_1b = W_1b/(2t)` and `p_0b = W_0b/[2(1-t)]`.
- Only the second active: exchange the two coordinates in the previous formula.
- Both active: `p_00 = p_11 = (W_00+W_11)/2` and `p_01 = p_10 = (W_01+W_10)/2`.

Determine the active set by exactly checking feasibility, nonnegative beta_1 and beta_2, and

```text
W_ab/p_ab = eta+beta_1 a+beta_2 b,
beta_1(p_10+p_11-1/2) = beta_2(p_01+p_11-1/2) = 0.
```

When these hold, concavity and complementarity show that no feasible p has a larger objective: they certify global surrogate maximization. Degenerate cases with zero group expectation and other values of the upper bound epsilon are outside this particular closed-form implementation. All four groups were positive in all 1040 updates across the 13 cases. This is not a restriction of the original general linearly constrained M-step to these special cases.

## 4. Finite-precision candidates and ascent certificates

Unrestricted rational iteration causes denominators to grow. At each step, the implementation first computes the exact M-step solution u*, then converts it to a rational candidate with 22 significant decimal digits. The candidate is normalized. Only if rounding slightly violates a tail condition is the minimally required rational mixture with the same fixed initial mass R used to make it feasible.

This constructs the new candidate before acceptance; it is not a post-hoc repair of an earlier out-of-class estimate. No positive mass floor is added to the model, and all original constraints are checked unchanged.

Accept u only after rational logarithm intervals certify

```text
Qbar(u|m)-Qbar(m|m) >= 0,
[ell(u)-ell(m)]/N >= 0,
[ell(u)-ell(m)]/N >= K/N [Qbar(u|m)-Qbar(m|m)].
```

Each comparison uses the appropriate lower and upper interval endpoints. The actual iteration is therefore a constrained generalized EM with certified increase, rather than exact-maximization EM.

Concavity supplies an upper bound on the surrogate-maximization error at each step:

```text
0 <= Qbar(u*|m)-Qbar(u|m)
  <= sum_h w_h(u*_h/u_h-1).
```

This can also be recomputed rationally. Eighty fixed updates were checked for each of 13 cases; no general convergence claim was made merely because 80 updates had been executed.

## 5. Separately constructed boundary candidates and constrained stopping

Independently of the 80-step GEM audit, maximize the same conditional likelihood by SLSQP from three deterministic starting points: the 80-step iterate, common R, and uniform mass on the 16 classes. Newton correction on the active boundary face, using the analytic Hessian, is used to propose candidates. This is not a GEM step requiring surrogate ascent. Convert each proposal to nonnegative rational mass and check the original constraints and observed-likelihood increase from the 80-step iterate. Every positive-count observed numerator must be positive; in particular, zero boundary masses are not replaced by pseudo-observations.

For the mean log-likelihood gradient `g = grad ell(m)/N`, `g'm = 0`. The quantity

```text
delta(m) = max_{v in K} g'(v-m)
```

is the largest first-order change allowed by the constraints. If nonnegative lambda and an unrestricted z satisfy

```text
D_1h lambda_1 + D_2h lambda_2 + z >= g_h
```

for every mass column, then

```text
0 <= delta(m) <= z+(lambda_1+lambda_2)/2.
```

After converting the LP-proposed lambda to nonnegative rationals, reset `z = max_h(g_h-lambda_1 D_1h-lambda_2 D_2h)` and check every column exactly. This constructs a dual certificate; it does not change the parent class.

Every one of the 13 candidates had a first-order residual upper bound below 1e-8; the largest observed bound was approximately 5.8471e-14. In contrast, the largest unconstrained mean positive score at the same candidates was approximately 0.032846. An unconstrained score counts directions that are not allowed at the boundary and must not be reused as a constrained stopping criterion.

This first-order residual is not the difference from the global optimum of a nonconcave likelihood. Agreement between several numerical starting points is not a global proof either. The global optimality proved at this stage is that of each surrogate M-step, not of 13 global MLEs for the original conditional likelihood.

## 6. Cases not resolved by invisible mass alone

Removing the completely invisible mass theta from an old candidate leaves the observation ratios unchanged, but changes inclusion probabilities to `q_j/(1-theta)`. Only one of the 13 cases enters the tail class through that operation. In the other 12, even removal of all invisible mass leaves `min q_j < 1/2`; scale adjustment alone cannot restore class membership. The audit starts all cases from the same feasible initial distribution rather than silently renormalizing the old candidates and reporting them as new estimates.

Under a constrained M-step, even a completely invisible cell can follow

```text
u*_h = m_h/(eta+beta_1+beta_2).
```

The original conservation rule returns in the inactive case, but preservation of the initial invisible mass generally no longer holds when constraints are active. This reflects maximization within an additional parent class, not newly observed information about D. Invisible mass near zero in these candidates does not imply identification of the entire distribution.

## 7. Interpreting this experiment

All 13 existing out-of-class cases were included, with 80 updates each and 1040 total. The same cases were rerun into another output directory, and certificate comparison excluded time. No new parent samples, limits, learning predictions, confidence sets, or coverage evaluations were added. The original 719/720 acceptance result and 54 population-identification projections were unchanged.

Reference RMSE decreased in all 13 cases, but this analysis selected class-violation cases only. It is not a general comparison of mean risk or coverage. A likelihood lower than that of an earlier unconstrained candidate is not contradictory because the feasible sets differ. Monotone increase is compared between feasible iterates in the same constrained class.

The conclusion at this stage is that an analytic M-step, certificates of feasibility and increase, and a constrained first-order stopping diagnostic were connected for point estimation under the declared tail conditions. Many general tail constraints, all zero-group degeneracies, and a global proof for the original likelihood remained separate tasks.
