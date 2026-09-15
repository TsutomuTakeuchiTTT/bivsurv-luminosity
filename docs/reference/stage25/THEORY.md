# Stage 25: global likelihood gaps for constrained candidates

> English translation of [`reference_experiments/stage25/THEORY_ja.md`](../../../reference_experiments/stage25/THEORY_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

## 1. Verification target

Keep the 13 Stage 24 values of `independent_candidate.mass` completely unchanged and verify their gaps from the global maximum of the same all-data conditional likelihood. Preserve the 16 observation-equivalence classes covering the full real plane, public count tables, and tail upper bound 1/2. This calculation does not refit the 1040 Stage 24 GEM updates. Learning predictions, the validation-likelihood threshold, the parent confidence sets, and coverage of the truth are also unchanged.

The original parent class is `K = {m >= 0, 1'm = 1, V_j m >= 1/2, j=1,2,3}`. Define

```text
L(m) = product_i [(A_i m)/(V_i m)]^{n_i},
N = sum_i n_i,
f(m) = log L(m)/N.
```

Only positive-count records contribute factors to the product. Possible records absent from the sample are not removed from the observation model. K is compact, and L is continuous on K because denominators are at least 1/2. The saved candidate has positive likelihood, and the global maximum of L is attained. This does not conclude that the saved candidate attains it exactly.

## 2. Remove the completely invisible scale only for objective-value verification

Here the deep third layer's visibility row `eta = V_3` represents the union of all visible sets. The identity `eta_h = max_j V_jh` is checked on all 16 columns. Completely invisible columns have zero values for all A_i and V_j. Normalize the visible part using `r = eta'm > 0` and `p_h = m_h/r` for `eta_h = 1`. Then

```text
(A_i p)/(V_i p) = (A_i m)/(V_i m),
V_j p = (V_j m)/r >= V_j m >= 1/2,
V_3 p = 1.
```

Every original feasible parent has a visible representative with the same likelihood. Conversely, a normalized feasible visible distribution belongs to the original parent class. The two problems therefore have the same global maximum. This leaves 15 visible masses as optimization variables and only q_1 and q_2 as denominator-box coordinates.

This is an equivalent computational representation of the objective value, not a rewriting of the 13 original candidates or an assertion that invisible mass is truly identified as zero. Its feasibility does not extend unconditionally to models adding a positive lower bound on invisible mass or general moment constraints.

## 3. An upper bound covering one box

Write `q_j = V_j p`, `B = [l_1,u_1] x [l_2,u_2]`, a subset of `[1/2,1]^2`, and

```text
P_B = {p >= 0, 1'p = 1, l_j <= V_j p <= u_j},
omega_i = n_i/N,     lambda_j = n_j/N.
```

For l < u, the secant of -log q is

```text
s_j(q) = -(u-q) log(l)/(u-l) -(q-l) log(u)/(u-l),
```

and is -log l when l = u. Within the box, `-log q <= s_j(q)`, with error bounded by `(u-l)^2/(8l^2)`. Thus

```text
f_B(p) = sum_i omega_i log(A_i p) + sum_{j=1,2} lambda_j s_j(V_j p)
```

is a concave majorant of f(p). The third-layer denominator is 1 and contributes zero.

At a reference mass nu for which every positive-count numerator is positive, apply

```text
log(A_i p) <= log(A_i nu)+(A_i p)/(A_i nu)-1.
```

Together these give the affine majorant `f(p) <= T_B(p) = t'p`. The constant is added to all coefficients using `1'p = 1`. The reference need not lie in the box or maximize the objective.

Natural logarithms are enclosed in rational intervals. Use the upper endpoint for numerator logs and upper endpoints of -log l and -log u for denominator secants. If q is in the box, interpolation weights are nonnegative and preserve the upper bound. The same affine expression is evaluated at simplex basis vectors even when they lie outside the box; it is not required to bound the original likelihood at those outside basis points.

## 4. Full-column LP dual certificates

Collect box constraints as `G p <= b`, with rows V_1, -V_1, V_2, -V_2 and right-hand side u_1, -l_1, u_2, -l_2. For nonnegative y and unrestricted z, if

```text
G' y + z 1 >= t     on all 15 columns,
```

then `f(p) <= t'p <= b'y+z = U_B`. Convert the numerical LP's multipliers to nonnegative rationals, then reset `z = max_h[t_h-(G'y)_h]` to satisfy every column inequality exactly. Neither the LP success flag nor a tolerance on the dual residual is used as the guarantee. If multiplier proposals are unavailable, y = 0 still gives a valid, although possibly loose, bound.

Weak duality and concave tangents are established principles of convex analysis. The source note cites Boyd and Vandenberghe, *Convex Optimization*, Sections 3.1.3 and 5.1.3 ([book website](https://web.stanford.edu/~boyd/cvxbook/)). The added verification here applies them to the original conditional likelihood of the same 13 cases and connects full-domain coverage to a lower bound at the original candidate.

## 5. Covering the full domain and the final gap

Bisect the root `[1/2,1]^2` at coordinate midpoints and retain a tree whose two closed children cover each parent. Every leaf has its own bound U_B, so `U = max_leaf U_B` is a global upper bound. Evaluate the original saved 16-mass candidate m24 rationally as `f(m24) in [f_lo,f_hi]`; then

```text
0 <= f_star-f(m24) <= U-f_lo.
```

Exact equality of observation ratios with the normalized visible representative is checked as well. Stop when `U-f_lo <= 1e-6`, with a node limit of 501.

SLSQP is used only to propose tangent reference points for the concave relaxation. The internal proposal bound of 1e-12 on numerators is not imposed on the original model. Between tangents at the original candidate and at a proposed point, retain the smaller certified upper bound. No new point-estimate candidate is selected in this verification.

## 6. Main results

All 13 unchanged candidates achieved a certified mean-log-likelihood global gap of at most 1e-6. There were 537 nodes and 275 leaves. The largest certified gap was 9.320517628688491e-7; the largest total-log-likelihood gap was 2.7588732180917934e-4. This certifies approximate global maximization, not exact equality to the MLE. It also does not certify survival-probability error or coverage error of 1e-6.

If the saved mean-log-likelihood increase from the Stage 24 80-step GEM iterate to the original candidate is enclosed by [a,b], the global gap after 80 GEM updates is enclosed by `[a, b+(U-f_lo)]`. All cases have a positive lower bound. The smallest lower bound was 2.8312709804092922e-5, and the largest upper bound across cases was 2.49651765810805e-4. The fixed-80-step result and the separately obtained high-accuracy candidate serve different purposes.

## 7. Budget/accuracy controls and the distinction from identification

Restricting the first case to one node left a gap of approximately 0.0202101 and did not meet the tolerance. Requesting 1e-8 in that same case with 501 nodes retained a valid gap upper bound of 1.7421461047523004e-8 but did not attain the request. Finite termination at arbitrary precision or practical complexity for arbitrary models is not claimed.

For `coarse_d3__lambda0_N400_r010`, use the visible representative Q to form

```text
P_theta = (1-theta)Q + theta R_I,     R_I(I) = 1.
```

These parents belong to the same tail class for `0 <= theta <= approximately 0.006823529411754671`, and their conditional likelihoods are exactly equal. Three rational parents were saved, at zero, the midpoint, and the upper endpoint. They share the same approximate global guarantee but have different values of `S(d,d) = theta`. Optimization certification and unique distributional identification are separate. This family is not asserted to equal the full MLE family.

## 8. Verification scope

The independent verifier reconstructs the model from the original public records and recomputes candidate invariance, the invisible-scale identity, the log-likelihood lower bound, coverage by the root and all child boxes, and 8055 dual column inequalities over all nodes. Rather than invoking the search-side tangent-coefficient builder or LP, it directly constructs the upper-bound value at each basis vector. Foundational rational-logarithm arithmetic and the interval-observation model are shared with ancestor code. No LP, nonlinear optimization, EM refitting, root search, or new sample generation is performed by this verifier. Seven types of harmful alteration are rejected.

For another 312 deterministic feasible parent masses, 4848 ratio-scale identities and 936 inclusion-probability increases were checked exactly. Their likelihoods were also checked against the returned global upper bounds. All 13 searches were rerun into a separate output directory, and certificate fields agreed apart from time. The previous stage's certificates for 1040 updates were rechecked too, but that GEM search was not rerun. No new coverage rates, parent samples, learning predictions, or confidence projections were added.
