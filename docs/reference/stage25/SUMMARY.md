# Stage 25 results

> English translation of [`reference_experiments/stage25/SUMMARY_ja.md`](../../../reference_experiments/stage25/SUMMARY_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

Keeping the same 13 constrained Stage 24 candidates unchanged, global upper bounds were constructed for their differences from the optimum of the original conditional likelihood over the entire parent class.

| Quantity | Result |
|---|---:|
| Target candidates | 13; original masses unchanged |
| Parent model | 16 full-real-plane observation classes, tail upper bound 1/2 |
| Main requested accuracy | Mean-log-likelihood gap 1e-6 |
| Accuracy attained | 13/13 |
| Largest certified gap | 9.320517628688491e-7 |
| Largest total-log-likelihood gap | 2.7588732180917934e-4 |
| Total box nodes | 537 |
| Final leaves | 275 |
| Independently checked column inequalities | 8055 |
| Rejected harmful alterations | Seven types |

Normalization away from completely invisible mass was used only as a computational representation preserving the optimal likelihood value. The original saved candidate masses were not rewritten to set invisible probability to zero. The visible representative has 15 masses and two nonconstant denominators.

Stage 24's largest first-order residual was 5.85e-14. That number was not used as a global-optimum error; separate box, tangent, and dual certificates were constructed. The result certifies global near-optimality to 1e-6, not exact equality to the MLE or proximity of the survival function to its truth within 1e-6.

Combining the result with saved likelihood-increase certificates also gave a positive lower bound and finite upper bound on the global gap of every 80-step GEM candidate. All 13 lower bounds were positive; 80 fixed updates are not treated as completed optimization.

Budget control: limiting the first case to one node gave a gap of approximately 0.0202101 and did not meet the tolerance. Stricter-accuracy control: requesting 1e-8 in the same case with 501 nodes gave approximately 1.74215e-8, also not attaining the request. Valid bounds and attainment of requested precision were reported separately.

Additional checks covered 312 deterministic feasible candidates, 4848 ratio identities, and 936 inclusion-probability monotonicity comparisons. Three equal-likelihood parents differing in invisible mass were also constructed, preserving the distinction between global near-optimality and identification. All 13 searches were rerun, and all certificate fields agreed except execution time.

No new GEM refits, parent samples, learning predictions, parent confidence sets, coverage estimates, or survival-confidence projections were generated. Independent rechecking of Stage 24's 1040 updates was saved separately in `outputs/stage24_recheck.json`.
