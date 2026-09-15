# Stage 23 results

> English translation of [`reference_experiments/stage23/SUMMARY_ja.md`](../../../reference_experiments/stage23/SUMMARY_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

The 120 Stage 22 parent samples, containing 120000 objects, were reobserved under six designs, keeping the same latent values, layer assignments, and training splits. The designs combine depths d = 2, 3, 4 with either coarse intervals or fine intervals preserving the target boundaries and old thresholds. No new random numbers were generated; this was not an experiment with 720 independent samples.

## Sample results

Three-target reference survival RMSE for lambda = 0 and 1600 parent objects, averaged over 20 replicates:

| Depth | Coarse records | Boundary-preserving fine records |
|---|---:|---:|
| 2 | 0.108534 | 0.108005 |
| 3 | 0.033827 | 0.026702 |
| 4 | 0.022774 | 0.014930 |

All 120 baseline inputs and learning predictions agreed exactly with Stage 22. The maximum point-estimate difference was 1.11e-16. The original EM stopped numerically for all 720 conditions, but this does not certify global optimality. The initialized completely invisible masses 1/9, 1/36, and 1/81 were preserved in their respective designs.

The true parent belonged to the confidence set in 719 evaluations. The only rejection was `aligned_d2/lambda2_N400_r007`. These are observed results from 20 replicates per condition, not exact coverage probabilities. Thirteen unconstrained-EM candidates lay outside the tail class and were recorded without repair.

## Population identified ranges

For lambda = 0 and target S(2,2):

| Depth | Coarse records | Boundary-preserving fine records |
|---|---|---|
| 2 | [0,2/5] | [0,2/5] |
| 3 | [2/15,1/2] | [1/5,2/5] |
| 4 | [3/16,1/2] | [1/4,2/5] |

These are identified ranges when the population probabilities of all records are known exactly, not confidence intervals. With coarse records, increasing depth does not necessarily preserve information about the old thresholds. Even with fine records, invisible candidate mass remains because the true support is not imposed on candidates. The deepest design retains parent-range width 3/20 even though its true D probability is zero.

Rational primal/dual certificates enclosed all 54 projections; the maximum endpoint error was approximately 3.70e-16. The 27 fine-record projections were also checked against the closed form and 54 separate explicit endpoint parent distributions.

## Checks

A verifier that does not invoke LP or nonlinear optimization checked 720 inputs, 720 predictions, 720 true-parent membership decisions, 54 projection certificates, 108 parent distributions, and 2268 full-column dual inequalities. Additional checks included 510 direct-integration comparisons, 5760 record-membership comparisons, 480 selection-set inclusions, 41 initial-mass aggregations, and 484 cell comparisons of the recorded-information maps. Four types of harmful certificate alteration were rejected.

All six sample calculations and all 54 projections were rerun. The 720 input JSON files, 720 fit JSON files, 54 certificate JSON files, and aggregate CSV files were byte-identical.

No new set of 720 sharp sample-confidence projections was computed. This validation separates observation depth, interval precision, initial conditions, and identifiable targets.
