# Stage 22 initial integrated-pilot results

> English translation of [`reference_experiments/stage22/SUMMARY_ja.md`](../../../reference_experiments/stage22/SUMMARY_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

A total of 120 samples and 120000 objects were generated from continuous parent distributions. The true parent was accepted in all 120 samples. There were 20 replicates per setting; this does not assert a coverage probability of 1.

| Lambda | Parent sample size | Mean TV | Reference parent-survival RMSE | Parent-set acceptance |
|---|---:|---:|---:|---:|
| -3/4 | 400 | 0.059385 | 0.076564 | 20/20 |
| -3/4 | 1600 | 0.029562 | 0.073520 | 20/20 |
| 0 | 400 | 0.063484 | 0.111616 | 20/20 |
| 0 | 1600 | 0.033199 | 0.108534 | 20/20 |
| 3/4 | 400 | 0.064485 | 0.148211 | 20/20 |
| 3/4 | 1600 | 0.027549 | 0.146064 | 20/20 |

Completely invisible mass remains at the EM initial value 1/9. Its true values are 13/64, 1/4, and 19/64. Improving the observational fit must not be equated with recovering the entire parent distribution.

The target (1,1) was projected for the first sample of 400 parent objects from each parent distribution. The unchanged Stage 21 engine achieved endpoint brackets of width at most 1e-3 in all three cases, using 23, 13, and 15 nodes. Sharp projections were not computed for the other 117 samples.

Independent verification rechecked all 120 public inputs, 120 predictions and true-parent membership decisions, the CDF identity of the generation formula for 120000 objects, and three projection certificates. The sample audit was rerun and 240 JSON files, two CSV files, and 600 latent arrays agreed. This original summary does not record a rerun of the projection searches.

The later [resumption record](RESUMPTION.md) documents the additional searches and checks separately; it does not alter the scope reported in this original summary.
