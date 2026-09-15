# Stage 22 recovery and additional rerun record

> English translation of [`reference_experiments/stage22/RESUMPTION_ja.md`](../../../reference_experiments/stage22/RESUMPTION_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

The Stage 21 LaTeX record and Stage 22 results were recovered from files saved before an interrupted response. Recovered original outputs remain in `outputs/`; the additional rerun is separated into `resumption_checks/sampling_rerun/`.

## Operations performed on resumption

1. The 20 source files in the saved package were checked against their SHA-256 hashes in `PROVENANCE.json`.
2. The independent verifier rechecked the 120 samples, 120000 objects, 120 learning predictions, 120 true-parent confidence-set membership decisions, and projection certificates for three representative samples.
3. With the same seeds and input design, all 120 samples were regenerated through individual-object generation, observation selection, the two-step learning EM prediction, the all-data EM point estimate, and confidence-set membership.
4. For the three predeclared samples, one at each lambda with 400 parent objects and replicate 0, both the initial Stage 20 projection search and the improved Stage 21 projection search were rerun.
5. All rerun outputs were checked by the independent verifier.
6. Comparisons with the original outputs confirmed byte identity for 120 public-input JSON files, 120 fit JSON files, and two CSV files, and content identity for 600 latent arrays. The six projection certificates and aggregate JSON agreed apart from execution time.

`outputs/reproducibility.json` is the reproducibility record saved before the interruption. The additional comparison on resumption is recorded in `resumption_checks/resumption_reproducibility.json`. These execution scopes are distinct.

## Distinguishing the results

The true parent belonged to the confidence set in all 120 samples. This is an observed result with 20 replicates per setting, not a claim that coverage probability is 1 and not exact coverage obtained by enumerating every count table. Sharp survival projections were computed for only the three predeclared samples and one target per sample; the remaining 117 samples were not fully projected.

All three projection calculations achieved endpoint brackets of width at most 1e-3. They required 23, 13, and 15 nodes. In total, the checks covered 51 nodes, 54 dual certificates, 594 column inequalities, 82 tangents, and six inner parent distributions. The initial 15-node calculations that did not meet the requested precision were also preserved.

Point-estimation EM stopped numerically for every one of the 120 samples, but no global-MLE certificates were added for those fits. The estimated mass of the completely invisible region `(2, infinity)^2` retains its initialized value 1/9. The true masses at lambda = -3/4, 0, 3/4 are 13/64, 1/4, 19/64, respectively. Improving the fit to the observation law is distinct from point recovery of the entire parent distribution.

## Additional files

- `stage21_record_reissued.tex`: the Stage 21 Japanese record reissued on resumption.
- `resumption_checks/resumption_reproducibility.json`: agreement checks for the additional rerun.
- `resumption_checks/projection_endpoint_table.csv`: rational projection endpoints and outward-rounded eight-decimal representations.
- `resumption_checks/sampling_rerun/independent_verification.json`: independent verification of the additional rerun.
- `resumption_checks/stage21_latex_lint.json`: static checks of LaTeX environments, duplicate labels, and internal references. Compilation was not performed in that check.

The earlier `THEORY_ja.md` and `SUMMARY_ja.md` were retained as originally saved. In particular, the statement in the theory note that the three projection searches were not rerun refers to work completed before the interruption. The additional rerun described here did rerun those three projection searches.
