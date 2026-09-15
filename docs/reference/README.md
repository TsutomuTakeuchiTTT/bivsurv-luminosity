# Archived reference experiments: English reading guide

These notes describe the preserved implementations and evidence in [`reference_experiments/`](../../reference_experiments/). The English translations are deliberately stored outside those snapshots so that the original files, source hashes, and certificate provenance are unchanged. Folder names such as `vendor`, `prior`, and `reference_stage24` identify earlier snapshots from this project, as recorded in [`docs/provenance.json`](../provenance.json).

| Stage | Purpose | English documents | Original snapshot |
|---|---|---|---|
| 22 | Continuous-parent pilot: 120 samples, membership checks, three certified survival projections | [README](stage22/README.md), [theory](stage22/THEORY.md), [summary](stage22/SUMMARY.md), [resumption](stage22/RESUMPTION.md) | [stage22](../../reference_experiments/stage22/) |
| 23 | Six paired depth/resolution designs on the same 120 samples; separate population-identification projections | [README](stage23/README.md), [theory](stage23/THEORY.md), [summary](stage23/SUMMARY.md) | [stage23](../../reference_experiments/stage23/) |
| 24 | Constrained M-step and first-order diagnostics for the 13 Stage 23 tail-class failures | [README](stage24/README.md), [theory](stage24/THEORY.md), [summary](stage24/SUMMARY.md) | [stage24](../../reference_experiments/stage24/) |
| 25 | Global likelihood-gap certificates for the same 13 Stage 24 candidates, without changing their masses | [README](stage25/README.md), [theory](stage25/THEORY.md), [summary](stage25/SUMMARY.md) | [stage25](../../reference_experiments/stage25/) |

## Read the stages in their historical scope

Stage 22's original summary and theory were saved before a later resumption. They explicitly say the projection searches were not rerun in that earlier verification. The separate resumption note then records that both generations of those searches were rerun. Both statements are preserved with their temporal scope; neither is silently substituted for the other.

Stage 24 certifies individual surrogate M-steps and constrained first-order conditions. It does not claim a global optimum for the original likelihood. Stage 25 adds a global **near-optimality** bound for those same 13 candidates. Its 1e-6 mean-log-likelihood guarantee is not a bound on statistical survival-probability error.

The new Gaussian luminosity visualization experiment is a separate experiment. Its 240 high-resolution D-completion fits do not receive Stage 25's certificates merely because both implementations use D-completion. For those figures, read the root [README](../../README.md), [results](../RESULTS.md), and [figure guide](../FIGURES.md).

## Running a preserved stage

Use the working directory specified by that stage's English README. Original commands and output names are retained so they can be matched to the code. To protect distributed evidence, run regeneration or verification in a disposable copy; several verifiers write a report even though they do not rerun numerical searches. See the [reproducibility guide](../REPRODUCIBILITY.md) for the distinction.

The nested code is supplied for standalone reproducibility, not as a single polished public API. `reference_experiments/stage22/resumption_checks/compare_runs.py` retains assumptions about an earlier local working directory and is not a portable entry point. That archival helper has not been rewritten as part of translation.

## Historical Japanese LaTeX records

`stage21_record.tex`, `stage21_record_reissued.tex`, `stage22_record.tex`, `stage23_record.tex`, and `stage24_record.tex` are preserved Japanese manuscript/work-note fragments. They are not installed software manuals and are not required to follow the English execution instructions. The English theory and result documents explain each included Stage 22–25 experiment. The complete Japanese paper and those LaTeX records have not been declared translated by this documentation edition.
