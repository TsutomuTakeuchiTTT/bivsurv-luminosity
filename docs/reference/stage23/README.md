# Stage 23: paired observation-depth and recording-resolution comparisons on saved parent samples

> English translation of [`reference_experiments/stage23/README_ja.md`](../../../reference_experiments/stage23/README_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision. Paths and commands refer to the Stage 23 snapshot. The disposable-copy instruction is an added packaging safeguard, not a change to the original experiment.

The saved execution used Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, and SymPy. The inputs are exact copies of the 120 Stage 22 parent samples. No new random numbers are generated. Reobserving the all-data samples gives 720 design evaluations, but there are still only 120 independent samples.

Run the following commands from `reference_experiments/stage23/` in a disposable working copy. They reproduce the original commands and write to that copy's `outputs/`.

```bash
python depth_study.py --outdir outputs
python population_fibers.py --outdir outputs
python audit_controls.py --outdir outputs
python verify_stage23.py --outdir outputs
```

`depth_study.py` computes learning predictions, point estimates, and true-parent membership. `population_fibers.py` solves the separate population-identification problems and saves rational primal/dual certificates. `verify_stage23.py` checks the results without LP, nonlinear optimization, random generation, or numerical EM refitting. The two learning-EM steps are recomputed from their definition.

- `outputs/summary.csv`: aggregates over 20 replicates for each of 36 settings.
- `outputs/trials.csv`: 720 observation-design evaluations and 2160 reference survival values.
- `outputs/population_ranges.csv`: 54 population identified ranges.
- `outputs/population_certificates/`: 108 endpoint parent distributions and full-column dual certificates.
- `outputs/independent_verification.json`: independent verification results.
- `outputs/structural_controls.json`: preservation/loss of recorded information and violations of the tail constraints.
- `outputs/reproducibility.json`: comparison of the two executions recorded at this stage.
- [THEORY.md](THEORY.md): English translation of the observation design, initial measure, proofs of the main results, and scope of application.
- `stage22_record.tex`: Japanese LaTeX insertion for the preceding stage; environments and references were statically checked, but it was not compiled at this stage.

The true support, density, D counts, and total parent count are used for evaluation, not for fitting the conditional likelihood. The original EM does not enforce the tail constraints. The 13 candidates outside the class are recorded without repair. The 719 acceptances of the true parent are a finite-replicate outcome, not an exact coverage probability. Population identified ranges are not confidence intervals. Sharp sample-confidence projections for all 720 evaluations were not added at this stage.

See the English [results summary](SUMMARY.md).
