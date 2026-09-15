# Stage 24 reproducibility package

> English translation of [`reference_experiments/stage24/README_ja.md`](../../../reference_experiments/stage24/README_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision.

This stage addresses the 13 design/sample cases outside the explicit tail class in Stage 23. They arise from eight existing parent samples. No new random samples are generated.

## Execution

The saved calculations were verified with Python 3.13.5, NumPy 2.3.5, and SciPy 1.17.0, also using the standard-library `fractions.Fraction` and `decimal` modules. Run from `reference_experiments/stage24/`:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit_stage24.py --outdir new_outputs
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_stage24.py --outdir new_outputs
```

The default is 80 GEM updates per case, each with rational verification. Original results are in `outputs/`; rerun agreement is recorded in `outputs/reproducibility.json`. Each case's 80 updates and separately constructed boundary candidate are stored in `outputs/certificates/`.

`verify_stage24.py` does not call optimization. It reconstructs the record cells and the definition of the two-step learning EM. It does not invoke the M-step proposal function: it reconstructs and checks the solution from the KKT multipliers.

## Distinguishing the outputs

- `GEM_steps`: certificates of 80 likelihood-increasing updates from a feasible initial measure. Convergence after exactly 80 steps is not claimed.
- `independent_candidate`: a candidate obtained by separately optimizing the same constrained likelihood. Its feasibility, likelihood increase relative to the 80-step iterate, and first-order feasible-direction residual bound are checked. It is neither a single GEM update nor certified here as a global MLE.
- `comparison.csv`: diagnostics for the old out-of-class candidate, the 80-step GEM candidate, and the separately obtained constrained candidate.

Observed inputs and old fits are preserved unchanged in `inputs/`, with original file hashes in `PROVENANCE.json`. The three source files in `vendor/` are unchanged copies from earlier validation packages. This directory is self-contained for the commands above; executing all earlier stages is not required.

`stage23_record.tex` is the Japanese LaTeX record prepared for insertion at this stage. The complete manuscript and BibTeX database were not changed. Read the English [theory](THEORY.md) and [summary](SUMMARY.md) for the documented method and results.
