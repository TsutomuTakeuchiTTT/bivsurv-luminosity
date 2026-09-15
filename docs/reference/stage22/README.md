# Stage 22: initial integrated pilot with a continuous parent distribution

> English translation of [`reference_experiments/stage22/README_ja.md`](../../../reference_experiments/stage22/README_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision. Paths and commands are relative to the archived Stage 22 directory, not to this documentation directory. The disposable-copy instruction is a packaging safeguard added here.

The snapshot contains generation, point estimation, and true-parent confidence-set membership for 120 samples, together with multi-limit survival projections for three representative samples. A Japanese LaTeX record of Stage 21 is also included.

## Rerunning

Run the following in `reference_experiments/stage22/` **in a disposable working copy**. The commands reproduce the original defaults and may write to that copy's `outputs/`.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit_stage22.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python project_stage22.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python refine_projections_stage22.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_stage22.py
```

On Windows, set the same environment variables using the syntax of the chosen shell. The scripts use Python, NumPy, and SciPy. Network access and external data are not required for the computations once the dependencies are installed. Version information is recorded in `outputs/stage22_sampling_results.json`. Numerical solver proposals can depend on the environment.

- `audit_stage22.py`: individual-object generation, public records, learning predictions, point estimates, and evaluation.
- `project_stage22.py`: initial certificates using the unchanged Stage 20 engine.
- `refine_projections_stage22.py`: refinement using the unchanged Stage 21 engine.
- `verify_stage22.py`: read-only inspection of saved arrays, integrals, predictions, and certificates.
- `outputs/observed_inputs/`: public inputs supplied to the estimator.
- `outputs/latent_evaluation_only/`: latent samples, including D objects, used only for evaluation.
- `outputs/fits/`: raw EM histories, predictions, and certified intervals for true-parent membership.
- `outputs/refined_projection_certificates/`: final certificates for the three representative samples.

Distinguish membership in the parent confidence set for all 120 samples from the sharp projections computed for only three. The original unconstrained EM preserves the initialized invisible parent mass. Numerical stopping is not interpreted as a global-MLE guarantee.

Read the English [theory](THEORY.md), [summary](SUMMARY.md), and [resumption record](RESUMPTION.md). The resumption record documents additional projection searches performed after the original theory and summary were saved.
