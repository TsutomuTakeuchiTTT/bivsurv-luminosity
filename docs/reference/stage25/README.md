# Stage 25: certificates of the global likelihood gap

> English translation of [`reference_experiments/stage25/README_ja.md`](../../../reference_experiments/stage25/README_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision. The disposable-copy instruction is an added packaging safeguard. Refer to the root reproducibility guide for a read-only check of saved certificates.

The English [theory](THEORY.md) gives the derivation and conditions; the [summary](SUMMARY.md) gives the results. The archived `stage24_record.tex` is the Japanese LaTeX insertion recording the preceding stage.

## Rerunning

Python 3, NumPy, and SciPy are used. All inputs are included. Run these commands in `reference_experiments/stage25/` **in a disposable copy**: the paths below are the original output paths, and some of the commands recreate saved certificates.

```bash
OPENBLAS_NUM_THREADS=1 python global_bound.py --tol 1/1000000 --max-nodes 501 --outdir outputs/main_certificates
OPENBLAS_NUM_THREADS=1 python global_bound.py --tol 1/1000000 --max-nodes 501 --outdir outputs/rerun_certificates
OPENBLAS_NUM_THREADS=1 python global_bound.py --case coarse_d3__lambda0_N400_r004 --tol 1/1000000 --max-nodes 1 --outdir outputs/budget_control
OPENBLAS_NUM_THREADS=1 python global_bound.py --case coarse_d3__lambda0_N400_r004 --tol 1/100000000 --max-nodes 501 --outdir outputs/strict_control
OPENBLAS_NUM_THREADS=1 python verify_stage25.py
OPENBLAS_NUM_THREADS=1 python audit_stage25.py
```

All 13 main cases attained a mean-log-likelihood gap of at most 1e-6. Failure to meet the requested accuracy is the correct saved result for the budget control and stricter-accuracy control. Public inputs and snapshots of the original candidates are in `inputs/`. Ancestor code, original inputs, and original certificates are in `reference_stage24/`. The new search does not replace the original candidate masses.

## Main outputs

- `outputs/main_certificates/`: global upper bounds and covering trees for all 13 cases.
- `outputs/global_gap_comparison.csv`: first-order residuals, global gaps, and gaps after 80 GEM updates.
- `outputs/independent_verification.json`: verification through a separate representation, without optimization.
- `outputs/audit_results.json`: agreement of the repeated searches for all cases and arithmetic controls.
- `outputs/same_likelihood_family.json`: three parent distributions with the same likelihood.
- `outputs/stage24_recheck.json`: rechecking of the preceding stage's 1040 updates.

The rational-logarithm implementation `reference_stage24/vendor/exact.py` and observation-geometry code are reused unchanged. The independent verifier does not call the new search's `tangent_coeff`, `polytope`, LP, or SLSQP. It does, however, share the foundational rational-logarithm implementation and the definition of observation geometry; those were not independently reinvented. Altered inputs or certificates can cause verification to fail. The result is not certification of exact equality to an MLE or a statistical error of 1e-6.
