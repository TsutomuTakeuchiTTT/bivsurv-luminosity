# Reproducibility guide

This guide organizes the commands already supplied by the package and its stage READMEs. It distinguishes **file-integrity checks**, **refitting public records**, **regeneration of simulations and figures**, and **checking archived mathematical certificates**. A success in one category is not a success in every category.

## 1. Inspect and test the distributed snapshot

From the repository root:

```bash
python -m pip install -e ".[test]"
python scripts/check_release.py
python scripts/check_documentation.py
python -m pytest -q
```

`check_release.py` checks the manifest and package metadata. `check_documentation.py` checks the English-documentation inventory, source-translation references, local Markdown links, and language boundaries; it does not prove translation equivalence or statistical correctness. The existing unit tests cover the implemented numerical and observation conventions. These commands do not regenerate archived benchmark outputs.

Dependencies must already be installed or made available by the user's package-management environment. The library declares Python >=3.10. The earlier benchmark environment is recorded in [`requirements-pinned.txt`](../requirements-pinned.txt); it is not a guarantee of identical results on all versions or platforms. This English edition's executed checks and runtime are recorded in [the new verification report](english_documentation_verification.json).

## 2. Refit one catalogue without its latent truth

The following command reads only the geometry, public counts, and initial mass used by the estimator:

```bash
python scripts/fit_counts.py --model outputs/visual/models/deep.json --input outputs/visual/public_counts/deep/rho0.55_N16000_r000.json --initial outputs/visual/models/deep_arrays.npz --output outputs/my_fit.json
```

[`quickstart.ipynb`](../quickstart.ipynb) demonstrates the same public-count route and compares the fitted mass to the saved result after fitting. Latent arrays, true Gaussian parameters, discarded D values, and the parent total are not inputs to `fit_counts.py`. The full record alphabet, including zero-count categories, remains in the model.

## 3. Regenerate the visual experiment in a separate directory

Do not overwrite the distributed `outputs/visual/` tree when creating a fresh comparison:

```bash
python scripts/run_benchmark.py --output outputs/rerun
python scripts/verify_visual.py --output outputs/rerun
python scripts/make_figures.py --output outputs/rerun
python scripts/check_initialization.py --output outputs/rerun
```

The full configuration comprises 120 independent parent catalogues, 240 paired design evaluations, and 480 EM fits. The initialization command adds 12 separate sensitivity fits. Figures use the predeclared representative parent size 16000 and replicate 0, so a smoke-only run is not enough to generate the complete figure set.

For a smaller generation-and-fitting smoke test:

```bash
python scripts/run_benchmark.py --replicates 1 --sizes 1000 --output outputs/smoke
python scripts/verify_visual.py --output outputs/smoke
```

This gives four observation-design/correlation evaluations, not a rerun of the full study.

## 4. Compare a full rerun without modifying archived verification reports

The supplied command `python scripts/verify_visual.py --compare-rerun outputs/rerun` uses `outputs/visual/` as the reference but also rewrites its `verification.json`. For an immutable distribution, perform that comparison in a disposable copy. Alternatively, the same existing verifier can use the fresh rerun as the output tree and the saved tree as the comparison tree:

```bash
python scripts/verify_visual.py --output outputs/rerun --compare-rerun outputs/visual
```

The latter invocation writes its report under `outputs/rerun/`, leaving the distributed results untouched. This is an operational clarification based on the supplied script's output behaviour, not a change to the estimator or comparison implementation. Exact array/JSON agreement is the criterion implemented by that script; platform-dependent differences should be investigated rather than silently changing scientific tolerances or refreshing saved hashes.

## 5. Check older certificates without rerunning the searches

Use a disposable copy for these commands because the verifiers write summary reports. The Stage 25 verifier also uses a temporary tamper-test file under its own output directory. The commands below are existing entry points; they do not invoke the corresponding numerical fitting or global-bound searches.

```bash
cd reference_experiments/stage22
python verify_stage22.py --outdir outputs
cd ../stage23
python verify_stage23.py --outdir outputs
cd ../stage24
python verify_stage24.py --outdir outputs
cd ../stage25
python verify_stage25.py
```

The stage READMEs document their respective dependencies; Stage 23 additionally records SymPy. These independent verifiers may reconstruct record geometry, integrate the declared truth, or recompute the finite learning-prediction definition. “No numerical refit” does not mean no arithmetic is performed. Stage-specific qualifications and counts are in the [English reference notes](reference/README.md).

For full regeneration of each older stage, follow that stage's translated README in a disposable copy. Do not assume that all nested snapshots are called by the main visual-benchmark driver. The archival Stage 22 `compare_runs.py` helper retains old local paths and is not a portable release entry point.

## 6. Preserve the distinction between checks

The high-resolution visual experiment has floating-point checks of records, updates, likelihood, scores, and diagnostics. The earlier finite-model stages have rational certificates for their explicitly stated tasks. None of those checks, and no file hash, turns a representative parent curve into a uniquely identified distribution or an empirical replicate band into a confidence band.

Run `check_release.py --write` only after intended, reviewed edits to the release contents. Its original exclusions cover `outputs/rerun/`, `outputs/smoke/`, caches, build products, and the example `outputs/my_fit.json`; arbitrary new output paths may be reported as unlisted. That is why independent-output runs or disposable copies are recommended. Do not update the manifest simply to conceal an unexplained mismatch.
