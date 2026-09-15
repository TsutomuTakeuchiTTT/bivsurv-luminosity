# Bivariate survival estimation for two-band luminosity catalogues

**bivsurv-luminosity 0.1.0rc1 — BSD-3-Clause.**

This is the licensed publication-preparation copy. Creating the remote repository
and publishing its first prerelease remain separate steps; see
[`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md).

This bundle adds a reproducible, data-driven luminosity-function visualization
benchmark to the D-completion methodology. It also preserves the previous
Stage 22–25 implementations, inputs and validation results in
`reference_experiments/`. These are not the rejected product-integral/weighting
implementation from the old repository.

## English documentation

Start with the [documentation index](docs/README.md). English [results](docs/RESULTS.md),
[figure interpretation](docs/FIGURES.md), [reproduction instructions](docs/REPRODUCIBILITY.md),
and [Stage 22–25 theory and result notes](docs/reference/README.md) are supplied.
The original Japanese notes remain optional provenance records; their English
counterparts are outside the archived snapshots so no original evidence is changed.
The complete scientific manuscript is not translated by this software-documentation edition.
See the [language policy and source map](docs/DOCUMENTATION_POLICY.md).

## Quick start

From this directory, using Python 3.10 or newer:

```bash
python -m pip install -e ".[test]"
python scripts/check_release.py
python scripts/check_documentation.py
python -m pytest -q
```

The commands above check the distributed snapshot and unit tests without
regenerating the archived numerical outputs. The full benchmark, validation,
and figures can be regenerated in a separate directory:

```bash
python scripts/run_benchmark.py --output outputs/rerun
python scripts/verify_visual.py --output outputs/rerun
python scripts/make_figures.py --output outputs/rerun
```

The full new benchmark consists of 120 independent parent catalogues, paired
across two observation designs (240 design evaluations). Each has a
D-completion fit, a deliberately misspecified no-D ablation, and an uncorrected
A-only histogram. There are 480 EM fits, plus 12 optional initialization checks.
All default figure data and results are already in `outputs/visual/`.

A smaller numerical smoke test is:

```bash
python scripts/run_benchmark.py --replicates 1 --sizes 1000 --output outputs/smoke
python scripts/verify_visual.py --output outputs/smoke
```

The figure command uses the predeclared N=16000, replicate=0 catalogue; run the
full configuration before making the full figure set. No successful catalogue
is selected after looking at performance. Jupyter entry points are illustrated
in `quickstart.ipynb`.

## What is estimated and what the plots mean

`X_b = log10(L_b/L_b,0)`, `Z_b = -X_b`. Detection means `Z_b <= C_b`
(including equality). A, B and C are retained; D is absent. Detected luminosities
are interval records, not noiseless point-density observations. The likelihood
is `prod P(E_i)/P(D_i^c)`. Its D-completion update is implemented without
post-hoc clipping, CDF repair, smoothing or a probability floor.

Finite bins and both overflow intervals cover **the whole real plane**. The
Gaussian truth is used solely by the sampler and evaluator, never to delete
candidate cells. The benchmark grid and reference distribution are specified
before sampling. The reference is independent logistic in the declared
centered/scaled plotting coordinates; it shares those predeclared coordinate
scales with the mock, but does not use the true correlation or impose a
Gaussian parent model. Shifted reference distributions are tested separately.

Because permanent non-detection prevents full-parent normalization from being
identified, the primary joint maps show the observable-normalized distribution
`Q = P(. | U)`, where U is the union of the declared detection regions. The
Gaussian truth uses exactly the same conditioning. The maps mask regions where
at least one luminosity remains unresolved even at the deepest limit; they do
not extrapolate a density through those regions. In resolved bins the displayed
quantity is `Q(bin)/(bin width in dex)^2`, not a fitted continuous density.
All overflow and unresolved masses remain in the numerical model and its
normalization.

The primary error is RMSE of `Q(X1 >= x1, X2 >= x2)` on a fixed 13x13 threshold
grid in standardized luminosities [-0.5, 2.5]. The additional resolved-bin L1
error and record-law TV error are also saved. The full-parent reference error
is reported separately and depends on invisible initialization. Marginal plots
condition on BOTH standardized luminosities >= -0.5; they are not extrapolated
full-parent marginals.

The `same_likelihood` figures show three parent distributions with the same
fitted conditional record law. **They are not confidence bands.** The shaded
regions on repeat-error figures are empirical 10th–90th percentiles of 20
realizations, **not** 95% confidence intervals. The raw parent number is used
for designing/evaluating simulations, never as an estimator input.

## Comparisons and their interpretation

* `D_completion`: the manuscript's nonnegative mass update, with numerical
  full-column score diagnostics.
* `no_D_ablation`: retains A/B/C censoring information but omits the selection
  denominator and missing-D completion. A deliberate ablation, not a valid
  competing method for the same selection experiment.
* `A_only_histogram`: drops B/C and performs no truncation correction. This is
  **not** Lynden–Bell, 1/Vmax, or the correctly truncated A-only NPMLE.

No claim of superiority to all historical estimators is made. These controls
isolate the information loss from deleting limits and omitting D correction.

## Algorithm scope and numerical guarantees

The new fast `GridModel` specializes the exact interval representation to
thresholds aligned with recording edges. It **rejects** non-aligned thresholds
rather than silently quantizing them. The full arbitrary-interval geometry
remains in the recovered `reference_experiments` code. A singleton at a
recording edge has the same observation signature as the left adjacent
Z interval under the stated endpoint convention; this equivalence is tested
against the earlier full singleton/open-cell implementation.

The new EM routine maximizes over the simplex. A supplied tail bound is
**diagnosed, not imposed**. All 240 new D-completion candidates passed the
declared half-tail diagnostic in this experiment. For an input that fails,
do not call it a tail-constrained estimate: use the constrained M-step in the
reference implementation or adapt the general constrained formulation.

The new high-resolution fits have numerical score/likelihood checks;
**global-optimality certificates were not computed for these fits**. This
must not be confused with Stage 25's certificates for its separate 13 cases.
The preserved Stage 25 proof verifier was rerun successfully for those cases.

Truth bin masses use deterministic conditional-Gaussian one-dimensional
quadrature, including infinite tails. Its numerical error estimates and
tightening/marginal checks are recorded; these are not rigorous interval bounds.

## Reproduce one fit without access to the latent truth

```bash
python scripts/fit_counts.py --model outputs/visual/models/deep.json --input outputs/visual/public_counts/deep/rho0.55_N16000_r000.json --initial outputs/visual/models/deep_arrays.npz --output outputs/my_fit.json
```

`fit_counts.py` receives no latent array, population parameter, parent count,
non-detection value, or truth mass. `latent_evaluation_only/` is reserved for
simulation checks. Model JSON files contain the complete possible-record map;
zero-count record categories are retained.

## Reproducibility and provenance

Run in a second output directory and compare:

```bash
python scripts/run_benchmark.py --output outputs/rerun
python scripts/verify_visual.py --output outputs/rerun --compare-rerun outputs/visual
```

The comparison above writes its new report to `outputs/rerun/`, not to the archived
`outputs/visual/`. This is an invocation-only documentation adjustment; the verifier
implementation is unchanged.

Numerical arrays and public record JSONs matched exactly in the tested
environment; metrics matched apart from wall-clock time. Portable reproducibility
on another NumPy/SciPy/platform combination should be judged with suitable
numerical tolerances, not promised byte identity for PDFs or ZIP timestamps.

`validation_reference/` contains unmodified earlier code used in the unit
cross-check. `reference_experiments/` contains the complete recovered Stage
22–25 packages (including their nested dependency snapshots). File hashes and
changes to packaging are recorded in `docs/provenance.json`. New results do not
replace the earlier paper's experiments or their archived evidence.

## Licensing and publication status

The repository name **bivsurv-luminosity** and the **BSD-3-Clause** license have
been approved. The complete license is in [`LICENSE`](LICENSE); the scope and
preserved-reference provenance are described in
[`docs/LICENSE_SCOPE.md`](docs/LICENSE_SCOPE.md).

The release-candidate version remains **0.1.0rc1** because this preparation
has not yet been publicly released in this workflow. The latest revision adds
English documentation and updates its reading links only.
No estimator, benchmark settings, saved results, figures, or archived reference
implementation have been changed. The scientific manuscript is not relicensed
by this software-license decision.

Use [`CITATION.cff`](CITATION.cff) for software citation metadata. A repository
URL, release date, and archival DOI must be added only after they are verified.
None has been invented for this local copy. See
[`RELEASE_NOTES.md`](RELEASE_NOTES.md) for the proposed first-prerelease notes.
The old private repository is not reused or modified.

A CI workflow is included for unit tests and a small smoke experiment on Python
3.13. Its GitHub-hosted execution is pending the initial push. The earlier local preparation
checks remain in `docs/publication_preparation_verification.json`; checks run for this
English-documentation revision are separate in
[`docs/english_documentation_verification.json`](docs/english_documentation_verification.json).

## Snapshot integrity

`python scripts/check_release.py` checks every path listed in `SHA256SUMS`, the
approved license, and version agreement. It does not execute the estimator or
claim new statistical validation. After intentionally editing metadata before a
release, refresh the manifest with `python scripts/check_release.py --write`,
then rerun the check. Do not refresh hashes merely to hide an unexpected mismatch.
A Python wheel contains the installable library; use this complete repository
bundle for the benchmark data, notebooks, and earlier certificates.
