# The new bivariate luminosity visualization experiment

> English translation of [`docs/RESULTS_ja.md`](RESULTS_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision. The final section distinguishes the original pre-approval status from the later licensed preparation; it does not rewrite the historical results.

## Scope of the work performed

New Gaussian log-luminosity samples were generated with means (10.5, 10.2), standard deviations (0.55, 0.65), and correlations 0.55 and 0.90. Parent sample sizes were 1000, 4000, and 16000, with 20 replicates for every combination. There were 120 independent parent samples and 840000 parent objects. The same latent values, observation-layer assignments, and training assignments were used in the shallow and deep designs, giving 240 design evaluations. The calculations comprised 480 fits with or without D-completion, 240 uncorrected A-only histograms, and 12 additional initialization fits.

## Main results

[`summary.csv`](../outputs/visual/summary.csv) contains means and sample standard deviations calculated after evaluating the errors separately for each replicate. The D-completion RMSE for observable-normalized integrated probabilities decreased from approximately 0.0086–0.0111 with 1000 parent objects to approximately 0.0024–0.0030 with 16000 parent objects. The no-D and simple A-only alternatives deliberately omit parts of the observation mechanism; they are not comparisons against correctly formulated historical estimators.

The main density maps show bin-average densities of `Q_U = P(. | U)`. Regions in which both luminosities cannot be resolved simultaneously are masked and are not filled by extrapolation. Invisible mass retains its initial value, so this result must be distinguished from recovery of the full parent distribution. Bands in the replicate plots are empirical 10th–90th percentiles, not confidence bands. Curves for the three equal-likelihood candidates are not confidence bands either.

## Checks recorded when the visual benchmark was prepared

All 11 unit tests passed, including comparisons with the earlier rational update and checks of the observation map at boundaries. All 240 public inputs were reconstructed from the latent arrays, and all initial updates were compared with a dense-matrix implementation. The complete new experiment was rerun in a separate output directory: 240 input JSON files and 1808 arrays agreed, and the metric tables agreed apart from execution time. These results are saved in [`verification.json`](../outputs/visual/verification.json). Every new candidate satisfied the numerical score condition; no global-optimality proof was computed for the new 1156-mass candidates.

The independent Stage 25 verifier was also rerun for its earlier 13 candidates. Those are saved certificates for a different, smaller model; they do not confer global guarantees on the new high-resolution figures. The original Stage 22–25 implementations were retained in [`reference_experiments/`](../reference_experiments/), but not all of their earlier sample calculations were rerun during the visual-benchmark preparation.

## Figures and manuscript material

All 41 separate plots were saved as PDF and PNG. Fifteen of them were assembled into six manuscript figures, accompanied by a 12-row table. Other manuscript text, appendices, bibliography files, and settings were left unchanged. The added Japanese passage is [`manuscript/visual_benchmark_ja.tex`](../manuscript/visual_benchmark_ja.tex). The restriction of the fast API to limits aligned with recording boundaries is stated in the root README.

## Publication status recorded in the original results note

At the time this results note was written, the package had not been published on GitHub. It was a candidate awaiting confirmation of the repository name and licence. The state of the old repository had not been changed.

**Subsequent preparation status:** the name `bivsurv-luminosity` and BSD-3-Clause licence were approved in the supplied licensed package. This documentation revision does not itself publish a repository or release. The current next steps are in the English [release checklist](RELEASE_CHECKLIST.md).
