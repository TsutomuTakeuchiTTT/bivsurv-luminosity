# Documentation language and preservation policy

The publication-facing documentation is in English. Japanese source notes remain as optional historical records. This is an additive documentation revision: the English notes are outside the archived experiment directories, and no original numerical evidence is relabelled as a new run.

## Scope of translation

Every Japanese Markdown document in the supplied licensed package has an English counterpart. The results note, release checklist, and thirteen Stage 22–25 Markdown notes have complete translations, including their equations, commands, tables, qualifications, and historical execution scope. The root Japanese README already has a more extensive English counterpart; its navigation has been updated.

The six Japanese LaTeX fragments are retained unchanged as manuscript/work records. This is not a translation of the complete scientific article, which remains a separate task after the Japanese content is fixed. The English documentation provides the required software and reproducibility reading path without those fragments. The main and archived Python sources and machine-readable result records already use English identifiers and text; they have not been rewritten for this edition.

## Source-to-English map

Each English translation links to its unchanged original and is listed below. SHA-256 source hashes and the same mapping are recorded in [`english_documentation_map.json`](english_documentation_map.json).

| Original Japanese Markdown | English counterpart |
|---|---|
| [`docs/RESULTS_ja.md`](RESULTS_ja.md) | [`docs/RESULTS.md`](RESULTS.md) |
| [`docs/RELEASE_CHECKLIST_ja.md`](RELEASE_CHECKLIST_ja.md) | [`docs/RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md) |
| [`reference_experiments/stage22/README_ja.md`](../reference_experiments/stage22/README_ja.md) | [`docs/reference/stage22/README.md`](reference/stage22/README.md) |
| [`reference_experiments/stage22/RESUMPTION_ja.md`](../reference_experiments/stage22/RESUMPTION_ja.md) | [`docs/reference/stage22/RESUMPTION.md`](reference/stage22/RESUMPTION.md) |
| [`reference_experiments/stage22/SUMMARY_ja.md`](../reference_experiments/stage22/SUMMARY_ja.md) | [`docs/reference/stage22/SUMMARY.md`](reference/stage22/SUMMARY.md) |
| [`reference_experiments/stage22/THEORY_ja.md`](../reference_experiments/stage22/THEORY_ja.md) | [`docs/reference/stage22/THEORY.md`](reference/stage22/THEORY.md) |
| [`reference_experiments/stage23/README_ja.md`](../reference_experiments/stage23/README_ja.md) | [`docs/reference/stage23/README.md`](reference/stage23/README.md) |
| [`reference_experiments/stage23/SUMMARY_ja.md`](../reference_experiments/stage23/SUMMARY_ja.md) | [`docs/reference/stage23/SUMMARY.md`](reference/stage23/SUMMARY.md) |
| [`reference_experiments/stage23/THEORY_ja.md`](../reference_experiments/stage23/THEORY_ja.md) | [`docs/reference/stage23/THEORY.md`](reference/stage23/THEORY.md) |
| [`reference_experiments/stage24/README_ja.md`](../reference_experiments/stage24/README_ja.md) | [`docs/reference/stage24/README.md`](reference/stage24/README.md) |
| [`reference_experiments/stage24/SUMMARY_ja.md`](../reference_experiments/stage24/SUMMARY_ja.md) | [`docs/reference/stage24/SUMMARY.md`](reference/stage24/SUMMARY.md) |
| [`reference_experiments/stage24/THEORY_ja.md`](../reference_experiments/stage24/THEORY_ja.md) | [`docs/reference/stage24/THEORY.md`](reference/stage24/THEORY.md) |
| [`reference_experiments/stage25/README_ja.md`](../reference_experiments/stage25/README_ja.md) | [`docs/reference/stage25/README.md`](reference/stage25/README.md) |
| [`reference_experiments/stage25/SUMMARY_ja.md`](../reference_experiments/stage25/SUMMARY_ja.md) | [`docs/reference/stage25/SUMMARY.md`](reference/stage25/SUMMARY.md) |
| [`reference_experiments/stage25/THEORY_ja.md`](../reference_experiments/stage25/THEORY_ja.md) | [`docs/reference/stage25/THEORY.md`](reference/stage25/THEORY.md) |
| [`README_ja.md`](../README_ja.md) | [`README.md`](../README.md) |

## Explicit editorial additions

The English documentation index, reproducibility guide, figure-file map, and stage index organize information already present in the supplied scripts, configuration, notebook, and original notes. They are guides, not additional scientific results. Translation banners state that executions belong to the historical stage described. Stage 22's later resumption remains a separate record; Stage 24's first-order-only statement is not silently rewritten using Stage 25's subsequent result.

The English results note retains its original “awaiting licence confirmation” publication status as historical wording, followed by the approved status documented in the licensed package. Commands requiring regeneration are marked for a disposable copy or a separate output directory. In particular, the supplied visual verifier writes `verification.json` under its `--output` path; the main English comparison example now puts the rerun tree first to avoid rewriting the archived report. No verifier or estimator implementation was changed to make that adjustment.

The `quickstart.ipynb` edit changes Markdown instructions and links only. Original executable cells, stored outputs, execution counts, and metadata are retained. The BSD-3-Clause licence, citation metadata, version 0.1.0rc1, original provenance, and prior validation reports remain unchanged. The refreshed root `SHA256SUMS` describes this documentation edition; older stage manifests and provenance remain untouched.

## Checks and limits of automated checks

[`english_documentation_verification.json`](english_documentation_verification.json) records the checks executed for this edition, separately from earlier experiment and publication-preparation reports. [`scripts/check_documentation.py`](../scripts/check_documentation.py) audits the translation inventory, original-source hashes, local Markdown destinations, notebook Markdown destinations, and unexpected Japanese/CJK text outside preserved sources. Such checks detect omissions and broken references; they do not prove semantic equivalence of translations or new statistical guarantees. The English translations have also been reviewed against their source notes.

No GitHub repository, tag, release, DOI, publication date, or remote CI success is created by this documentation revision. Those steps follow the [English release checklist](RELEASE_CHECKLIST.md).
