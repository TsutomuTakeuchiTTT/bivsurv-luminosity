# bivsurv-luminosity 0.1.0rc1

Draft notes for the first public **prerelease**. This file is not evidence that a
GitHub release has already been published.

## Included

- D-completion bivariate survival estimation for two-band union-selected,
  deterministically interval-recorded catalogues.
- Fixed Gaussian benchmark settings, saved inputs and outputs, figure-generation
  code, and an executed quickstart notebook.
- Archived Stage 22–25 implementation snapshots and their separate certificates.
- BSD-3-Clause license, citation metadata, file manifest, and a small CI workflow.

## Scope

The fast API requires limits aligned with recording edges. The new high-resolution
fits have numerical diagnostics, not global-optimality certificates. The separate
Stage 25 certificates do not certify the new visual fits. The estimator does not
infer an absolute galaxy number density, and the unresolved/invisible regions
are not filled in as recovered density.

## Publication preparation

The author approved the new repository name and BSD-3-Clause licensing. Estimator
code, benchmark configuration, saved numerical results, plots, and the reference
snapshots are unchanged from the earlier local 0.1.0rc1 candidate. The version is
retained because no earlier public 0.1.0rc1 release was created in this workflow.
No repository URL, publication date, DOI, or remote CI result is claimed here.

## English-documentation edition

English results, reproduction and release instructions, a figure guide, and full
Markdown translations of the Stage 22–25 READMEs, theory notes, summaries, and
Stage 22 resumption record are included under `docs/`. Original Japanese notes,
scientific code, configuration, stored outputs, and certificates are unchanged.
The README and quickstart now link to the English documentation. A documentation-only
checker and refreshed root manifest are included. The full Japanese manuscript is
not translated by this software-documentation edition. No remote publication
operation has been performed. See `docs/english_documentation_verification.json`
for the checks actually executed for this edition.
