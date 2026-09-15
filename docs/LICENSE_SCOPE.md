# License scope and source provenance

## Approved software license

The BSD 3-Clause license in the root `LICENSE` applies to the original software,
software-use documentation, configuration files, tests, example notebooks, and
accompanying original synthetic benchmark products in this distribution, except
where an explicit separate notice applies. The copyright notice is:

`Copyright (c) 2026, Tsutomu T. Takeuchi`.

The choice was approved by the project author in the current publication workflow.
The previous proposed-license file has been replaced by the operative root license.
This license does not add a requirement to cite the software; `CITATION.cff`
expresses the requested scholarly attribution.

## Scientific manuscript

The manuscript text in `manuscript/` and manuscript `.tex` fragments in the
reference snapshots are provided to document the experiments. They are not
relicensed by this software licensing step. No journal publication agreement or
article reuse license is created or modified here. References to published works
do not grant rights to those works.

## Historical source snapshots and external dependencies

`reference_experiments/` and `validation_reference/` preserve earlier implementation
snapshots from this same project. Directory names such as `vendor`, `vendor19`,
and `prior` refer to the project's earlier stages, as recorded in
`docs/provenance.json`; they are not copies of NumPy or SciPy. The earlier files
and their provenance records were not rewritten when adding the root license.

NumPy, SciPy, Matplotlib, pytest, and build tools are installed separately and
retain their respective licenses. No third-party package binary, font file, or
license override is distributed in this preparation bundle. A filename/header
inspection found no separate third-party copyright/license notice in the source
snapshot. This inspection is not an independent legal provenance audit; any
subsequently identified upstream notice must be preserved rather than replaced.

Canonical license identifier and text:
https://spdx.org/licenses/BSD-3-Clause.html
