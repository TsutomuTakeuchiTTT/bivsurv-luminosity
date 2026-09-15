# bivsurv-luminosity: first-publication checklist

> English translation of [`docs/RELEASE_CHECKLIST_ja.md`](RELEASE_CHECKLIST_ja.md). The Japanese source is retained unchanged. Execution statements and numerical results refer to the stage or preparation described in that source, not to a new run during this documentation revision. Execution commands are translated from the supplied checklist; no remote commands have been executed as part of this revision.

## Confirmed decisions

The new repository name is `bivsurv-luminosity`, the software licence is `BSD-3-Clause`, and the copyright holder is `Tsutomu T. Takeuchi`. The proposed first public tag is `v0.1.0rc1`, to be marked as a **prerelease** on GitHub. The old private repository's history, visibility, and files must not be changed.

No remote creation, push, GitHub Actions execution, GitHub Release publication, or Zenodo registration was performed when preparing the supplied licensed package. The repository owner and URL must be confirmed through the actual GitHub connection or repository-creation screen. The English-documentation revision does not change this publication status.

## 1. Check the distribution

Extract the ZIP into a new directory and enter the `bivsurv-luminosity` directory containing `LICENSE` and `pyproject.toml`. Publish the extracted source tree, not just a ZIP file placed in a repository.

```bash
python -m pip install -e ".[test]"
python scripts/check_release.py
python -m pytest -q
```

The saved numerical results do not need to be recomputed before they can be inspected. To regenerate the full experiments, use a separate output directory as documented in the [README](../README.md). The main library declares Python 3.10 or newer. The earlier publication preparation was tested with Python 3.13.5. `requirements-pinned.txt` records the environment of the saved experiment and is not necessarily compatible with every other Python version. Checks of this English-documentation package are recorded separately in [documentation verification](english_documentation_verification.json).

## 2. Create a new public repository

Confirm the correct owner in GitHub's new-repository screen. Set the name to `bivsurv-luminosity` and the visibility to **Public**. If a repository with that name already exists, stop and inspect it rather than overwriting it. The description can be:

```text
D-completion bivariate survival estimation for two-band luminosity catalogues, with reproducible benchmarks and validation code.
```

Do not ask GitHub to add a README, `.gitignore`, or licence: these files are already supplied. Commit and push the complete directory with Git or GitHub Desktop. The package contains several thousand files, so do not rely on dragging every file into the browser at once.

The following example uses GitHub CLI. Execute it only in a freshly extracted directory and first verify the authenticated account. If `gh repo create` fails, for example because the repository already exists, do not automatically push to a different remote.

```bash
gh auth status
git init -b main
git add .
git commit -m "Prepare licensed reproducibility release 0.1.0rc1"
gh repo create bivsurv-luminosity --public --source=. --remote=origin --push --description "D-completion bivariate survival estimation for two-band luminosity catalogues, with reproducible benchmarks and validation code."
```

If Git author information has not been configured, set the name and email locally to match the author's GitHub settings before committing. Do not paste authentication tokens or private keys into files that will be committed. These instructions apply only to the new repository; do not reuse the old repository's remote.

## 3. Align metadata with the actual publication

After repository creation, add the verified URL to `repository-code` in [`CITATION.cff`](../CITATION.cff) and to the README where appropriate. Update pre-publication wording to reflect the actual state. Add `date-released` only when the real release date is known, and add a DOI only after it has been issued. Do not invent citation identifiers.

After editing documentation, inspect the intended differences, regenerate the checksums, and commit and push the metadata update.

```bash
python scripts/check_release.py --write
python scripts/check_release.py
python -m pytest -q
git add CITATION.cff README.md README_ja.md SHA256SUMS
git commit -m "Record verified publication metadata"
git push
```

The command above reproduces the original checklist's staging list. Also stage any other documentation files that were intentionally edited. `--write` refreshes integrity information after reviewed changes; it is not a way to suppress unexpected mismatches.

## 4. Check CI, then publish the first prerelease

In GitHub Actions, confirm that **Tests and smoke benchmark** succeeds. The workflow tests Python 3.13 and a small numerical experiment; it is not a rerun of all 480 fits or a new global-optimality proof. Resolve failures before creating a release.

In Releases, use the new tag `v0.1.0rc1`, target `main`, and select **Prerelease**. Use [`RELEASE_NOTES.md`](../RELEASE_NOTES.md) as the release description. Do not move or delete an already public tag of the same name; choose the next version instead. Select a commit containing all required metadata updates.

With GitHub CLI, create and push the local tag at the verified commit before publishing:

```bash
git tag -a v0.1.0rc1 -m "First public release candidate"
git push origin v0.1.0rc1
gh release create v0.1.0rc1 --verify-tag --prerelease --title "bivsurv-luminosity 0.1.0rc1" --notes-file RELEASE_NOTES.md
```

## 5. Add code-availability information to the paper

Add the paper's Code availability statement only after verifying that the repository and release actually exist. If using a persistent archive, record the version-specific DOI once it has been issued. Do not yet insert a statement asserting that the software has been published or assigned a DOI.

## Notes on the reference implementations

`reference_experiments/stage22/resumption_checks/compare_runs.py` is an archival comparison helper that assumes the original working directories. Use `scripts/verify_visual.py` for the new reproducibility comparisons. The fast API remains specialized to aligned finite recording boundaries; it has not become a general API for non-aligned limits or measurement errors. Keep the new experiment's numerical stopping diagnostics separate from the existing Stage 25 global certificates for its different set of 13 candidates.

English [reference-experiment documentation](reference/README.md) explains the saved stages without requiring the Japanese notes. Original Japanese notes and LaTeX fragments remain as historical materials.

## Official procedures cited by the original checklist

- [Create a repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository)
- [GitHub CLI repository creation](https://cli.github.com/manual/gh_repo_create)
- [Citation files](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)
- [GitHub CLI release creation](https://cli.github.com/manual/gh_release_create)
- [Licence text](https://spdx.org/licenses/BSD-3-Clause.html)
