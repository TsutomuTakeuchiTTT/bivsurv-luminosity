#!/usr/bin/env python3
"""Check this release snapshot. Use --write only after intentional changes.

Standard library only (Python >=3.10). This is a packaging/integrity check,
not a scientific-validation or global-optimality certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
SKIP_PARTS = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", "build", "dist"}
SKIP_OUTPUTS = {"outputs/rerun", "outputs/smoke"}


def included_files(root: Path) -> list[Path]:
    result = []
    for path in root.rglob("*"):
        rel = path.relative_to(root)
        name = rel.as_posix()
        if any(p in SKIP_PARTS or p.endswith(".egg-info") for p in rel.parts):
            continue
        if any(name == p or name.startswith(p + "/") for p in SKIP_OUTPUTS):
            continue
        if name in {"SHA256SUMS", "outputs/my_fit.json", ".DS_Store"}:
            continue
        if path.suffix == ".pyc":
            continue
        if path.is_symlink():
            raise ValueError(f"Symlink must be reviewed before packaging: {name}")
        if path.is_file():
            if path.name == ".env" or path.name.startswith(".env."):
                raise ValueError(f"Environment file must not be packaged: {name}")
            result.append(path)
    return sorted(result, key=lambda p: p.relative_to(root).as_posix())


def metadata(root: Path) -> dict:
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    cff = (root / "CITATION.cff").read_text(encoding="utf-8")
    init = (root / "bivsurv/__init__.py").read_text(encoding="utf-8")
    patterns = [(text, r'^version\s*=\s*"([^"\n]+)"'),
                (cff, r'^version:\s*[\'"]?([^\'"\s]+)'),
                (init, r"^__version__\s*=\s*['\"]([^'\"]+)['\"]")]
    versions = []
    for source, pattern in patterns:
        match = re.search(pattern, source, re.MULTILINE)
        if match is None:
            raise ValueError("Missing version metadata")
        versions.append(match.group(1))
    if len(set(versions)) != 1:
        raise ValueError(f"Version mismatch: {versions}")
    if not re.search(r'^license\s*=\s*"BSD-3-Clause"$', text, re.MULTILINE):
        raise ValueError("pyproject license must be BSD-3-Clause")
    if not re.search(r'^license:\s*BSD-3-Clause$', cff, re.MULTILINE):
        raise ValueError("CFF license must be BSD-3-Clause")
    license_text = (root / "LICENSE").read_text(encoding="utf-8")
    required = ["BSD 3-Clause License", "Copyright (c) 2026, Tsutomu T. Takeuchi",
                "1. Redistributions of source code", "2. Redistributions in binary form",
                "3. Neither the name", '"AS IS"']
    if not all(part in license_text for part in required):
        raise ValueError("Incomplete license text")
    if (root / "docs/LICENSE_BSD3_PROPOSED.txt").exists():
        raise ValueError("Obsolete license proposal remains")
    return {"version": versions[0], "license": "BSD-3-Clause"}


def check(root: Path) -> dict:
    info = metadata(root)
    failures = []
    seen = set()
    checked = 0
    for number, line in enumerate((root / "SHA256SUMS").read_text(encoding="utf-8").splitlines(), 1):
        if not line:
            continue
        try:
            expected, name = line.split("  ", 1)
        except ValueError as exc:
            raise ValueError(f"Malformed manifest line {number}") from exc
        parts = PurePosixPath(name)
        if not re.fullmatch(r"[0-9a-f]{64}", expected) or parts.is_absolute() or ".." in parts.parts or "\\" in name:
            raise ValueError(f"Unsafe or invalid manifest line {number}")
        if name in seen:
            raise ValueError(f"Duplicate manifest path: {name}")
        seen.add(name)
        path = root / parts
        if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
            failures.append(name + ": missing or unsafe")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            failures.append(name + ": SHA-256 mismatch")
        checked += 1
    unlisted = [p.relative_to(root).as_posix() for p in included_files(root)
                if p.relative_to(root).as_posix() not in seen]
    if failures or unlisted:
        raise ValueError(json.dumps({"failures": failures, "unlisted_files": unlisted}, indent=2))
    return {"passed": True, "files_checked": checked, **info,
            "scope": "packaging and file integrity only"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="replace hashes after reviewed, intentional edits")
    args = parser.parse_args()
    if args.write:
        metadata(ROOT)
        lines = [hashlib.sha256(p.read_bytes()).hexdigest() + "  " + p.relative_to(ROOT).as_posix()
                 for p in included_files(ROOT)]
        (ROOT / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
    try:
        report = check(ROOT)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Release check failed: {exc}\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
