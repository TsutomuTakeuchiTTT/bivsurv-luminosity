#!/usr/bin/env python3
"""Audit English documentation and preserved-source links (standard library only).

This checks inventory, source hashes, language markers, and local link targets.
It does not certify translation semantics, external web links, or science results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
CJK = re.compile(r"[\u3040-\u30ff\u3400-\u9fff]")
SKIP = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", "build", "dist"}
INLINE_LINK = re.compile(r"!?\[[^\]\n]*\]\(([^\s)]+)(?:\s+\"[^\"]*\")?\)")


def retained_file(path: Path) -> bool:
    rel = path.relative_to(ROOT)
    if any(part in SKIP or part.endswith('.egg-info') for part in rel.parts):
        return False
    return not (rel.as_posix().startswith(('outputs/rerun/', 'outputs/smoke/')))


def local_target_errors(path: Path, text: str) -> tuple[list[str], int]:
    """Check local paths, not URL availability or Markdown anchor generation."""
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    errors: list[str] = []
    checked = 0
    for match in INLINE_LINK.finditer(text):
        target = match.group(1).strip('<>')
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc or not parsed.path:
            continue
        dest = (path.parent / unquote(parsed.path)).resolve()
        checked += 1
        if not dest.is_relative_to(ROOT.resolve()):
            errors.append(f"{path.relative_to(ROOT)}: target outside root: {target}")
        elif not dest.exists():
            errors.append(f"{path.relative_to(ROOT)}: missing target: {target}")
    return errors, checked


def check() -> dict:
    data = json.loads((ROOT / 'docs/english_documentation_map.json').read_text(encoding='utf-8'))
    errors: list[str] = []
    sources = set()
    translations = set()
    for row in data['translations']:
        src, dst = ROOT / row['source'], ROOT / row['english']
        sources.add(src.resolve())
        translations.add(dst.resolve())
        if not src.is_file():
            errors.append(f"Missing Japanese source: {row['source']}")
        elif hashlib.sha256(src.read_bytes()).hexdigest() != row['source_sha256']:
            errors.append(f"Original source hash changed: {row['source']}")
        if not dst.is_file():
            errors.append(f"Missing English counterpart: {row['english']}")
        elif CJK.search(dst.read_text(encoding='utf-8')):
            errors.append(f"Unexpected Japanese/CJK characters in English counterpart: {row['english']}")
    preserved_tex = { (ROOT / p).resolve() for p in data['preserved_japanese_latex'] }
    for path in preserved_tex:
        if not path.is_file():
            errors.append(f'Missing preserved manuscript fragment: {path}')
    english_md = []
    original_md = []
    links = 0
    for path in sorted(ROOT.rglob('*.md')):
        if not retained_file(path):
            continue
        text = path.read_text(encoding='utf-8')
        if CJK.search(text):
            original_md.append(path)
            if path.resolve() not in sources:
                errors.append(f'Japanese Markdown without mapped English counterpart: {path.relative_to(ROOT)}')
        else:
            english_md.append(path)
            found, count = local_target_errors(path, text)
            errors.extend(found)
            links += count
    unexpected_japanese = []
    scan_count = 0
    extensions = {'.py', '.json', '.ipynb', '.tex', '.cff', '.toml', '.yml', '.yaml', '.txt'}
    for path in ROOT.rglob('*'):
        if not path.is_file() or not retained_file(path) or path.suffix not in extensions:
            continue
        scan_count += 1
        text = path.read_text(encoding='utf-8', errors='strict')
        if path.suffix in {'.json', '.ipynb'}:
            text = json.dumps(json.loads(text), ensure_ascii=False)
        if CJK.search(text) and path.resolve() not in preserved_tex:
            unexpected_japanese.append(path.relative_to(ROOT).as_posix())
    errors.extend('Unexpected Japanese in public text: '+s for s in unexpected_japanese)
    notebook = json.loads((ROOT/'quickstart.ipynb').read_text(encoding='utf-8'))
    markdown_cells = 0
    for cell in notebook['cells']:
        if cell['cell_type'] == 'markdown':
            text = ''.join(cell['source'])
            found, count = local_target_errors(ROOT/'quickstart.ipynb', text)
            errors.extend(found)
            links += count
            markdown_cells += 1
    report = {'passed':not errors, 'english_markdown_files':len(english_md),
              'preserved_japanese_markdown_with_english_counterparts':len(original_md),
              'mapped_japanese_markdown_sources':len(sources), 'english_counterparts':len(translations),
              'preserved_japanese_latex_fragments':len(preserved_tex),
              'local_link_destinations_checked':links, 'notebook_markdown_cells_checked':markdown_cells,
              'additional_text_files_scanned':scan_count, 'errors':errors,
              'scope':'Documentation inventory, original-source hashes, local link paths and language scan. '
                      'Not semantic translation certification, external URL checking, or scientific validation.'}
    return report


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Optional JSON report; prefer a path outside the release snapshot.')
    args=parser.parse_args()
    try:
        report=check()
    except (OSError,ValueError,KeyError) as exc:
        parser.exit(1,f'Documentation check failed: {exc}\n')
    encoded=json.dumps(report,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(encoded,encoding='utf-8')
    print(encoded,end='')
    if not report['passed']:
        raise SystemExit(1)

if __name__=='__main__':
    main()
