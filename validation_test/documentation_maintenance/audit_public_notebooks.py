"""Read-only corpus checks; optional generated JSON evidence, never rerun claims."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'packages/radia-mcp/src'))
from radia_mcp.bibliography._bibparse import read_bib_file
from radia_mcp.document_meta.tools import (
    document_meta_notebook_citation_audit,
    document_meta_notebook_result_audit,
)


def audit(live_bibliography: Path | None = None) -> dict:
    citations = document_meta_notebook_citation_audit(str(ROOT), max_items=10000)
    results = document_meta_notebook_result_audit(str(ROOT), max_items=10000)
    rows = []
    cited = set()
    for result in results['notebooks']:
        rel = result['path']
        path = ROOT / rel
        nb = json.loads(path.read_bytes())
        old_bytes = subprocess.check_output(['git', 'show', 'HEAD:' + rel], cwd=ROOT)
        old = json.loads(old_bytes)
        code = lambda n: [c for c in n['cells'] if c['cell_type'] == 'code']
        meta = nb.get('metadata', {}).get('radia', {})
        cited.update(meta.get('citation_keys', []))
        code_cells = code(nb)
        rich = [i for i, c in enumerate(code_cells) if any(
            'application/vnd.webgui' in str(o.get('data', {}))
            or 'webgui' in str(o.get('data', {})).lower()
            for o in c.get('outputs', []))]
        rows.append({
            'path': rel,
            'notebook_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'baseline_notebook_sha256': hashlib.sha256(old_bytes).hexdigest(),
            'complete_code_cells_equal_head': code(nb) == code(old),
            'citation_keys': meta.get('citation_keys', []),
            'notebook_role': meta.get('notebook_role'),
            'webgui_required': meta.get('webgui_required'),
            'webgui_field_required': meta.get('webgui_field_required'),
            'webgui_rich_code_cells_heuristic': rich,
            'webgui_applicability_review_required': meta.get('webgui_required') is not True,
        })
    tracked = {r['path'] for r in rows}
    present = {p.relative_to(ROOT).as_posix() for p in (ROOT / 'docs').rglob('*.ipynb')
               if '.ipynb_checkpoints' not in p.parts}
    report = {
        'schema': 'radia.docs.scholarship_audit.v1',
        'generated_at_utc': datetime.now(timezone.utc).isoformat(),
        'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'scope': 'Tracked public notebooks present in this worktree; metadata/prose audit, not solver acceptance.',
        'numerical_reexecution': False,
        'citations': citations,
        'saved_results': results,
        'notebooks': rows,
        'other_notebooks_present_not_in_tracked_audit': sorted(present - tracked),
        'webgui_applicability_review_count': sum(r['webgui_applicability_review_required'] for r in rows),
        'acceptance_limits': [
            'Citation resolution does not establish correct scientific attribution or a complete derivation.',
            'Existing saved outputs are not a rerun against current dirty solver sources.',
            'WebGUI metadata gates do not cover notebooks without explicit applicability declarations.',
            'AXIFEM catalog relocation still needs a genuine public result-bearing replacement.',
            'Wakao analytical-formula reports Part 1-9 still need individually verified bibliographic metadata.',
        ],
    }
    if live_bibliography:
        main = ROOT / 'packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib'
        a = {e.key: e.fields for e in read_bib_file(main)}
        b = {e.key: e.fields for e in read_bib_file(live_bibliography)}
        report['live_bibliography_comparison'] = {
            'path': str(live_bibliography), 'worktree_entries': len(a), 'live_entries': len(b),
            'shared_keys': len(a.keys() & b.keys()),
            'cited_keys_missing_in_live': sorted(cited - b.keys()),
            'shared_keys_with_field_differences': sorted(k for k in a.keys() & b.keys() if a[k] != b[k]),
            'status': 'equal' if a == b else 'divergent_requires_reconciliation',
            'note': 'This reads the live-reported source file; it does not prove a process reload.',
        }
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--live-bibliography', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = audit(args.live_bibliography)
    if args.output:
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({
        'citations': report['citations']['summary'],
        'saved_results': report['saved_results']['summary'],
        'code_and_outputs_unchanged': all(r['complete_code_cells_equal_head'] for r in report['notebooks']),
        'webgui_applicability_review_count': report['webgui_applicability_review_count'],
        'other_notebooks': len(report['other_notebooks_present_not_in_tracked_audit']),
        'live_bibliography': report.get('live_bibliography_comparison'),
    }, ensure_ascii=False))
