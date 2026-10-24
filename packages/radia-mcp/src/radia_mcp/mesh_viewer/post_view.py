"""One-way VOL post views through Radia's production exporter, never a solver input."""
from __future__ import annotations

import asyncio
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from ._worker_process import run_worker, transaction_deadline


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _load_exporter():
    from radia import gmsh_post_export as exporter
    if (not hasattr(exporter, 'PayloadV42')
            or not callable(getattr(getattr(exporter, 'GmshPostExport', None), 'write_mesh', None))):
        raise RuntimeError('Installed radia lacks the production MSH 4.2 exporter; update radia')
    return exporter


def _require_exporter_capability(*, deadline=None) -> None:
    asyncio.run(run_worker(
        [sys.executable, '-B', '-m', 'radia_mcp.mesh_viewer.post_view', '--check-exporter'],
        env=os.environ, deadline=deadline if deadline is not None else transaction_deadline(),
        phase='VOL post-view exporter preflight'))

def _publish_mesh(staged: Path, target: Path) -> str:
    """Publish only a completed verified copy; linking never overwrites a target."""
    fd, name = tempfile.mkstemp(prefix='.' + target.name + '-', suffix='.partial', dir=target.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, 'wb') as destination, staged.open('rb') as original:
            for block in iter(lambda: original.read(1024 * 1024), b''):
                destination.write(block)
        digest = sha256(temporary)
        if digest != sha256(staged):
            raise RuntimeError('Post-view byte-copy hash mismatch')
        # Requires hard-link support; an unsupported output filesystem fails
        # explicitly. There is no overwrite or partial-copy fallback.
        os.link(temporary, target)
        return digest
    finally:
        temporary.unlink(missing_ok=True)


def validate_post_request(source: Path, output: str, checked_report: str,
                          source_sha256: str, group_scope: str, *, probe_exporter=True, deadline=None) -> tuple[Path, dict]:
    """Require a prior check-vol receipt plus the caller's checked-input digest."""
    if source.suffix.lower() != '.vol' or group_scope not in {'materials', 'boundaries'}:
        raise ValueError('VOL post-view requires group_scope materials or boundaries')
    target = Path(output).expanduser().resolve()
    if target.suffix != '.msh' or not target.parent.is_dir():
        raise ValueError('Post output must be a new .msh in an existing directory')
    if target.parent == source.resolve().parent or target.exists():
        raise ValueError('Post output must not exist or be beside/over the source VOL')
    if any(p.exists() for p in (Path(str(target) + '.opt'), target.with_suffix('.opt'))):
        raise ValueError('Post output must have no adjacent option script')
    if len(source_sha256) != 64 or sha256(source) != source_sha256.lower():
        raise ValueError('Checked source SHA256 mismatch')
    handoff = json.loads(Path(checked_report).read_text(encoding='utf-8-sig'))
    if (handoff.get('schema') != 'radia.checked-vol-handoff.v1'
            or handoff.get('source_sha256') != source_sha256.lower()):
        raise ValueError('Checked handoff must bind the check-vol report to this source SHA256')
    report = handoff.get('report', {})
    if (report.get('schema') != 'cubit-mesh-export.vol-check.v1'
            or report.get('passed') is not True
            or report.get('labels', {}).get('passed') is not True
            or report.get('labels', {}).get('schema') != 'radia.vol-label-contract.v1'):
        raise ValueError('A passing check-vol report with versioned label contract is required')
    if report.get('mesh', {}).get('curve_order') != 1:
        raise ValueError('Curved VOL post-view needs a retained-geometry protocol; unsupported here')
    scratch = Path(os.environ.get('RADIA_MESH_VIEWER_SCRATCH', str(target.parent))).resolve()
    if not scratch.is_dir() or not str(scratch).isascii() or scratch == source.resolve().parent:
        raise ValueError('Set RADIA_MESH_VIEWER_SCRATCH to an existing ASCII job directory outside the source directory')
    if probe_exporter:
        _require_exporter_capability(deadline=deadline)
    return target, report


def convert_post_view(source: Path, output: str, checked_report: str,
                      source_sha256: str, group_scope: str) -> dict:
    """Keep exporter companions in staging; copy only identical checked MSH bytes."""
    target, report = validate_post_request(source, output, checked_report,
                                           source_sha256, group_scope, probe_exporter=False)
    exporter = _load_exporter()
    scratch = Path(os.environ.get('RADIA_MESH_VIEWER_SCRATCH', str(target.parent))).resolve()
    staging = Path(tempfile.mkdtemp(prefix='.radia-post-', dir=scratch))
    staged_source = staging / 'source.vol'
    with staged_source.open('xb') as copy, source.open('rb') as original:
        for block in iter(lambda: original.read(1024 * 1024), b''):
            copy.write(block)
    if sha256(staged_source) != source_sha256.lower():
        raise ValueError('Staged VOL byte-copy hash mismatch')
    import ngsolve as ng
    mesh = ng.Mesh(staged_source.as_posix())
    if mesh.dim not in (2, 3) or mesh.ngmesh.GetCurveOrder() != 1:
        raise ValueError('Only checked linear 2D/3D VOL maps are supported')
    expected_mesh = report['mesh']
    actual = (mesh.dim, mesh.ne, mesh.nv)
    expected = (expected_mesh['dimension'], expected_mesh['n_elements'], expected_mesh['n_points'])
    if actual != expected:
        raise ValueError(f'Checked mesh inventory mismatch: actual={actual}, checked={expected}')
    for key, names in (('materials', mesh.GetMaterials()), ('boundaries', mesh.GetBoundaries())):
        if set(names) != set(report['labels'][key]):
            raise ValueError(f'Checked {key} differ from loaded VOL')
    boundary = group_scope == 'boundaries'
    if mesh.dim == 2 and boundary:
        raise ValueError('2D boundary-line export is unsupported')
    dimension = mesh.dim - int(boundary)
    names = mesh.GetBoundaries() if boundary else mesh.GetMaterials()
    expected_groups = {(dimension, name) for name in names}
    if not expected_groups or any(not name for _, name in expected_groups):
        raise ValueError('Selected region scope has no named physical groups')
    staged = staging / 'view.msh'
    exporter.GmshPostExport(mesh, boundary=boundary).write_mesh(str(staged))
    from radia_mcp.gmsh.runtime import require_msh42
    from radia_mcp.gmsh.msh_inspect import inspect_msh
    require_msh42(staged)
    inspection = inspect_msh(staged)
    observed_groups = {(p['dim'], p['name']) for p in inspection['physical_names'] or []}
    if observed_groups != expected_groups:
        raise ValueError(f'Exporter physical groups mismatch: {observed_groups} != {expected_groups}')
    selected_count = sum(1 for _ in mesh.Elements(ng.BND if boundary else ng.VOL))
    actual_count = inspection.get('elements', {}).get('count')
    if (inspection.get('parse_errors') or not inspection.get('elements', {}).get('parse_complete')
            or actual_count != selected_count):
        raise ValueError(f'Exporter selected element count mismatch: {actual_count} != {selected_count}')
    if sha256(source) != source_sha256.lower():
        raise ValueError('Source VOL changed during conversion')
    # Filesystem copy is not a second MSH writer.
    digest = _publish_mesh(staged, target)
    return {'schema': 'radia.mesh-post-view.v1', 'output': str(target), 'sha256': digest, 'staging': str(staging),
            'source_sha256': source_sha256.lower(), 'group_scope': group_scope,
            'physical_groups': sorted(observed_groups), 'solver_input': False}


async def prepare_post_view_async(source: Path, output: str, checked_report: str,
                                  source_sha256: str, group_scope: str, *,
                                  python=None, env=None, deadline=None) -> dict:
    """Await verified conversion without importing native libraries into the caller."""
    target = Path(output).expanduser().resolve()
    if source.suffix.lower() != '.vol' or target.suffix != '.msh' or not target.parent.is_dir():
        raise ValueError('VOL post-view requires an existing source and new .msh output directory')
    if target.exists() or target.parent == source.resolve().parent:
        raise ValueError('Post output must not exist or be beside/over the source VOL')
    scratch = Path(os.environ.get('RADIA_MESH_VIEWER_SCRATCH', str(target.parent))).resolve()
    if not scratch.is_dir() or not str(scratch).isascii() or scratch == source.resolve().parent:
        raise ValueError('Set RADIA_MESH_VIEWER_SCRATCH to an existing ASCII job directory outside the source directory')
    receipt_root = Path(tempfile.mkdtemp(prefix='.radia-post-receipt-', dir=scratch))
    receipt = receipt_root / 'receipt.json'
    argv = [str(python or sys.executable), '-B', '-m', 'radia_mcp.mesh_viewer.post_view',
            str(source), str(target), checked_report, source_sha256, group_scope,
            '--receipt', str(receipt)]
    await run_worker(argv, env=env if env is not None else os.environ,
                     deadline=deadline if deadline is not None else transaction_deadline(),
                     phase=f'VOL conversion {source}')
    if receipt.stat().st_size > 1024 * 1024:
        raise RuntimeError('Conversion receipt exceeds 1 MiB; verified artifact retained at ' + str(target))
    data = json.loads(receipt.read_text(encoding='utf-8'))
    import re
    if (not target.is_file() or data.get('schema') != 'radia.mesh-post-view.v1' or data.get('output') != str(target)
            or data.get('source_sha256') != source_sha256.lower()
            or data.get('group_scope') != group_scope or data.get('solver_input') is not False
            or not isinstance(data.get('staging'), str) or not data.get('physical_groups')
            or not re.fullmatch('[0-9a-f]{64}', data.get('sha256', ''))):
        raise RuntimeError('Unbound conversion receipt; no GUI launched: ' + str(receipt))
    return {**data, 'worker_receipt': str(receipt)}


def prepare_post_view(source: Path, output: str, checked_report: str,
                      source_sha256: str, group_scope: str, *, deadline=None) -> dict:
    return asyncio.run(prepare_post_view_async(source, output, checked_report,
                                               source_sha256, group_scope, deadline=deadline))

def main() -> int:
    import argparse
    from contextlib import redirect_stdout
    if sys.argv[1:] == ['--check-exporter']:
        try:
            _load_exporter()
            return 0
        except Exception as error:
            print(str(error), file=sys.stderr)
            return 1
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output')
    parser.add_argument('checked_report')
    parser.add_argument('source_sha256')
    parser.add_argument('group_scope', choices=['materials', 'boundaries'])
    parser.add_argument('--receipt', type=Path, required=True)
    args = parser.parse_args()
    try:
        with redirect_stdout(sys.stderr):
            result = convert_post_view(args.source, args.output, args.checked_report,
                                       args.source_sha256, args.group_scope)
        payload = json.dumps(result, indent=2)
        with args.receipt.open('x', encoding='utf-8') as stream:
            stream.write(payload)
        return 0
    except Exception as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
