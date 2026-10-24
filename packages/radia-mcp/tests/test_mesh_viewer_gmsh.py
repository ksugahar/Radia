"""Guards: strict post-view formats, no option-script execution, no source writes.
Assumes: native GUI and production-exporter proof is collected separately.
"""
import hashlib
import asyncio
import json
from pathlib import Path
import sys
from types import SimpleNamespace as NS

import pytest
from radia_mcp.mesh_viewer import launcher, post_view

EXPORTER_PROBE = post_view._require_exporter_capability


def mesh_text(version='4.2', name='iron'):
    return f'$MeshFormat\n{version} 0 8\n$EndMeshFormat\n$PhysicalNames\n1\n3 1 "{name}"\n$EndPhysicalNames\n'


def test_msh_version_and_adjacent_script_preflight(tmp_path):
    path = tmp_path / '日本語 space ;$.msh'
    path.write_text(mesh_text('4.1'), encoding='utf-8')
    with pytest.raises(ValueError, match='4.2'): launcher.validate_mesh_path(str(path))
    path.write_text(mesh_text(), encoding='utf-8')
    assert launcher.validate_mesh_path(str(path)) == path
    Path(str(path) + '.opt').write_text('untrusted', encoding='utf-8')
    with pytest.raises(ValueError, match='option'): launcher.validate_mesh_path(str(path))


def test_supplied_extension_cannot_be_changed_by_resolution(tmp_path, monkeypatch):
    target = tmp_path / 'actual.msh'; target.write_text(mesh_text())
    with pytest.raises(ValueError, match='existing .vol or .msh'):
        launcher.validate_mesh_path(str(tmp_path / 'alias.txt'))
    monkeypatch.setattr(Path, 'resolve', lambda self, **kwargs: target)
    with pytest.raises(ValueError, match='must match'):
        launcher.validate_mesh_path(str(tmp_path / 'alias.vol'))


def test_gmsh_uses_no_configuration_or_qualification_writer(tmp_path, monkeypatch):
    calls = []
    option = NS(setNumber=lambda *a: calls.append(a),
                getNumber=lambda _: 4.2, getString=lambda _: '5.0.0-git-test')
    fake = NS(__version__='5.0.0', option=option,
              initialize=lambda *a, **k: calls.append(('initialize', a, k)),
              open=lambda p: calls.append(('open', p)),
              fltk=NS(run=lambda: calls.append(('run',))),
              finalize=lambda: calls.append(('finalize',)))
    monkeypatch.setitem(sys.modules, 'gmsh', fake)
    path = tmp_path / '日本語 [view] ;$.msh'
    launcher.open_mesh(path)
    assert calls[0][2] == {'readConfigFiles': False}
    assert ('General.SaveOptions', 0) in calls and ('General.SaveSession', 0) in calls
    assert ('open', str(path)) in calls and calls[-1] == ('finalize',)


def test_missing_gmsh_names_declared_extra(monkeypatch):
    from radia_mcp.gmsh import runtime
    def absent(_): raise runtime.importlib.metadata.PackageNotFoundError('gmsh')
    monkeypatch.setattr(runtime.importlib.metadata, 'version', absent)
    with pytest.raises(RuntimeError, match=r'radia-mcp\[gmsh\]') as error:
        runtime.verify_runtime_identity()
    assert 'pip install --pre --upgrade -i https://gmsh.info/python-packages-dev gmsh' in str(error.value)
    assert 'then pip install "radia-mcp[gmsh]"' in str(error.value)


@pytest.fixture
def post_case(tmp_path, monkeypatch):
    monkeypatch.setattr(post_view, '_require_exporter_capability', lambda: None)
    srcdir = tmp_path / 'inputs'; srcdir.mkdir()
    outdir = tmp_path / 'outputs'; outdir.mkdir()
    source = srcdir / '日本語 space.vol'; source.write_bytes(b'checked fixture')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    report = {'schema': 'cubit-mesh-export.vol-check.v1', 'passed': True,
              'mesh': {'dimension': 3, 'curve_order': 1, 'n_elements': 6, 'n_points': 8},
              'labels': {'schema': 'radia.vol-label-contract.v1', 'passed': True,
                         'materials': ['iron'], 'boundaries': ['wall']}}
    handoff = {'schema': 'radia.checked-vol-handoff.v1', 'source_sha256': digest, 'report': report}
    receipt = tmp_path / 'check.json'; receipt.write_text(json.dumps(handoff))
    mesh = NS(dim=3, ne=6, nv=8, ngmesh=NS(GetCurveOrder=lambda: 1),
              GetMaterials=lambda: ('iron',), GetBoundaries=lambda: ('wall',),
              Elements=lambda region: range(6 if region == 'VOL' else 12))
    monkeypatch.setitem(sys.modules, 'ngsolve', NS(Mesh=lambda _: mesh, VOL='VOL', BND='BND'))
    from radia_mcp.gmsh import msh_inspect
    monkeypatch.setattr(msh_inspect, 'inspect_msh', lambda _: {
        'physical_names': [{'dim': 3, 'name': 'iron'}],
        'elements': {'count': 6, 'parse_complete': True}, 'parse_errors': []})
    writes = []
    class Exporter:
        def __init__(self, mesh, boundary): self.boundary = boundary
        def write_mesh(self, filename):
            writes.append(filename)
            Path(filename).write_text(mesh_text())
            Path(filename + '.opt').write_text('never execute')
    monkeypatch.setitem(sys.modules, 'radia', NS(gmsh_post_export=NS(PayloadV42=object, GmshPostExport=Exporter)))
    return source, outdir / 'view.msh', receipt, digest, writes


def test_post_preflight_rejects_missing_exporter_before_writes(post_case, monkeypatch):
    source, output, receipt, digest, writes = post_case
    monkeypatch.setattr(post_view, '_require_exporter_capability', EXPORTER_PROBE)
    def probe(argv, **kwargs):
        assert argv[-1] == '--check-exporter'
        assert kwargs['timeout'] == 20 and 'shell' not in kwargs
        return NS(returncode=1, stderr='Installed radia lacks the production MSH 4.2 exporter; update radia')
    monkeypatch.setattr(post_view.subprocess, 'run', probe)
    with pytest.raises(RuntimeError, match='Installed radia lacks'):
        post_view.validate_post_request(source, str(output), str(receipt), digest, 'materials')
    assert not output.exists() and writes == []


def test_exporter_check_requires_the_production_api(monkeypatch):
    monkeypatch.setitem(sys.modules, 'radia', NS(gmsh_post_export=NS()))
    with pytest.raises(RuntimeError, match='MSH 4.2 exporter'):
        post_view._load_exporter()


def test_mcp_returns_preflight_error_instead_of_launch_acceptance(post_case, monkeypatch):
    from radia_mcp.mesh_viewer import mcp as adapter
    source, output, receipt, digest, _ = post_case
    monkeypatch.setattr(adapter, 'require_interactive_desktop', lambda: {})
    monkeypatch.setattr(adapter.subprocess, 'CREATE_NO_WINDOW', 0, raising=False)
    async def communicate():
        return b'', b'Installed radia lacks the production MSH 4.2 exporter; update radia'
    async def spawn_probe(*args, **kwargs):
        assert '--check' in args and str(output) in args
        return NS(returncode=1, communicate=communicate)
    monkeypatch.setattr(adapter.asyncio, 'create_subprocess_exec', spawn_probe)
    monkeypatch.setattr(adapter.subprocess, 'Popen', lambda *a, **k: pytest.fail('must not launch GUI'))
    with pytest.raises(RuntimeError, match='Viewer preflight failed: Installed radia lacks'):
        asyncio.run(adapter.open_mesh_viewer(str(source), str(output), str(receipt), digest, 'materials'))


def test_post_copy_is_identical_and_keeps_companions_in_staging(post_case):
    source, output, receipt, digest, writes = post_case
    result = post_view.convert_post_view(source, str(output), str(receipt), digest, 'materials')
    assert len(writes) == 1
    assert output.read_bytes() == Path(writes[0]).read_bytes()
    assert output.parent != source.parent and not Path(str(output) + '.opt').exists()
    assert post_view.sha256(source) == digest and result['solver_input'] is False
    with pytest.raises(ValueError, match='exist'):
        post_view.convert_post_view(source, str(output), str(receipt), digest, 'materials')


def test_post_fails_on_unpreserved_groups_and_never_creates_final_output(post_case):
    source, output, receipt, digest, _ = post_case
    with pytest.raises(ValueError, match='physical groups'):
        post_view.convert_post_view(source, str(output), str(receipt), digest, 'boundaries')
    assert not output.exists()


def test_post_rejects_source_directory_and_digest_mismatch(post_case):
    source, output, receipt, digest, _ = post_case
    with pytest.raises(ValueError, match='beside'):
        post_view.validate_post_request(source, str(source.with_suffix('.msh')), str(receipt), digest, 'materials')
    with pytest.raises(ValueError, match='SHA256'):
        post_view.validate_post_request(source, str(output), str(receipt), '0' * 64, 'materials')


def test_post_rejects_a_handoff_for_different_checked_bytes(post_case):
    source, output, receipt, digest, _ = post_case
    handoff = json.loads(receipt.read_text())
    handoff['source_sha256'] = '0' * 64
    receipt.write_text(json.dumps(handoff))
    with pytest.raises(ValueError, match='bind'):
        post_view.convert_post_view(source, str(output), str(receipt), digest, 'materials')
    assert not output.exists()


def test_post_rejects_partial_element_export(post_case, monkeypatch):
    source, output, receipt, digest, _ = post_case
    from radia_mcp.gmsh import msh_inspect
    monkeypatch.setattr(msh_inspect, 'inspect_msh', lambda _: {
        'physical_names': [{'dim': 3, 'name': 'iron'}],
        'elements': {'count': 5, 'parse_complete': True}, 'parse_errors': []})
    with pytest.raises(ValueError, match='element count'):
        post_view.convert_post_view(source, str(output), str(receipt), digest, 'materials')
    assert not output.exists()


def test_unicode_output_uses_an_explicit_ascii_job_scratch(post_case, monkeypatch):
    source, output, receipt, digest, writes = post_case
    unicode_dir = output.parent / '日本語 出力'; unicode_dir.mkdir()
    output = unicode_dir / 'view.msh'
    with pytest.raises(ValueError, match='ASCII job'):
        post_view.validate_post_request(source, str(output), str(receipt), digest, 'materials')
    monkeypatch.setenv('RADIA_MESH_VIEWER_SCRATCH', str(output.parent.parent))
    result = post_view.convert_post_view(source, str(output), str(receipt), digest, 'materials')
    assert Path(result['staging']).parent == output.parent.parent
    assert Path(writes[0]).read_bytes() == output.read_bytes()


def test_conversion_process_is_separate_and_uses_a_file_receipt(post_case, monkeypatch):
    source, output, receipt, digest, _ = post_case
    calls = []
    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        Path(argv[-1]).write_text(json.dumps({'output': str(output), 'solver_input': False}))
        return NS(returncode=0, stdout='native log noise', stderr='')
    monkeypatch.setattr(post_view.subprocess, 'run', run)
    monkeypatch.setitem(sys.modules, 'ngsolve', None)
    result = post_view.prepare_post_view(source, str(output), str(receipt), digest, 'materials')
    assert calls[0][0][1:4] == ['-B', '-m', 'radia_mcp.mesh_viewer.post_view']
    assert 'shell' not in calls[0][1] and calls[0][1]['stdin'] == post_view.subprocess.DEVNULL
    assert result['output'] == str(output) and Path(result['worker_receipt']).is_file()


def test_failed_worker_does_not_leave_an_empty_json_receipt(post_case, monkeypatch):
    source, output, receipt, digest, _ = post_case
    receipts = []
    def fail(argv, **kwargs):
        receipts.append(Path(argv[-1]))
        assert not receipts[-1].exists()
        return NS(returncode=1, stderr='export failed', stdout='')
    monkeypatch.setattr(post_view.subprocess, 'run', fail)
    with pytest.raises(RuntimeError, match='export failed'):
        post_view.prepare_post_view(source, str(output), str(receipt), digest, 'materials')
    assert len(receipts) == 1 and not receipts[0].exists()
    assert receipts[0].parent.is_dir() and not output.exists()
