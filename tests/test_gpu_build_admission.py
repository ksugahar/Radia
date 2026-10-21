"""A new GPU worker requires exact host/source/task admission, not an env switch."""
import hashlib
import importlib.util
import sys
import json
from pathlib import Path
import pytest

path = Path(__file__).parents[1] / 'packages/radia-simulink/tools/matlab_cuda_toolchain.py'
sys.path.insert(0, str(path.parent))
spec = importlib.util.spec_from_file_location('admission_toolchain', path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_unadmitted_worker_and_hash_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(mod.socket, 'gethostname', lambda: 'new-worker')
    monkeypatch.setenv('RADIA_ALLOW_GPU_BUILD', '1')
    with pytest.raises(RuntimeError, match='admission'):
        mod.require_build_admission(tmp_path, 'hmatrix', None, None)
    file = tmp_path / 'admission.json'; file.write_text('{}')
    root = tmp_path / 'source'; root.mkdir()
    with pytest.raises(RuntimeError, match='hash'):
        mod.require_build_admission(root, 'hmatrix', file, '0'*64)


@pytest.mark.parametrize('change', [{}, {'lane':'projection'}, {'source_commit':'other'},
                                  {'execution_host':'other'}, {'review_accepted':False},
                                  {'task':'unknown'}, {'task':[]},
                                  {'limits':dict(cpu_count=4, job_commit_bytes=8*2**30, wall_seconds=1800)},
                                  {'limits':dict(cpu_count=4, job_commit_bytes=16*2**30, wall_seconds=901)}])
def test_exact_reviewed_binding(tmp_path, monkeypatch, change):
    monkeypatch.setattr(mod.socket, 'gethostname', lambda: 'new-worker')
    monkeypatch.setattr(mod.subprocess, 'check_output',
                        lambda args, **kwargs: '' if 'status' in args else 'abc123\n')
    record = dict(schema='radia.gpu-build-admission.v1', draft_only=False,
                  review_accepted=True, package='radia-simulink', task='hdiv-mmm', lane='hmatrix',
                  worker_role='gpu-build-worker', execution_host='new-worker', source_commit='abc123',
                  msvc_toolset='14.44', cuda_route='MATLAB-bundled-12.8',
                  limits=dict(cpu_count=4, job_commit_bytes=16*2**30, wall_seconds=900))
    record.update(change)
    file=tmp_path/'admission.json'; file.write_text(json.dumps(record))
    digest=hashlib.sha256(file.read_bytes()).hexdigest()
    root=tmp_path/'source'; root.mkdir()
    if change:
        with pytest.raises(RuntimeError, match='mismatch'):
            mod.require_build_admission(root,'hmatrix',file,digest)
    else:
        assert mod.require_build_admission(root,'hmatrix',file,digest)['admission_sha256']==digest

@pytest.mark.parametrize('change', [
    {}, {'task':'hdiv-mmm'}, {'task':'unknown'}, {'lane':'hmatrix'},
    {'task':'hdiv-mmm',
     'limits':dict(cpu_count=4, job_commit_bytes=16*2**30, wall_seconds=900)},
    {'limits':dict(cpu_count=5, job_commit_bytes=8*2**30, wall_seconds=1800)},
    {'limits':dict(cpu_count=4, job_commit_bytes=16*2**30, wall_seconds=1800)},
    {'limits':dict(cpu_count=4, job_commit_bytes=8*2**30, wall_seconds=1801)},
    {'limits':dict(cpu_count=4.0, job_commit_bytes=8*2**30, wall_seconds=1800)},
    {'limits':dict(cpu_count=True, job_commit_bytes=8*2**30, wall_seconds=1800)},
    {'limits':dict(cpu_count=4, job_commit_bytes=8*2**30)},
    {'limits':dict(cpu_count=4, job_commit_bytes=8*2**30, wall_seconds=1800, extra=1)},
])
def test_ih_declared_envelope_does_not_weaken_hdiv(tmp_path, monkeypatch, change):
    monkeypatch.setattr(mod.socket, 'gethostname', lambda: 'new-worker')
    monkeypatch.setattr(mod.subprocess, 'check_output',
                        lambda args, **kwargs: '' if 'status' in args else 'abc123\n')
    record = dict(schema='radia.gpu-build-admission.v1', draft_only=False,
                  execution_authorized=True, authorization_kind='explicit-user',
                  package='radia-simulink', task='ih', lane='ih',
                  worker_role='gpu-build-worker', execution_host='new-worker', source_commit='abc123',
                  msvc_toolset='14.44', cuda_route='MATLAB-bundled-12.8',
                  limits=dict(cpu_count=4, job_commit_bytes=8*2**30, wall_seconds=1800))
    record.update(change)
    file = tmp_path/'admission.json'; file.write_text(json.dumps(record))
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    root = tmp_path/'source'; root.mkdir()
    if change:
        with pytest.raises(RuntimeError, match='mismatch'):
            mod.require_build_admission(root, 'ih', file, digest)
    else:
        result = mod.require_build_admission(root, 'ih', file, digest)
        assert result['admission']['limits'] == record['limits']
        assert result['admission_sha256'] == digest


def test_explicit_ci_route_has_no_host_bypass(tmp_path, monkeypatch):
    monkeypatch.delenv('RADIA_GPU_CI_BUILD', raising=False)
    monkeypatch.setenv('GITHUB_ACTIONS', 'true')
    with pytest.raises(RuntimeError, match='admission'):
        mod.require_build_admission(tmp_path, 'hmatrix', None, None)
    monkeypatch.setenv('RADIA_GPU_CI_BUILD', '1')
    assert mod.require_build_admission(tmp_path, 'hmatrix', None, None)['route']=='explicit-ci-build'
    monkeypatch.setenv('GITHUB_ACTIONS', 'false')
    with pytest.raises(RuntimeError, match='admission'):
        mod.require_build_admission(tmp_path, 'hmatrix', None, None)


@pytest.mark.parametrize("script,lane", [
    ("build_ams_cuda.py", "ams"), ("build_cuda_mex.py", "csr-ic"),
    ("build_field_cuda.py", "field"), ("build_fmm_cuda.py", "fmm"),
    ("build_ih_cuda.py", "ih"),
])
def test_gateway_build_requires_unconditional_explicit_admission(script, lane):
    import ast
    tree = ast.parse((path.parent / script).read_text(encoding="utf-8"))
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Assign)
             and any(isinstance(target, ast.Name) and target.id == "admission"
                     for target in node.targets)]
    assert len(calls) == 1
    call = calls[0].value
    assert isinstance(call, ast.Call)
    assert isinstance(call.func, ast.Name) and call.func.id == "require_build_admission"
    assert isinstance(call.args[1], ast.Constant) and call.args[1].value == lane
    # Host identity may be recorded in a private receipt, never authorize a build.
    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            assert not any(isinstance(part, ast.Attribute) and part.attr == "gethostname"
                           for part in ast.walk(node.test))
