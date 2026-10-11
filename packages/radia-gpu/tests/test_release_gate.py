"""Publication must use a successful main build and exact installed acceptance."""
import importlib.util
import json
from pathlib import Path
import pytest


def test_release_rejects_foreign_ci_and_mismatched_installed_artifact(tmp_path,monkeypatch):
    root=Path(__file__).resolve().parents[1]
    monkeypatch.syspath_prepend(str(root))
    spec=importlib.util.spec_from_file_location('gpu_release_gate',root/'release_gate.py')
    gate=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    run=dict(id=42,status='completed',conclusion='success',event='push',head_branch='main',
             path='.github/workflows/radia-gpu.yml',head_sha='a'*40,
             head_repository={'full_name':'owner/repository'})
    wheel=dict(passed=True,version='0.1.0a1',source_commit='a'*40,wheel_sha256='b'*64,files=49,
               manifest_sha256='d'*64)
    monkeypatch.setattr(gate,'verify',lambda path:wheel)
    (tmp_path/'radia_gpu-0.1.0a1-py3-none-win_amd64.whl').write_bytes(b'checked by wheel verifier')
    record=dict(schema='radia.gpu.installed-acceptance.v1',passed=True,ltspice_checked=True,
                source_commit='a'*40,wheel_sha256='b'*64,distribution_version='0.1.0a1',
                manifest_sha256='d'*64,
                simulink_repeated_runs=2,simulink_failure_cleanup=True)
    evidence=tmp_path/'result.json'
    evidence.write_text(json.dumps(record))
    qualify=lambda candidate=run:gate.qualify(candidate,'owner/repository',42,tmp_path,'b'*64)
    assert qualify()['tag']=='radia-gpu-v0.1.0a1'
    for overrides in ({'conclusion':'failure'},{'head_branch':'feature'},{'event':'pull_request'},
                      {'head_repository':{'full_name':'foreign/repository'}},
                      {'path':'.github/workflows/radia-optuna.yml'}):
        with pytest.raises(ValueError,match='successful repository main'):
            qualify(dict(run,**overrides))
    record['wheel_sha256']='c'*64
    evidence.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='exact release wheel'):
        qualify()
    wheel['version']='0.1.0'
    with pytest.raises(ValueError,match='independent review'):
        qualify()


def test_stable_release_requires_review_and_same_wheel_build_host_parity(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[1]
    monkeypatch.syspath_prepend(str(root))
    spec = importlib.util.spec_from_file_location('gpu_stable_gate', root/'release_gate.py')
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    wheel = dict(passed=True, version='0.1.0', source_commit='a'*40,
                 wheel_sha256='b'*64, files=115, cuda_solvers=['csr_pcg', 'ams'],
                 manifest_sha256='d'*64)
    monkeypatch.setattr(gate, 'verify', lambda path: wheel)
    artifact = tmp_path/'artifact'
    artifact.mkdir()
    (artifact/'radia_gpu-0.1.0-py3-none-win_amd64.whl').write_bytes(b'checked separately')
    installed = dict(schema='radia.gpu.installed-acceptance.v1', passed=True,
                     source_commit=wheel['source_commit'], wheel_sha256=wheel['wheel_sha256'],
                     distribution_version=wheel['version'], ltspice_checked=True,
                     manifest_sha256=wheel['manifest_sha256'],
                     simulink_repeated_runs=2, simulink_failure_cleanup=True)
    (artifact/'result.json').write_text(json.dumps(installed))
    run = dict(id=42, status='completed', conclusion='success', event='push',
               head_branch='main', path='.github/workflows/radia-gpu.yml',
               head_sha=wheel['source_commit'], head_repository={'full_name':'owner/repository'})
    qualification = dict(schema='radia.gpu.stable-qualification.v1', version=wheel['version'],
                         source_commit=wheel['source_commit'], wheel_sha256=wheel['wheel_sha256'],
                         independent_review=dict(source_commit=wheel['source_commit'],
                             outcome='accepted', unresolved_findings=[], reviewer='independent Codex',
                             report='Completed read-only review of the exact release source.'),
                         build_host_acceptance=dict(installed, execution_role='cuda-build-host', device='Quadro RTX 5000',
                             cpu_parity_checked=True, cpu_parity={
                                 'csr_pcg':dict(relative_l2=2e-10, true_relative_residual=1e-10),
                                 'ams':dict(relative_l2=3e-9, true_relative_residual=2e-10)}))
    qualify = lambda: gate.qualify(run, 'owner/repository', 42, artifact, 'b'*64, tmp_path)
    path = tmp_path/'0.1.0.json'
    with pytest.raises(ValueError, match='independent review'):
        qualify()
    path.write_text(json.dumps(qualification))
    assert qualify()['prerelease'] is False
    mutations = [
        ('source_commit', 'c'*40),
        ('independent_review.source_commit', 'c'*40),
        ('independent_review.outcome', 'pending'),
        ('independent_review.unresolved_findings', ['P2 still open']),
        ('independent_review.report', ''),
        ('build_host_acceptance.wheel_sha256', 'c'*64),
        ('build_host_acceptance.manifest_sha256', 'c'*64),
        ('build_host_acceptance.device', 'NVIDIA A100'),
        ('build_host_acceptance.execution_role', 'cpu-distribution-host'),
        ('build_host_acceptance.ltspice_checked', False),
        ('build_host_acceptance.cpu_parity_checked', False),
        ('build_host_acceptance.simulink_failure_cleanup', False),
        ('build_host_acceptance.simulink_repeated_runs', float('nan')),
        ('build_host_acceptance.simulink_repeated_runs', float('inf')),
        ('build_host_acceptance.cpu_parity', {}),
        ('build_host_acceptance.cpu_parity', {'csr_pcg':dict(relative_l2=0, true_relative_residual=0)}),
        ('build_host_acceptance.cpu_parity.csr_pcg.relative_l2', 1e-5),
        ('build_host_acceptance.cpu_parity.csr_pcg.true_relative_residual', float('nan')),
    ]
    for key, value in mutations:
        candidate = json.loads(json.dumps(qualification))
        current = candidate
        parts = key.split('.')
        for part in parts[:-1]:
            current = current[part]
        current[parts[-1]] = value
        path.write_text(json.dumps(candidate))
        with pytest.raises(ValueError):
            qualify()
    path.write_text(json.dumps(qualification))
    for value in (float('nan'), float('inf'), True, 2.5):
        invalid = dict(installed, simulink_repeated_runs=value)
        (artifact/'result.json').write_text(json.dumps(invalid))
        with pytest.raises(ValueError, match='exact release wheel'):
            qualify()
