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
    wheel=dict(passed=True,version='0.1.0a1',source_commit='a'*40,wheel_sha256='b'*64,files=49)
    monkeypatch.setattr(gate,'verify',lambda path:wheel)
    (tmp_path/'radia_gpu-0.1.0a1-py3-none-win_amd64.whl').write_bytes(b'checked by wheel verifier')
    record=dict(schema='radia.gpu.installed-acceptance.v1',passed=True,ltspice_checked=True,
                source_commit='a'*40,wheel_sha256='b'*64,distribution_version='0.1.0a1',
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
