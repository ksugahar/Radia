"""Admission tests using non-executable synthetic payloads, never GPU evidence."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sysconfig
import sys
import zipfile

import pytest

ROOT = Path(__file__).parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def selected(tmp_path):
    repo, build, package = (tmp_path / name for name in ('repo', 'build', 'package'))
    build.mkdir()
    name = '_resident_amg' + sysconfig.get_config_var('EXT_SUFFIX')
    (build / name).write_bytes(b'not-a-native-library')
    sources = {}
    paths = ['src/core/gpu/resident_amg_' + suffix
             for suffix in ('cuda.cu', 'python.cpp', 'host.hpp')]
    paths += ['packages/radia-simulink/tools/' + name for name in (
        'build_resident_amg_python.py', 'build_cuda_mex.py',
        'build_csr_python.py', 'native_source_paths.py',
        'matlab_cuda_toolchain.py', 'vs2022_toolchain.py')]
    paths += ['Build.ps1']
    for relative in paths:
        path = repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(relative.encode())
        sources[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    receipt = dict(schema='radia.resident-amg-python-build.v1', source_commit='1' * 40,
                   source_dirty=False, binary=name,
                   extension_suffix=sysconfig.get_config_var('EXT_SUFFIX'),
                   binary_sha256=hashlib.sha256((build / name).read_bytes()).hexdigest(),
                   source_hashes=sources, fault_injection=False, device_fault_injection=False,
                   lifecycle_diagnostic=False, matlab_engine_user_dependency=False,
                   admission={'private': 'not-published'}, runtime_import_checked=False,
                   redistribution_accepted=False)
    runtime = dict(schema='radia.resident-amg-runtime.v1', source_commit='1' * 40,
                   binary_sha256=receipt['binary_sha256'], accepted=True,
                   runtime_import_checked=True, redistribution_accepted=True,
                   matlab_engine_loaded=False, matlab_dll_search_excluded=True,
                   matlab_dlls_loaded=False, runtime_provider='nvidia-pip',
                   loaded_runtime_sha256='2' * 64, redistribution_reference_sha256='3' * 64,
                   lifecycle={key: True for key in (
                       'repeated_solve', 'close_idempotent', 'solve_after_close_rejected',
                       'invalid_input_recovery')},
                   cases=[dict(backend='gpu', cpu_oracle_independent=True,
                               relative_l2=1e-12, true_relative_residual=1e-13)])
    return repo, build, package, receipt, runtime


@pytest.mark.parametrize('mutation', [
    'none', 'fault-missing', 'fault-integer', 'diagnostic', 'source', 'runtime-sha',
    'redistribution', 'provider', 'lifecycle', 'no-cases', 'cpu-oracle', 'residual',
    'dll-search-missing', 'dll-search-false', 'dll-search-integer',
    'dll-loaded-missing', 'dll-loaded-true', 'dll-loaded-integer',
])
def test_selected_payload_fails_closed_before_copy(selected, mutation):
    repo, build, package, receipt, runtime = selected
    if mutation.startswith('dll-'):
        key = 'matlab_dll_search_excluded' if mutation.startswith('dll-search-') else 'matlab_dlls_loaded'
        if mutation.endswith('missing'):
            runtime.pop(key)
        else:
            runtime[key] = {'false': False, 'true': True, 'integer': 0}[mutation.rsplit('-', 1)[1]]
    elif mutation == 'fault-missing':
        receipt.pop('fault_injection')
    elif mutation == 'fault-integer':
        receipt['fault_injection'] = 0
    elif mutation == 'diagnostic':
        receipt['lifecycle_diagnostic'] = True
    elif mutation == 'source':
        (repo / 'Build.ps1').write_text('changed')
    elif mutation == 'runtime-sha':
        runtime['binary_sha256'] = '0' * 64
    elif mutation == 'redistribution':
        runtime['redistribution_accepted'] = False
    elif mutation == 'provider':
        runtime['runtime_provider'] = 'matlab'
    elif mutation == 'lifecycle':
        runtime['lifecycle'].pop('close_idempotent')
    elif mutation == 'no-cases':
        runtime['cases'] = []
    elif mutation == 'cpu-oracle':
        runtime['cases'][0]['cpu_oracle_independent'] = False
    elif mutation == 'residual':
        runtime['cases'][0]['true_relative_residual'] = float('nan')
    (build / 'resident-amg-build.json').write_text(json.dumps(receipt))
    (build / 'resident-amg-result.json').write_text(json.dumps(runtime))
    stage = load('resident_amg_payload').stage
    if mutation != 'none':
        with pytest.raises(RuntimeError):
            stage(repo, package, '1' * 40, build)
        assert not package.exists()
        return
    files = stage(repo, package, '1' * 40, build)
    assert files == {receipt['binary'], 'resident-amg-build.json'}
    public = json.loads((package / 'resident-amg-build.json').read_bytes())
    assert 'admission' not in public
    assert public['runtime_import_checked'] is True
    tag = f'cp{sys.version_info.major}{sys.version_info.minor}'
    wheel = repo / f'radia_simulink-0.1.5-{tag}-{tag}-win_amd64.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        for name in files:
            archive.write(package / name, 'radia_simulink/' + name)
    manifest = dict(source_commit='1' * 40, files={receipt['binary']: receipt['binary_sha256']})
    with zipfile.ZipFile(wheel) as archive:
        load('verify_wheel').verify_python_binding(wheel, archive.namelist(), archive, manifest)
