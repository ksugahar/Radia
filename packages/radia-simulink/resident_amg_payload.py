"""Stage only explicitly qualified, exact-source resident AMG Python artifacts."""
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sysconfig


def stage(repo, package, commit, build):
    repo, package, build = Path(repo), Path(package), Path(build)
    name = '_resident_amg' + sysconfig.get_config_var('EXT_SUFFIX')
    receipt = json.loads((build / 'resident-amg-build.json').read_bytes())
    binary = build / name
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    if (receipt.get('schema') != 'radia.resident-amg-python-build.v1'
            or receipt.get('source_commit') != commit
            or receipt.get('source_dirty') is not False
            or receipt.get('binary') != name or receipt.get('binary_sha256') != digest
            or receipt.get('extension_suffix') != sysconfig.get_config_var('EXT_SUFFIX')
            or any(receipt.get(key) is not False for key in (
                'fault_injection', 'device_fault_injection', 'lifecycle_diagnostic',
                'matlab_engine_user_dependency'))):
        raise RuntimeError('Exact clean non-diagnostic resident AMG build required')
    required_sources = {
        'src/core/gpu/resident_amg_cuda.cu',
        'src/core/gpu/resident_amg_python.cpp',
        'src/core/gpu/resident_amg_host.hpp',
        'packages/radia-simulink/tools/build_resident_amg_python.py',
        'packages/radia-simulink/tools/build_cuda_mex.py',
        'packages/radia-simulink/tools/build_csr_python.py',
        'packages/radia-simulink/tools/native_source_paths.py',
        'packages/radia-simulink/tools/matlab_cuda_toolchain.py',
        'packages/radia-simulink/tools/vs2022_toolchain.py', 'Build.ps1',
    }
    sources = receipt.get('source_hashes', {})
    if set(sources) != required_sources or any(
            hashlib.sha256((repo / name).read_bytes()).hexdigest() != value
            for name, value in sources.items()):
        raise RuntimeError('Resident AMG source provenance mismatch')
    runtime_path = build / 'resident-amg-result.json'
    runtime = json.loads(runtime_path.read_bytes())
    if (runtime.get('schema') != 'radia.resident-amg-runtime.v1'
            or runtime.get('source_commit') != commit
            or runtime.get('binary_sha256') != digest
            or any(runtime.get(key) is not True for key in (
                'accepted', 'runtime_import_checked', 'redistribution_accepted',
                'matlab_dll_search_excluded'))
            or runtime.get('matlab_engine_loaded') is not False
            or runtime.get('matlab_dlls_loaded') is not False
            or runtime.get('runtime_provider') != 'nvidia-pip'
            or any(not isinstance(runtime.get(key), str) or not re.fullmatch(
                r'[0-9a-f]{64}', runtime[key]) for key in (
                    'loaded_runtime_sha256', 'redistribution_reference_sha256'))):
        raise RuntimeError('Resident AMG runtime and distribution qualification required')
    lifecycle = runtime.get('lifecycle', {})
    if any(lifecycle.get(key) is not True for key in (
            'repeated_solve', 'close_idempotent', 'solve_after_close_rejected',
            'invalid_input_recovery')):
        raise RuntimeError('Resident AMG lifecycle qualification required')
    cases = runtime.get('cases', [])
    if not isinstance(cases, list) or not cases:
        raise RuntimeError('Independent resident AMG numerical qualification required')
    for case in cases:
        if case.get('backend') != 'gpu' or case.get('cpu_oracle_independent') is not True:
            raise RuntimeError('Independent CPU oracle required for resident AMG')
        for key, limit in (('relative_l2', 1e-8), ('true_relative_residual', 1e-6)):
            value = case.get(key)
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not math.isfinite(value) or not 0 <= value <= limit):
                raise RuntimeError('Resident AMG numerical qualification failed')
    public = {key: receipt[key] for key in (
        'schema', 'source_commit', 'source_dirty', 'binary', 'binary_sha256',
        'extension_suffix', 'source_hashes', 'fault_injection',
        'device_fault_injection', 'lifecycle_diagnostic', 'matlab_engine_user_dependency')}
    public.update({key: runtime[key] for key in (
        'runtime_import_checked', 'redistribution_accepted',
        'matlab_dll_search_excluded', 'matlab_dlls_loaded',
        'loaded_runtime_sha256', 'redistribution_reference_sha256')})
    public['runtime_validation_sha256'] = hashlib.sha256(runtime_path.read_bytes()).hexdigest()
    package.mkdir(parents=True, exist_ok=True)
    shutil.copy2(binary, package / name)
    (package / 'resident-amg-build.json').write_text(json.dumps(public, indent=2) + '\n')
    return {name, 'resident-amg-build.json'}
