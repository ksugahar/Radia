"""Stage an explicitly selected research IH artifact with exact provenance."""
from pathlib import Path
import hashlib
import json
import shutil

FEATURES = {'element-temperature-bh', 'complete-flat-p1-reaction',
            'accepted-element-impedance', 'transactional-temperature-range'}
NATIVE_SOURCES = ('radia_cuda_ih.cu', 'ih_coupled_state.cuh', 'ih_complex_solver.cuh',
                  'ih_bicgstab_kernels.cuh', 'ih_thermal_kernels.cuh',
                  'ih_resident_kernels.cuh', 'ih_preconditioner.cuh',
                  'ih_esim_coupling.cuh', 'ih_temperature_material.cuh',
                  'ih_esim_cell.cuh', 'ih_peec_bem_assembly.cuh')
CORE_SOURCES = ('rad_bem_galerkin.cpp', 'rad_bem_galerkin.h',
                'rad_biot_savart_filaments.cpp', 'rad_biot_savart_filaments.h',
                'rad_host_parallel.h')
REQUIRED_SOURCES = ({'packages/radia-gpu/native/' + n for n in NATIVE_SOURCES}
                    | {'src/core/' + n for n in CORE_SOURCES}
                    | {'packages/radia-gpu/tools/build_ih_cuda.py',
                       'packages/radia-gpu/ih_payload.py', 'Build.ps1'})

MATLAB_FILES = ('radia_ih_cuda_sfun.m', '+radia/+simulink/ihResidentSFunction.m',
                '+radia/+simulink/makeCompletePEECBEMConfig.m',
                '+radia/+simulink/buildTemperatureBHResearchModel.m')
SAMPLES = ('makeTemperatureBHConfig.m', 'runTemperatureBHResearch.m')


def stage(root, package, commit, folder):
    root, package, folder = Path(root), Path(package), Path(folder)
    receipt = json.loads((folder / 'build.json').read_text(encoding='utf-8'))
    if (receipt.get('accepted') is not True or receipt.get('source_dirty') is not False
            or receipt.get('source_commit') != commit
            or receipt.get('taskmanager_regions') is not False
            or receipt.get('python_runtime_required') is not False
            or receipt.get('cpu_fallback') is not False
            or not FEATURES <= set(receipt.get('features', []))):
        raise RuntimeError('IH requires clean exact-source, native-only research provenance')
    sources = receipt.get('source_sha256', {})
    required = REQUIRED_SOURCES
    if not required <= set(sources):
        raise RuntimeError('Incomplete IH source receipt')
    for relative, digest in sources.items():
        path = Path(relative)
        if path.is_absolute() or '..' in path.parts:
            raise RuntimeError('Invalid IH source path')
        if hashlib.sha256((root / path).read_bytes()).hexdigest() != digest:
            raise RuntimeError('IH source receipt mismatch: ' + relative)
    binary = folder / 'radia_cuda_ih.mexw64'
    if hashlib.sha256(binary.read_bytes()).hexdigest() != receipt.get('binary_sha256'):
        raise RuntimeError('IH binary receipt mismatch')
    files = {'matlab/radia_cuda_ih.mexw64': binary}
    for name in MATLAB_FILES:
        files['matlab/' + name] = root / 'packages/radia-gpu/matlab' / name
    for name in SAMPLES:
        files['matlab/samples/ih-temperature-bh/' + name] = root / 'validation_test/ih_cuda' / name
    files['IH_TEMPERATURE_BH.md'] = root / 'packages/radia-gpu/docs/IH_TEMPERATURE_BH.md'
    files['IH_EXPERIMENT.md'] = root / 'validation_test/ih_cuda/TEMPERATURE_BH_EXPERIMENT.md'
    if any(not p.is_file() for p in files.values()):
        raise RuntimeError('Incomplete IH MATLAB/sample/documentation payload')
    for relative, source in files.items():
        target = package / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    public = {key: receipt[key] for key in
              ('source_commit', 'source_dirty', 'source_sha256', 'binary_sha256',
               'matlab', 'cuda_targets', 'taskmanager_regions',
               'python_runtime_required', 'cpu_fallback', 'features')}
    public.update(schema='radia.gpu.ih-research-payload.v1', research_only=True,
                  experimental_validation=False, independent_fem_validation=False)
    (package / 'ih-build.json').write_text(json.dumps(public, indent=2)+'\n', encoding='utf-8')
    return set(files) | {'ih-build.json'}
