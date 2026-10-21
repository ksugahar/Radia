"""Stage explicitly selected, exact-commit bounded HDiv CUDA artifacts."""
from pathlib import Path
import hashlib
import json
import shutil

WRAPPERS = ('solveHDivElements.m', 'affineHDivQuadrature.m', 'curvedHDivQuadrature.m',
            'HMatrixOperator.m', 'HDivBHMaterial.m', 'PreparedHDivNonlinear.m',
            'makeHDivBHMaterial.m', 'solveNonlinearHDivElements.m')


SIMULINK_WRAPPERS = ('matlab/radia_hdiv_gpu_sfun.m',
    'matlab/+radia/+simulink/hdivGPUWorkflowSFunction.m',
    'matlab/+radia/+simulink/validateHDivGPUWorkflowConfig.m',
    'matlab/+radia/+simulink/HDivGPUWorkflow.md',
    'matlab/radia_hdiv_gpu_workflow.slx')

DOCUMENTATION_FILES = ('THIRD_PARTY_NOTICES.md', 'licenses/PYBIND11-LICENSE.txt',
                       'licenses/NGSOLVE-LICENSE.txt', 'licenses/NETGEN-LICENSE.txt',
                       'docs/python-gpu-runtime.md')


def stage_documentation(source, package):
    source, package = Path(source), Path(package)
    for name in DOCUMENTATION_FILES:
        target = package / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / name, target)
    return set(DOCUMENTATION_FILES)


def wheel_tag(python_selected, platform, version):
    if python_selected:
        if platform != 'win32':
            raise RuntimeError('Windows native binding wheel required')
        tag = f'cp{version[0]}{version[1]}'
        return tag, tag, 'win_amd64'
    return 'py3', 'none', 'win_amd64'


def runtime_contract(python_selected):
    interfaces = {
        'matlab_gpu': ['MATLAB R2026a', 'Parallel Computing Toolbox', 'NVIDIA GPU/driver'],
        'simulink_gpu': ['MATLAB R2026a', 'Parallel Computing Toolbox', 'Simulink', 'NVIDIA GPU/driver'],
        'python_cpu': ['Compatible Python, NumPy and NGSolve'],
        'ltspice': ['User-installed LTspice'],
    }
    if python_selected:
        interfaces['python_gpu'] = ['Qualified CPython/NGSolve ABI', 'NVIDIA GPU/driver',
                                    'Qualified NVIDIA CUDA runtime packages']
    return ['Windows x64'], interfaces


def stage_fmm(root, package, commit, folder):
    root, package, folder = Path(root), Path(package), Path(folder)
    receipt = json.loads((folder / 'build.json').read_text(encoding='utf-8'))
    if (not receipt.get('accepted') or receipt.get('source_dirty', True)
            or receipt.get('source_commit') != commit or receipt.get('taskmanager_used', True)
            or receipt.get('python_adapter_required', True)):
        raise RuntimeError('FMM requires clean exact-commit native provenance')
    for relative, digest in receipt['source_sha256'].items():
        if Path(relative).is_absolute() or '..' in Path(relative).parts:
            raise RuntimeError('Invalid FMM source path')
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError('FMM source mismatch: ' + relative)
    binary = folder / 'radia_cuda_fmm.mexw64'
    if hashlib.sha256(binary.read_bytes()).hexdigest() != receipt['binary_sha256']:
        raise RuntimeError('FMM binary mismatch')
    target = package / 'matlab' / binary.name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(binary, target)
    public = {key: receipt[key] for key in ('schema', 'source_commit', 'source_dirty',
              'source_sha256', 'binary_sha256', 'matlab', 'taskmanager_used', 'python_adapter_required')}
    (package / 'fmm-build.json').write_text(json.dumps(public, indent=2) + '\n')
    return {'matlab/' + binary.name, 'fmm-build.json'}


def stage(root, package, commit, hmatrix, projection):
    root, package = Path(root), Path(package)
    receipts = {}
    expected = set()
    for name, folder in [('hmatrix', hmatrix), ('projection', projection)]:
        folder = Path(folder)
        receipt = json.loads((folder / 'build.json').read_text(encoding='utf-8'))
        if receipt.get('fault_injection') is not False:
            raise RuntimeError('Diagnostic or unknown HDiv build cannot be staged for release: ' + name)
        if name == 'projection' and receipt.get('lifecycle_diagnostic') is not False:
            raise RuntimeError('Diagnostic or unknown projection lifecycle mode cannot be staged')
        if (not receipt.get('accepted') or receipt.get('source_dirty', True)
                or receipt.get('source_commit') != commit
                or receipt.get('taskmanager_started', True)):
            raise RuntimeError('HDiv requires clean exact-commit native provenance: ' + name)
        for relative, digest in receipt['source_sha256'].items():
            if Path(relative).is_absolute() or '..' in Path(relative).parts:
                raise RuntimeError('Invalid HDiv source path')
            source = root / relative
            if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != digest:
                raise RuntimeError('HDiv source mismatch: ' + relative)
        binary = folder / ('radia_cuda_' + name + '.mexw64')
        digest = hashlib.sha256(binary.read_bytes()).hexdigest()
        if digest != receipt['binary_sha256']:
            raise RuntimeError('HDiv binary mismatch: ' + name)
        target = package / 'matlab' / binary.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(binary, target)
        expected.add(target.relative_to(package).as_posix())
        receipts[name] = dict(source_commit=commit, source_dirty=False,
                              binary_sha256=digest, source_sha256=receipt['source_sha256'],
                              matlab=receipt['matlab'], taskmanager_started=False, fault_injection=False)
        if name == 'projection':receipts[name]['lifecycle_diagnostic']=False
    for name in WRAPPERS:
        relative = Path('matlab/+radia/+gpu') / name
        target = package / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / 'packages/radia-simulink' / relative, target)
        expected.add(relative.as_posix())
    for relative in SIMULINK_WRAPPERS:
        target = package / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / 'packages/radia-simulink' / relative, target)
        expected.add(relative)
    record = dict(schema='radia.gpu.hdiv-payload.v1', artifacts=receipts,
                  scope='Experimental bounded affine/P2 tetrahedra, HDiv orders 1/2, linear isotropic material',
                  limitations=['Dense near-rule Gram; general adaptive high-order HACApK unsupported.',
                               'Images, internal interfaces, constraints and nonlinear material unaccepted.',
                               'User-installed compatible NGSolve native libraries and Python DLL required; no Python solver adapter.'])
    (package / 'hdiv-build.json').write_text(json.dumps(record, indent=2) + '\n')
    expected.add('hdiv-build.json')
    return expected


def stage_python(root, package, commit, folder, runtime_path):
    """Stage the optional shared HDiv extension only with exact accepted proofs."""
    import math
    import re
    import sysconfig
    root,package,folder=Path(root),Path(package),Path(folder)
    build=json.loads((folder/'build.json').read_text(encoding='utf-8'))
    binding=build.get('python_binding',{})
    runtime=json.loads(Path(runtime_path).read_text(encoding='utf-8'))
    if (build.get('accepted') is not True or build.get('source_dirty') is not False
            or build.get('source_commit')!=commit or build.get('fault_injection') is not False):
        raise RuntimeError('HDiv Python requires a clean exact-source release build')
    suffix=sysconfig.get_config_var('EXT_SUFFIX');name='_hdiv_cuda'+suffix
    if (binding.get('fault_injection') is not False or binding.get('device_fault_injection') is not False or binding.get('schema')!='radia.hdiv-python-build.v1' or binding.get('binary')!=name
            or binding.get('extension_suffix')!=suffix or not re.fullmatch(r'_hdiv_cuda\.cp\d+-win_amd64\.pyd',name)
            or binding.get('matlab_engine_user_dependency') is not False
            or not binding.get('ngsolve_version') or not binding.get('pybind11_version')
            or not re.fullmatch(r'[0-9a-f]{64}',binding.get('pybind_header_sha256',''))):
        raise RuntimeError('HDiv Python binding ABI mismatch')
    for relative,digest in build.get('source_sha256',{}).items():
        path=Path(relative)
        if path.is_absolute() or '..' in path.parts or hashlib.sha256((root/path).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('HDiv Python source provenance mismatch')
    binary=folder/name;digest=hashlib.sha256(binary.read_bytes()).hexdigest()
    if digest!=binding.get('binary_sha256'):
        raise RuntimeError('HDiv Python binary provenance mismatch')
    if (runtime.get('schema')!='radia.hdiv-python-runtime.v1' or runtime.get('accepted') is not True
            or runtime.get('source_commit')!=commit or runtime.get('binary_sha256')!=digest
            or runtime.get('ngsolve_version')!=binding.get('ngsolve_version')
            or runtime.get('runtime_import_checked') is not True
            or runtime.get('matlab_engine_user_dependency') is not False
            or runtime.get('redistribution_accepted') is not True
            or not isinstance(runtime.get('redistribution_reference'),str) or not runtime['redistribution_reference']
            or not re.fullmatch(r'[0-9a-f]{64}',runtime.get('runtime_identity_sha256',''))):
        raise RuntimeError('HDiv Python exact runtime and redistribution evidence required')
    cases=runtime.get('cases',[])
    if not cases:raise RuntimeError('HDiv Python numerical evidence required')
    for case in cases:
        if case.get('backend')!='gpu':raise RuntimeError('HDiv Python explicit GPU evidence required')
        for key,limit in [('coefficient_error',1e-10),('independent_original_residual',1e-11)]:
            value=case.get(key)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=limit:
                raise RuntimeError('HDiv Python numerical evidence failed')
    public={key:binding[key] for key in ['schema','binary','binary_sha256','extension_suffix','ngsolve_version','pybind11_version','pybind_header_sha256']}
    public.update(fault_injection=False,device_fault_injection=False,source_commit=commit,matlab_engine_user_dependency=False,runtime_import_checked=True,redistribution_accepted=True,
                  redistribution_reference_sha256=hashlib.sha256(runtime['redistribution_reference'].encode()).hexdigest(),runtime_identity_sha256=runtime['runtime_identity_sha256'],
                  runtime_validation_sha256=hashlib.sha256(Path(runtime_path).read_bytes()).hexdigest())
    package.mkdir(parents=True,exist_ok=True);shutil.copy2(binary,package/name)
    (package/'hdiv-python-build.json').write_text(json.dumps(public,indent=2)+'\n',encoding='utf-8')
    return {name,'hdiv-python-build.json'}
