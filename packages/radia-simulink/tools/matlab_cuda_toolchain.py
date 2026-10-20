"""Select VS2022 and MATLAB's own CUDA; never search external toolkits."""
from pathlib import Path
import hashlib
import os
import re
import json
import socket
import subprocess
from vs2022_toolchain import discover_vs, compiler_bin


def require_build_admission(root, lane, path, expected_sha256):
    """Bind a GPU build to a reviewed private host/task/source admission."""
    if path is None and expected_sha256 is None:
        if os.environ.get('RADIA_GPU_CI_BUILD') != '1' or os.environ.get('GITHUB_ACTIONS') != 'true':
            raise RuntimeError('Explicit GPU build admission or RADIA_GPU_CI_BUILD=1 in GitHub Actions required')
        return {'route': 'explicit-ci-build', 'resource_enforcement': 'CI runner ownership'}
    if path is None or not expected_sha256:
        raise RuntimeError('Explicit GPU build admission and hash required; no host/environment fallback')
    root = Path(root).resolve()
    path = Path(path).resolve(strict=True)
    if path.is_relative_to(root):
        raise RuntimeError('Private build admission must be outside the source checkout')
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected_sha256:
        raise RuntimeError('GPU build admission hash mismatch')
    record = json.loads(path.read_text(encoding='utf-8'))
    commit = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', '-C', str(root), 'status', '--porcelain'], text=True).strip():
        raise RuntimeError('Admitted GPU source must be clean')
    if (record.get('schema') != 'radia.gpu-build-admission.v1'
            or record.get('draft_only') is not False
            or not (record.get('review_accepted') is True or
                    (record.get('execution_authorized') is True and record.get('authorization_kind') == 'explicit-user'))
            or record.get('package') != 'radia-simulink' or record.get('task') != 'hdiv-mmm'
            or record.get('lane') != lane or record.get('worker_role') != 'gpu-build-worker'
            or record.get('execution_host', '').casefold() != socket.gethostname().casefold()
            or record.get('source_commit') != commit
            or record.get('msvc_toolset') != '14.44' or record.get('cuda_route') != 'MATLAB-bundled-12.8'
            or record.get('limits') != {'cpu_count': 4, 'job_commit_bytes': 16 * 2**30, 'wall_seconds': 900}):
        raise RuntimeError('GPU build host/task/source/toolchain admission mismatch')
    return dict(admission_sha256=digest, admission=record,
                resource_enforcement='Required external owned Job Object; this record is not an actuator')


def reject_overrides(environment=None):
    env = os.environ if environment is None else environment
    for name in ('NVCCFLAGS', 'NVCC_FLAGS', 'NVCC_PREPEND_FLAGS',
                 'NVCC_APPEND_FLAGS', 'RADIA_CUDA_ROOT', 'RADIA_GPU_CUDA_ROOT',
                 'RADIA_MEXCUDA_CONFIG', 'MW_NVCC_PATH', 'CUDA_LIB_PATH',
                 'CUDA_BIN_PATH'):
        if env.get(name, '').strip():
            raise RuntimeError(f'Unreviewed CUDA override: {name}')
    for name in ('CUDA_PATH', 'CUDA_ROOT', 'CUDA_HOME'):
        value = env.get(name, '')
        if value and not re.search(r'[/\\]MATLAB[/\\].*[/\\]sys[/\\]cuda[/\\]win64[/\\]cuda[/\\]?$', value, re.I):
            raise RuntimeError(f'External CUDA selection forbidden: {name}')
    for name, value in env.items():
        if name.startswith('CUDA_PATH_V') and value.strip():
            raise RuntimeError(f'External version-specific CUDA selection: {name}')


def validate_stock_profile(text):
    """Check only profile identity; MathWorks owns its command-template syntax.

    Installed profiles contain literal CMD operators such as && inside quoted
    attributes. They must be passed unchanged to mexcuda, not normalized through
    a generic XML parser. This is a bounded identity check, not an XML validator.
    """
    header = re.search(r'<config\s+([^>]+)>', text, re.S)
    if header is None:
        raise RuntimeError('Installed MathWorks profile header missing')
    attributes = dict(re.findall(r'(\w+)\s*=\s*"([^"]*)"', header.group(1)))
    required = {'Name': 'NVIDIA CUDA Compiler', 'ShortName': 'nvcc',
                'HostCompiler': 'MSVCPP170', 'Language': 'CUDA',
                'Location': r'$CUDA_ROOT\bin\nvcc.exe'}
    if any(attributes.get(key) != value for key, value in required.items()):
        raise RuntimeError('Unacceptable installed MathWorks compiler profile')
    if re.search(r'allow[-_]unsupported[-_]compiler', text, re.I):
        raise RuntimeError('Unsupported compiler override forbidden')
    return attributes


def configure_engine(engine):
    """Use the installed, unmodified MathWorks profile and bundled libraries."""
    reject_overrides()
    vs = discover_vs()
    root = Path(engine.matlabroot()).resolve()
    cuda = root / 'sys/cuda/win64/cuda'
    config = root / 'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
    # Reject redirection outside the selected installed MATLAB tree.
    config.resolve(strict=True).relative_to(root)
    if config.is_symlink():
        raise RuntimeError('Installed compiler profile must not be a redirected file')
    xml = config.read_text(encoding='utf-8')
    profile = validate_stock_profile(xml)
    nvcc = cuda / 'bin/nvcc.exe'
    version = subprocess.check_output([str(nvcc), '--version'], text=True)
    if not re.search(r'release 12\.8\b', version):
        raise RuntimeError('MATLAB bundled CUDA12.8 required')
    files = [config, nvcc, compiler_bin(vs) / 'cl.exe',
             cuda / 'include/crt/host_config.h']
    for name in ('cudart', 'cublas', 'cusparse', 'cusolver'):
        files.append(cuda / f'lib/x64/{name}.lib')
        runtime_major=11 if name=='cusolver' else 12
        files.append(root / f'bin/win64/{name}64_{runtime_major}.dll')
    hashes = {str(file): hashlib.sha256(file.read_bytes()).hexdigest() for file in files}
    engine.setenv('CUDA_PATH', str(cuda), nargout=0)
    engine.setenv('CUDA_ROOT', '', nargout=0)
    engine.setenv('CUDA_HOME', '', nargout=0)
    engine.setenv('CUDA_PATH_V12_8', '', nargout=0)
    engine.setenv('VS170COMNTOOLS', str(vs / 'Common7/Tools'), nargout=0)
    engine.setenv('NVCC_CCBIN', str(compiler_bin(vs)), nargout=0)
    path = str(engine.getenv('PATH'))
    entries = [part for part in path.split(';') if not re.search(
        r'NVIDIA GPU Computing Toolkit|(?:^|[/\\])cuda(?:[-_]?v?\d[\w.-]*)?(?:[/\\]|$)', part, re.I)]
    engine.setenv('PATH', ';'.join([str(root / 'bin/win64'), str(compiler_bin(vs)), *entries]), nargout=0)
    # The stock profile may declare 12.2; actual compiler identity above is decisive.
    receipt = {'msvc': str(compiler_bin(vs)), 'cuda_selection': 'MATLAB bundled only',
               'selector_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'profile_declared_cuda_version': profile.get('Version'),
               'profile_validation': 'Bounded header identity; MathWorks owns command-template validation',
               'nvcc_version': version, 'profile': str(config), 'sha256': hashes}
    return config, receipt
