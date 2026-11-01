"""One package-owned nvcc archive, Python first, then plain CMake/MSVC MEX.

No installation, device context, fallback toolchain or numerical validation.
Explicit verification receipt from validate_kernels.py is required before MEX.
"""
import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
GPU_SOURCE=ROOT/'native/optuna_cuda.cu'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def reject_compiler_overrides(environment):
    for key in ['CUDA_PATH','CUDA_BIN_PATH','CUDA_ROOT','CUDA_LIB_PATH','MW_NVCC_PATH',
                'MW_ALLOW_ANY_CUDA','MW_ALLOW_ANY_COMPILER','NVCC_CCBIN',
                'NVCC_PREPEND_FLAGS','NVCC_APPEND_FLAGS','NVCC_FLAGS',
                'CL','_CL_','CFLAGS','CXXFLAGS','CUDAFLAGS']:
        if environment.get(key):raise RuntimeError('Unsupported compiler override: '+key)
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matlab-root',type=Path,required=True)
    p.add_argument('--host-compiler',type=Path,required=True)
    p.add_argument('--build-dir',type=Path,required=True)
    p.add_argument('--stage',choices=['python','matlab'],required=True)
    p.add_argument('--python-verification',type=Path)
    a=p.parse_args();matlab=a.matlab_root.resolve();cl=a.host_compiler.resolve()
    source_sha=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
    if subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip():
        raise RuntimeError('Native build requires a clean exact source checkout')
    nvcc=matlab/'sys/cuda/win64/cuda/bin/nvcc.exe'
    if os.name!='nt' or not nvcc.is_file() or not cl.is_file():raise RuntimeError('Explicit Windows baseline paths required')
    if not re.search(r'/VC/Tools/MSVC/14[.]44[.][^/]+/bin/HostX64/x64/cl[.]exe$',cl.as_posix(),re.I):
        raise RuntimeError('MSVC14.44 x64 host path required')
    if metadata.version('pybind11')!='3.0.2':raise RuntimeError('pybind11==3.0.2 required by reviewed CI baseline')
    reject_compiler_overrides(os.environ)
    version=subprocess.check_output([str(nvcc),'--version'],text=True)
    if not re.search(r'release 12\.8[,\s]',version):raise RuntimeError('CUDA12.8 required')
    build=a.build_dir.resolve();build.mkdir(parents=True,exist_ok=True)
    native_sources={f.relative_to(ROOT).as_posix():sha(f) for f in sorted((ROOT/'native').glob('*')) if f.is_file()}
    recipes={name:sha(ROOT/name) for name in ['CMakeLists.txt','pyproject.toml','tools/build_native.py','tools/build_native.ps1']}
    cuda=nvcc.parent.parent;cuda_object=build/'optuna_cuda.obj'
    nvcc_command=[str(nvcc),'-v','-c','-std=c++20','--fmad=false',
        '--compiler-options=/Zp8,/GR,/W3,/EHs,/nologo,/MD,/O2,/Oy-,/DNDEBUG,/openmp',
        '-ccbin',str(cl),'--dont-use-profile','-ldir',str(matlab/'bin/win64'),
        '-I'+str(cuda/'include'),'-DNOMINMAX','-D_CRT_SECURE_NO_DEPRECATE',
        '-D_SCL_SECURE_NO_DEPRECATE','-D_SECURE_SCL=0',
        '-gencode=arch=compute_75,code=sm_75','-gencode=arch=compute_80,code=sm_80',
        '-gencode=arch=compute_80,code=compute_80','-o',str(cuda_object),str(GPU_SOURCE)]
    previous=None
    if a.stage=='matlab':
        if not a.python_verification:raise RuntimeError('Python-first verification receipt required')
        previous=json.loads((build/'python-build.json').read_text())
        proof=json.loads(a.python_verification.read_text())
        if (not proof['ok'] or proof['source_hashes']!=native_sources or
            proof['archive_sha256']!=previous['archive_sha256'] or
            previous['recipe_hashes']!=recipes or previous['source_sha']!=source_sha or
            proof.get('source_sha')!=source_sha or proof.get('recipe_hashes')!=recipes or
            proof.get('oracle_sha256')!=sha(ROOT/'oracles/optuna50_cuda_oracle.json') or
            previous['nvcc_command']!=nvcc_command or
            previous['objects'].get(str(cuda_object))!=sha(cuda_object)):
            raise RuntimeError('Python verification is missing, failed or from another source/archive')
    pybind_dir=subprocess.check_output([sys.executable,'-m','pybind11','--cmakedir'],text=True).strip()
    configure=['cmake','-S',str(ROOT),'-B',str(build),'-G','Ninja',
      '-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_COMPILER='+cl.as_posix(),
      '-DOPTUNA_CUDA_OBJECT='+cuda_object.as_posix(),'-DOPTUNA_CUDA_ROOT='+cuda.as_posix(),
      '-DPython3_EXECUTABLE='+Path(sys.executable).as_posix(),'-Dpybind11_DIR='+pybind_dir,
      '-DMatlab_ROOT_DIR='+matlab.as_posix(),'-DOPTUNA_BUILD_MATLAB='+('ON' if a.stage=='matlab' else 'OFF')]
    commands=[configure,['cmake','--build',str(build),'--parallel','2','--target',
                          'optuna_mex' if a.stage=='matlab' else '_optuna_kernels']]
    if a.stage=='python':commands.insert(0,nvcc_command)
    compiler_env=dict(os.environ)
    compiler_env['PATH']=str(cuda/'nvvm/bin')+os.pathsep+str(cuda/'bin')+os.pathsep+compiler_env['PATH']
    for backend in [cuda/'nvvm/bin/cicc.exe',cuda/'bin/ptxas.exe',cuda/'bin/fatbinary.exe']:
        if not backend.is_file():raise RuntimeError('Missing explicitly selected CUDA compiler backend: '+str(backend))
    for index,cmd in enumerate(commands):
        with (build/f'{a.stage}-{index}.log').open('w',encoding='utf-8') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,env=compiler_env)
    objects=[Path(x) for x in (build/'object-list.txt').read_text().splitlines() if x]
    archive=build/'optuna_kernels.lib'
    link_commands=subprocess.check_output(['ninja','-C',str(build),'-t','commands'],text=True)
    (build/f'{a.stage}-link-commands.txt').write_text(link_commands,encoding='utf-8',newline='\n')
    receipt={'schema':'radia.optuna.native-build.v1','stage':a.stage,'python':sys.version,
             'source_sha':source_sha,'recipe_hashes':recipes,'nvcc_command':nvcc_command,
             'libdevice_sha256':sha(matlab/'bin/win64/libdevice.10.bc'),
             'compiler':str(cl),'nvcc':str(nvcc),'nvcc_version':version,'source_hashes':native_sources,
             'archive':str(archive),'archive_sha256':sha(archive),
             'objects':{str(f):sha(f) for f in objects},'commands':commands,
             'compile_commands_sha256':sha(build/'compile_commands.json'),
             'qualification':'build only; numerical/ABI/lifecycle parity separate'}
    if previous and (previous['archive_sha256']!=receipt['archive_sha256'] or previous['objects']!=receipt['objects']):
        raise RuntimeError('MEX did not consume the identical Python-verified archive/objects')
    if a.stage=='matlab':
        receipt['python_verification_sha256']=sha(a.python_verification)
        mex=ROOT/'matlab/optuna_mex.mexw64';receipt['mex_sha256']=sha(mex)
    receipt['python_extension_sha256s']={f.name:sha(f) for f in (ROOT/'src/radia_optuna').glob('_optuna_kernels*.pyd')}
    (build/f'{a.stage}-build.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8',newline='\n')
    if a.stage=='matlab':
        # Public wheel provenance must not expose private host/checkout paths.
        public={k:receipt[k] for k in ['schema','stage','source_sha','source_hashes','recipe_hashes','libdevice_sha256','archive_sha256',
                'mex_sha256','python_extension_sha256s','python_verification_sha256']}
        public['objects']={Path(f).name:h for f,h in receipt['objects'].items()}
        public['build_route']='nvcc-shared-library'
        public['wrapper']='plain-msvc-mex'
        (ROOT/'src/radia_optuna/native_build.json').write_text(json.dumps(public,indent=2),encoding='utf-8',newline='\n')
if __name__=='__main__':main()
