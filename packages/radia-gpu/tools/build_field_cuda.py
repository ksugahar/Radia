"""Stage the experimental RadFld cloud CUDA module on the designated build host."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import socket
import subprocess
import tempfile
import ctypes


def main():
    from vs2022_toolchain import discover_vs, compiler_bin as supported_compiler_bin
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-clean-source', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    if socket.gethostname().lower() != 't7920':
        raise RuntimeError('Candidate CUDA builds belong to 100')
    commit = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    dirty = bool(subprocess.check_output(['git', '-C', str(root), 'status', '--porcelain'], text=True).strip())
    if args.require_clean_source and dirty:
        raise RuntimeError('Clean source required')
    core = root / 'src/core'
    package = root / 'packages/radia-gpu'
    sources = [core/'rad_hdiv_field_evaluator.cpp', core/'rad_hdiv_vim.cpp',
               core/'rad_hdiv_field_evaluator.h', core/'rad_hdiv_vim.h', core/'rad_field_cloud_kernel.h',
               core/'rad_field_sources.h', core/'rad_field_source_kernel.inc', core/'rad_hdiv_field_math.inc',
               core/'rad_hdiv_field_kernel_api.inc',
               package/'native/radia_cuda_field.cu', package/'native/radia_field_exact.cuh', Path(__file__).with_name('vs2022_toolchain.py'),Path(__file__), root/'Build.ps1']
    hashes = {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    base = root/'build-field-cuda'
    base.mkdir(exist_ok=True)
    out = Path(tempfile.mkdtemp(prefix='job-',dir=base))
    print(str(out),flush=True)
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
    vs = discover_vs()
    compiler_bin = supported_compiler_bin(vs)
    os.environ['NVCC_CCBIN'] = str(compiler_bin)
    os.environ['PATH'] = str(compiler_bin)+os.pathsep+os.environ['PATH']
    import matlab.engine
    eng = matlab.engine.start_matlab('-nodesktop -nosplash')
    try:
        pid=int(eng.feature('getpid'))
        subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        config = Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        if not config.is_file():
            raise RuntimeError('Supported MATLAB CUDA compiler configuration missing')
        objects = []
        for source in sources[:2]:
            obj = out/(source.stem+'.obj')
            command = ['cl', '/nologo', '/c', '/O2', '/MD', '/EHsc', '/std:c++17', '/fp:precise',
                       '/DRADIA_FIELD_SERIAL_ONLY', '/DNOMINMAX', f'/I{core}', f'/Fo{obj}', str(source)]
            batch = out/(source.stem+'.cmd')
            batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+
                             '" >nul\nif errorlevel 1 exit /b 1\n'+subprocess.list2cmdline(command)+
                             '\nexit /b %errorlevel%\n', encoding='utf-8')
            with (out/(source.stem+'.log')).open('w') as log:
                subprocess.run(['cmd', '/c', str(batch)], stdout=log, stderr=subprocess.STDOUT, check=True)
            objects.append(str(obj))
        # No fast math: match the CPU operation order, including compensated sums.
        flags = ('NVCCFLAGS=--fmad=false -gencode=arch=compute_75,code=sm_75 '
                 '-gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80')
        eng.mexcuda('-v', '-f', str(config), flags, f'-I{core}', '-outdir', str(out),
                    str(package/'native/radia_cuda_field.cu'), *objects, nargout=0)
        if hashes != {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}:
            raise RuntimeError('Source changed during build')
        if subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip() != commit:
            raise RuntimeError('Source commit changed during build')
        binary = out/'radia_cuda_field.mexw64'
        record = dict(schema='radia.field-cuda-build.v1', source_commit=commit, source_dirty=dirty,
                      source_sha256=hashes, binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                      matlab=eng.version(), cuda_targets=['sm_75', 'sm_80', 'compute_80 PTX'],
                      source_representation='existing retained point cloud only',
                      cpu_fallback=False, python_runtime_required=False, taskmanager_regions=False)
        (out/'build.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(record, indent=2))
    finally:
        eng.quit()


if __name__ == '__main__':
    main()
