"""Staged standalone AMG pybind module; invoked through Build.ps1 only."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import sysconfig

def build(package,out,admission):
    import matlab.engine
    import netgen
    from matlab_cuda_toolchain import configure_engine
    from vs2022_toolchain import discover_vs,compiler_bin
    out.mkdir(parents=True,exist_ok=True)
    native=package/'native'
    repo=package.parents[1]
    commit=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
    def require_clean():
        current=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
        dirty=subprocess.check_output(['git','-C',str(repo),'status','--porcelain'],text=True).strip()
        if current!=commit or dirty:
            raise RuntimeError('Resident AMG requires unchanged clean exact-source build')
    require_clean()
    engine=matlab.engine.start_matlab('-nodesktop -nosplash')
    try:
        config,toolchain=configure_engine(engine)
        flags='NVCCFLAGS=-gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80'
        engine.mexcuda('-v','-c','-f',str(config),flags,'-outdir',str(out),str(native/'resident_amg_cuda.cu'),nargout=0)
        obj=out/'resident_amg_cuda.obj';vs=discover_vs();compiler=compiler_bin(vs)
        binary=out/('_resident_amg'+sysconfig.get_config_var('EXT_SUFFIX'))
        root=Path(engine.matlabroot());include=Path(netgen.__file__).parent/'include'
        if not (include/'pybind11/pybind11.h').is_file():raise RuntimeError('NGSolve-matched pybind headers missing')
        command=[str(compiler/'cl.exe'),'/nologo','/LD','/EHsc','/MD','/O2','/fp:precise','/std:c++20','/utf-8',f'/I{include}',f'/I{sysconfig.get_path("include")}',f'/Fo{out/"resident_amg_python.obj"}',str(native/'resident_amg_python.cpp'),str(obj),'/link',f'/OUT:{binary}',f'/IMPLIB:{out/"resident_amg.lib"}',str(Path(sys.base_prefix)/'libs'/f'python{sys.version_info.major}{sys.version_info.minor}.lib')]+[str(root/f'sys/cuda/win64/cuda/lib/x64/{name}.lib') for name in ['cudart','cublas']]
        batch=out/'compile-resident-amg.cmd';batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\nif errorlevel 1 exit /b 1\n'+subprocess.list2cmdline(command)+'\nexit /b %errorlevel%\n')
        with (out/'resident-amg-host-build.log').open('w') as log:subprocess.run(['cmd','/c',str(batch)],stdout=log,stderr=subprocess.STDOUT,check=True)
        sources=[native/f'resident_amg_{suffix}' for suffix in ['cuda.cu','python.cpp','host.hpp']]
        sources += [Path(__file__), Path(__file__).with_name('build_cuda_mex.py'),
                    Path(__file__).with_name('matlab_cuda_toolchain.py'),
                    Path(__file__).with_name('vs2022_toolchain.py'), repo/'Build.ps1']
        require_clean()
        record=dict(schema='radia.resident-amg-python-build.v1', admission=admission,
                    toolchain=toolchain,source_commit=commit,source_dirty=False,
                    binary=binary.name,extension_suffix=sysconfig.get_config_var('EXT_SUFFIX'),
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    source_hashes={f.relative_to(repo).as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in sources},
                    matlab=engine.version(),python=sys.version,engine_runtime_required=False,
                    matlab_engine_user_dependency=False,fault_injection=False,
                    device_fault_injection=False,lifecycle_diagnostic=False,
                    runtime_import_checked=False,redistribution_accepted=False)
        (out/'resident-amg-build.json').write_text(json.dumps(record,indent=2)+'\n')
    finally:engine.quit()
