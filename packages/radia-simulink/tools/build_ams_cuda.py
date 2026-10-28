"""Build the experimental AMS bridge through Build.ps1 -AmsCudaOnly on an explicitly admitted build host."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import socket
import subprocess
import sys
import sysconfig
import ctypes

def main(args):
    from matlab_cuda_toolchain import require_build_admission
    from check_mex_error_boundaries import require_boundaries
    require_boundaries()
    from vs2022_toolchain import discover_vs, compiler_bin as supported_compiler_bin
    import netgen
    import matlab.engine
    package=Path(__file__).resolve().parents[1]
    root=package.parents[1]
    commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
    dirty=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip())
    if args.require_clean_source and dirty:raise RuntimeError('Clean source required for AMS build provenance')
    admission=require_build_admission(root,'ams',args.build_admission,args.build_admission_sha256)
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
    out=root/'build-ams-cuda';out.mkdir(exist_ok=True)
    net=Path(netgen.__file__).parent
    vs=discover_vs()
    vcvars=vs/'VC/Auxiliary/Build/vcvars64.bat'
    if not vcvars.is_file():raise RuntimeError('CUDA-compatible VS2022 toolchain required')
    os.environ['NVCC_CCBIN']=str(supported_compiler_bin(vs))
    os.environ['PATH']=os.environ['NVCC_CCBIN']+os.pathsep+os.environ['PATH']
    source=package/'native/ams_snapshot_module.cpp'
    sources=[source,package/'native/ams_cpu_reference.cpp',package/'native/ams_resident.cuh',package/'native/ams_nonlinear.cuh',package/'native/radia_cuda_ams.cu',package/'native/ams_matlab_setup.cpp',
             root/'src/core/gpu/ams_amg_kernels.cuh',root/'src/core/gpu/ams_core.cuh',
             root/'src/core/gpu/snapshot.hpp',package/'native/ams_setup_gpu.hpp',package/'native/ams_setup_gpu.cuh',
             package/'native/ams_fem_host.hpp',package/'native/ams_fem_host.cpp',package/'native/ams_fem.cuh',
             root/'src/ext/sparsesolv/include/sparsesolv/preconditioners/compact_ams.hpp',
             root/'src/ext/sparsesolv/include/sparsesolv/preconditioners/compact_amg.hpp',
             root/'src/ext/sparsesolv/include/sparsesolv/preconditioners/compact_parallel.hpp',
             root/'src/ext/sparsesolv/include/sparsesolv/ngsolve/ams_accelerator_snapshot.hpp',
             Path(__file__).with_name('matlab_cuda_toolchain.py'),Path(__file__).with_name('vs2022_toolchain.py'),Path(__file__),root/'Build.ps1']
    sources.append(root/'packages/radia-simulink/native/mex_error_boundary.hpp')
    sources.append(Path(__file__).with_name('check_mex_error_boundaries.py'))
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    obj=out/'snapshot.obj';binary=out/'_radia_ams_snapshot.pyd'
    python_lib=Path(sys.base_prefix)/'libs'/f'python{sys.version_info.major}{sys.version_info.minor}.lib'
    command=['cl','/nologo','/LD','/EHsc','/MD','/O2','/fp:precise','/std:c++20','/bigobj','/utf-8','/arch:AVX2',
             '/DWIN32','/DWNT','/DMSVC_EXPRESS','/DMAX_SYS_DIM=3','/DLAPACK','/DNDEBUG',
             '/DNOMINMAX','/DNGS_PYTHON','/DNETGEN_PYTHON','/DNG_PYTHON','/DHAVE_NETGEN_SOURCES',
             '/DHAVE_STRUCT_TIMESPEC','/DPYBIND11_SIMPLE_GIL_MANAGEMENT','/DPy_NO_LINK_LIB',
             '/DSPARSESOLV_SERIAL_AMS_SETUP','/DSPARSESOLV_OPENMP_AMS_SETUP','/openmp','/D_CRT_SECURE_NO_WARNINGS',
             f'/I{net / "include"}',f'/I{net / "include/include"}',f'/I{sysconfig.get_path("include")}',
             f'/I{root / "src/ext/sparsesolv/include"}',f'/Fo{obj}',str(source),'/link',f'/OUT:{binary}',
             f'/IMPLIB:{out / "_radia_ams_snapshot.lib"}',f'/PDB:{out / "_radia_ams_snapshot.pdb"}',
             str(python_lib),str(net/'lib/libngsolve.lib'),str(net/'lib/ngcore.lib'),str(net/'lib/nglib.lib')]
    batch=out/'compile.cmd'
    batch.write_text('@echo off\ncall "'+str(vcvars)+'" >nul\nif errorlevel 1 exit /b 1\n'+
                     subprocess.list2cmdline(command)+'\nexit /b %errorlevel%\n',encoding='utf-8')
    with (out/'native-build.log').open('w') as log:
        subprocess.run(['cmd','/c',str(batch)],stdout=log,stderr=subprocess.STDOUT,check=True)
    cpu_source=package/'native/ams_cpu_reference.cpp'
    cpu_command=[arg.replace(str(source),str(cpu_source)).replace('_radia_ams_snapshot','_radia_ams_cpu_reference').replace('snapshot.obj','cpu_reference.obj')
                 for arg in command if arg not in ('/DSPARSESOLV_SERIAL_AMS_SETUP','/DSPARSESOLV_OPENMP_AMS_SETUP','/openmp')]
    cpu_batch=out/'compile-cpu-reference.cmd'
    cpu_batch.write_text('@echo off\ncall "'+str(vcvars)+'" >nul\nif errorlevel 1 exit /b 1\n'+
                        subprocess.list2cmdline(cpu_command)+'\nexit /b %errorlevel%\n',encoding='utf-8')
    with (out/'cpu-reference-build.log').open('w') as log:
        subprocess.run(['cmd','/c',str(cpu_batch)],stdout=log,stderr=subprocess.STDOUT,check=True)
    eng=matlab.engine.start_matlab('-nodesktop -nosplash')
    try:
        pid=int(eng.feature('getpid'))
        subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        config=Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        cuda=package/'native/radia_cuda_ams.cu'
        # Compile shared NGSolve host code with the same supported MSVC as
        # mexcuda; do not send NGSolve/Python headers through the CUDA parser.
        host_obj=out/'ams_matlab_setup.obj'
        host_command=command[:command.index('/link')]
        host_command=[arg for arg in host_command if arg not in ('/LD',str(source)) and not arg.startswith('/Fo')]
        host_command+=['/c','/DMX_COMPAT_64','/DMATLAB_DEFAULT_RELEASE=R2017b','/DMATLAB_MEX_FILE',
                       f'/I{Path(eng.matlabroot()) / "extern/include"}',
                       f'/Fo{host_obj}',str(package/'native/ams_matlab_setup.cpp')]
        batch.write_text('@echo off\ncall "'+str(vcvars)+'" >nul\nif errorlevel 1 exit /b 1\n'+
                         subprocess.list2cmdline(host_command)+'\nexit /b %errorlevel%\n',encoding='utf-8')
        with (out/'host-build.log').open('w') as log:
            subprocess.run(['cmd','/c',str(batch)],stdout=log,stderr=subprocess.STDOUT,check=True)
        fem_obj=out/'ams_fem_host.obj'
        fem_command=[arg.replace(str(host_obj),str(fem_obj)).replace('ams_matlab_setup.cpp','ams_fem_host.cpp') for arg in host_command]
        batch.write_text('@echo off\ncall "'+str(vcvars)+'" >nul\nif errorlevel 1 exit /b 1\n'+
                         subprocess.list2cmdline(fem_command)+'\nexit /b %errorlevel%\n',encoding='utf-8')
        with (out/'fem-host-build.log').open('w') as log:
            subprocess.run(['cmd','/c',str(batch)],stdout=log,stderr=subprocess.STDOUT,check=True)
        flags='NVCCFLAGS=-gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80'
        eng.mexcuda('-v','-f',str(config),flags,'-lcublas','-outdir',str(out),str(cuda),str(host_obj),str(fem_obj),
                    str(python_lib),str(net/'lib/libngsolve.lib'),str(net/'lib/ngcore.lib'),str(net/'lib/nglib.lib'),nargout=0)
        if (hashes!={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
            or commit!=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()):
            raise RuntimeError('Source changed while building; candidate provenance rejected')
        record={'build_admission':admission,'schema':'radia.ams-cuda-build.v1','host':socket.gethostname(),'python':sys.executable,
                'matlab':eng.version(),'configuration':str(config),'compiler_bin':os.environ['NVCC_CCBIN'],
                'source_commit':commit,'source_dirty':dirty,
                'sha256':{p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (binary,out/'_radia_ams_cpu_reference.pyd',out/'radia_cuda_ams.mexw64')},
                'source_sha256':hashes,'no_cpu_solver_fallback':True}
        record['snapshot_setup_backend']='OpenMP row/chunk dispatch; default one thread; serial coarse factorization'
        record['matlab_setup_backend']='Shared native AMS; optional CUDA AP and Pt(AP); serial NGSolve sparsecholesky; no TaskManager regions'
        record['python_runtime_adapter_required']=False
        record['python_dll_indirect_dependency_permitted']=True
        (out/'build.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps(record,indent=2))
    finally:eng.quit()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--require-clean-source',action='store_true')
    p.add_argument('--build-admission',type=Path);p.add_argument('--build-admission-sha256')
    main(p.parse_args())
