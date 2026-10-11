"""Build the experimental AMS bridge through Build.ps1 -AmsCudaOnly on 100."""
from pathlib import Path
import hashlib
import json
import os
import socket
import subprocess
import sys
import sysconfig

def main():
    import netgen
    import matlab.engine
    root=Path(__file__).resolve().parents[2]
    if socket.gethostname().lower()!='t7920':
        raise RuntimeError('Candidate builds belong to 100/T7920')
    out=root/'build-ams-cuda';out.mkdir(exist_ok=True)
    net=Path(netgen.__file__).parent
    vs=Path('C:/BuildTools/VS2022-CUDA')
    vcvars=vs/'VC/Auxiliary/Build/vcvars64.bat'
    if not vcvars.is_file():raise RuntimeError('CUDA-compatible VS2022 toolchain required')
    os.environ['NVCC_CCBIN']=str(vs/'VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64')
    os.environ['PATH']=os.environ['NVCC_CCBIN']+os.pathsep+os.environ['PATH']
    source=Path(__file__).with_name('snapshot_module.cpp')
    sources=[source,root/'src/matlab/cuda/radia_cuda_ams.cu',
             root/'src/ext/sparsesolv/include/sparsesolv/preconditioners/compact_ams.hpp',
             root/'src/ext/sparsesolv/include/sparsesolv/preconditioners/compact_amg.hpp',
             root/'src/ext/sparsesolv/include/sparsesolv/ngsolve/ams_accelerator_snapshot.hpp',
             Path(__file__),root/'Build.ps1']
    hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    obj=out/'snapshot.obj';binary=out/'_radia_ams_snapshot.pyd'
    python_lib=Path(sys.base_prefix)/'libs'/f'python{sys.version_info.major}{sys.version_info.minor}.lib'
    command=['cl','/nologo','/LD','/EHsc','/MD','/O2','/fp:precise','/std:c++20','/bigobj','/utf-8','/arch:AVX2',
             '/DWIN32','/DWNT','/DMSVC_EXPRESS','/DMAX_SYS_DIM=3','/DLAPACK','/DNDEBUG',
             '/DNOMINMAX','/DNGS_PYTHON','/DNETGEN_PYTHON','/DNG_PYTHON','/DHAVE_NETGEN_SOURCES',
             '/DHAVE_STRUCT_TIMESPEC','/DPYBIND11_SIMPLE_GIL_MANAGEMENT','/DPy_NO_LINK_LIB',
             '/DSPARSESOLV_USE_NGSOLVE_TASKMANAGER','/D_CRT_SECURE_NO_WARNINGS',
             f'/I{net / "include"}',f'/I{net / "include/include"}',f'/I{sysconfig.get_path("include")}',
             f'/I{root / "src/ext/sparsesolv/include"}',f'/Fo{obj}',str(source),'/link',f'/OUT:{binary}',
             f'/IMPLIB:{out / "_radia_ams_snapshot.lib"}',f'/PDB:{out / "_radia_ams_snapshot.pdb"}',
             str(python_lib),str(net/'lib/libngsolve.lib'),str(net/'lib/ngcore.lib'),str(net/'lib/nglib.lib')]
    batch=out/'compile.cmd'
    batch.write_text('@echo off\ncall "'+str(vcvars)+'" >nul\nif errorlevel 1 exit /b 1\n'+
                     subprocess.list2cmdline(command)+'\nexit /b %errorlevel%\n',encoding='utf-8')
    with (out/'native-build.log').open('w') as log:
        subprocess.run(['cmd','/c',str(batch)],stdout=log,stderr=subprocess.STDOUT,check=True)
    eng=matlab.engine.start_matlab('-nodesktop -nosplash')
    try:
        config=Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        cuda=root/'src/matlab/cuda/radia_cuda_ams.cu'
        flags='NVCCFLAGS=-gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80'
        eng.mexcuda('-v','-f',str(config),flags,'-lcublas','-outdir',str(out),str(cuda),nargout=0)
        if hashes!={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}:
            raise RuntimeError('Source changed while building; candidate provenance rejected')
        record={'schema':'radia.ams-cuda-build.v1','host':socket.gethostname(),'python':sys.executable,
                'matlab':eng.version(),'configuration':str(config),'compiler_bin':os.environ['NVCC_CCBIN'],
                'source_commit':subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
                'source_dirty':bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip()),
                'sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (binary,out/'radia_cuda_ams.mexw64')},
                'source_sha256':hashes,'no_cpu_solver_fallback':True}
        (out/'build.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps(record,indent=2))
    finally:eng.quit()

if __name__=='__main__':main()
