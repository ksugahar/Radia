"""Build the isolated numeric Laplace FMM through Build.ps1 -FmmCudaOnly."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,socket,subprocess,tempfile,sys

def main():
    from vs2022_toolchain import discover_vs, compiler_bin as supported_compiler_bin
    parser=argparse.ArgumentParser();parser.add_argument('--require-clean-source',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parents[3]
    if socket.gethostname().lower()!='t7920':raise RuntimeError('Candidate builds belong to 100')
    dirty=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip())
    if args.require_clean_source and dirty:raise RuntimeError('Clean source required')
    source=root/'packages/radia-gpu/native/radia_cuda_fmm.cu'
    sources=[source,source.with_name('fmm_interpolation.cuh'),source.with_name('fmm_boundary_host.cpp'),root/'src/core/rad_bem_galerkin.cpp',root/'src/core/rad_bem_galerkin.h',root/'src/core/rad_host_parallel.h',Path(__file__).with_name('vs2022_toolchain.py'),Path(__file__),root/'Build.ps1']
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    base=root/'build-fmm-cuda';base.mkdir(exist_ok=True);out=Path(tempfile.mkdtemp(prefix='job-',dir=base));print(out,flush=True)
    receipt=dict(schema='radia.fmm-cuda-build.v1',accepted=False,source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),source_dirty=dirty,source_sha256=hashes,python_adapter_required=False,taskmanager_used=False)
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
    vs=discover_vs();cc=supported_compiler_bin(vs)
    os.environ['NVCC_CCBIN']=str(cc);os.environ['PATH']=str(cc)+os.pathsep+os.environ['PATH']
    engine=None
    try:
        import matlab.engine
        engine=matlab.engine.start_matlab('-nodesktop -nosplash')
        pid=int(engine.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],capture_output=True,check=True)
        config=Path(engine.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        if not config.is_file():raise RuntimeError('Supported MATLAB CUDA compiler configuration missing')
        import netgen
        net=Path(netgen.__file__).parent;prefix=Path(sys.base_prefix);objects=[]
        for host in [source.with_name('fmm_boundary_host.cpp'),root/'src/core/rad_bem_galerkin.cpp']:
            obj=out/(host.stem+'.obj');objects.append(str(obj))
            includes=[net/'include',net/'include/include',prefix/'include',root/'src/core',Path(engine.matlabroot())/'extern/include']
            command=['cl','/nologo','/c','/O2','/MD','/EHsc','/std:c++20','/DWIN32','/DNOMINMAX','/DRAD_MATLAB_HOST','/DMATLAB_MEXCMD_RELEASE=R2018a','/DUSE_MEX_CMD','/DMATLAB_MEX_FILE',*[f'/I{p}' for p in includes],f'/Fo{obj}',str(host)]
            batch=out/(host.stem+'.cmd');batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\nif errorlevel 1 exit /b 1\n'+subprocess.list2cmdline(command)+'\nexit /b %errorlevel%\n',encoding='utf-8')
            with (out/(host.stem+'.log')).open('w') as log:subprocess.run(['cmd','/c',str(batch)],stdout=log,stderr=subprocess.STDOUT,check=True)
        engine.mexcuda('-v','-f',str(config),'-R2018a','NVCCFLAGS=-std=c++17 --fmad=false -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80','-outdir',str(out),str(source),*objects,str(net/'lib/libngsolve.lib'),str(net/'lib/ngcore.lib'),str(net/'lib/nglib.lib'),str(prefix/'libs/python312.lib'),nargout=0)
        if hashes!={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}:raise RuntimeError('Sources changed during build')
        receipt.update(accepted=True,binary_sha256=hashlib.sha256((out/'radia_cuda_fmm.mexw64').read_bytes()).hexdigest(),matlab=engine.version(),compiler_bin=str(cc),configuration=str(config))
    except BaseException as e:receipt['error']=repr(e);raise
    finally:
        if engine is not None:engine.quit()
        (out/'build.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
