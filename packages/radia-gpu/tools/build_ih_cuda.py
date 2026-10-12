"""Isolated, low-priority persistent IH mexcuda build on 100."""
from pathlib import Path
import argparse,ctypes,hashlib,json,socket,subprocess

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--require-clean-source',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parents[3];package=root/'packages/radia-gpu';out=root/'build-ih-cuda';out.mkdir(exist_ok=True)
    record=dict(schema='radia.ih-cuda-build.v1',accepted=False,scope='Persistent fixed operators; dynamic ESIM physical frontend pending',a100_used=False)
    eng=None
    try:
        if socket.gethostname().lower()!='t7920':raise RuntimeError('IH candidate builds belong to 100')
        record['source_commit']=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
        record['source_dirty']=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip())
        if args.require_clean_source and record['source_dirty']:raise RuntimeError('Clean source required')
        sources=[package/'native'/name for name in ('radia_cuda_ih.cu','ih_coupled_state.cuh','ih_complex_solver.cuh','ih_bicgstab_kernels.cuh','ih_thermal_kernels.cuh','ih_resident_kernels.cuh','ih_preconditioner.cuh')]+[Path(__file__),root/'Build.ps1']
        hashes=lambda:{p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        record['source_sha256']=hashes()
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash');pid=int(eng.feature('getpid'))
        subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        compiler='C:/BuildTools/VS2022-CUDA/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64'
        eng.setenv('NVCC_CCBIN',compiler,nargout=0);eng.setenv('PATH',compiler+';'+eng.getenv('PATH'),nargout=0)
        config=Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        flags='NVCCFLAGS=--fmad=false -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80'
        eng.mexcuda('-v','-R2018a','-f',str(config),flags,'-I'+str(package/'native'),'-lcublas','-outdir',str(out),str(sources[0]),nargout=0)
        if hashes()!=record['source_sha256'] or subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()!=record['source_commit']:raise RuntimeError('Source changed during build')
        record.update(accepted=True,matlab=eng.version(),binary_sha256=hashlib.sha256((out/'radia_cuda_ih.mexw64').read_bytes()).hexdigest(),cuda_targets=['sm_75','sm_80','compute_80 PTX'],cpu_fallback=False,python_runtime_required=False,taskmanager_regions=False)
    except BaseException as error:record['error']=str(error);raise
    finally:
        if eng is not None:eng.quit()
        (out/'build.json').write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
