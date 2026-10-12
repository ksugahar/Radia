"""Isolated resident HACApK algebra CUDA build, not a physical solver release."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,socket,subprocess,tempfile


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--require-clean-source',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parents[3]
    if socket.gethostname().lower()!='t7920':raise RuntimeError('CUDA build belongs to 100')
    dirty=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip())
    if args.require_clean_source and dirty:raise RuntimeError('Clean source required')
    files=[root/'packages/radia-gpu/native/radia_cuda_hmatrix.cu',root/'packages/radia-gpu/native/hmatrix_resident.cuh',root/'packages/radia-gpu/native/hmatrix_hacapk_pack.cuh',root/'packages/radia-gpu/native/hmatrix_entry_oracle.cuh',root/'packages/radia-gpu/native/hmatrix_aca_callback.cuh',root/'src/ext/HACApK/cHACApK_calc_entry_ij.h',root/'src/ext/HACApK/cHACApK_base.h',root/'packages/radia-gpu/native/analytic_matrix_build.cuh',root/'packages/radia-gpu/native/radia_field_exact.cuh',root/'src/core/rad_charge_quadrature.h',root/'src/core/rad_hdiv_vim.h',root/'src/core/rad_field_cloud_kernel.h',root/'src/core/rad_hdiv_field_math.inc',root/'src/core/rad_hdiv_field_kernel_api.inc',root/'src/core/rad_field_source_kernel.inc',Path(__file__),root/'Build.ps1']
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    base=root/'build-hmatrix-cuda';base.mkdir(exist_ok=True);out=Path(tempfile.mkdtemp(prefix='job-',dir=base));print(out,flush=True)
    record=dict(schema='radia.hmatrix-cuda-build.v1',accepted=False,source_dirty=dirty,source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),source_sha256=hashes,a100_used=False,scope='Resident numeric HACApK leaf application only; assembly/physical solver not accepted')
    eng=None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        cc=Path('C:/BuildTools/VS2022-CUDA/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64');os.environ['NVCC_CCBIN']=str(cc)
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
        pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        eng.setenv('NVCC_CCBIN',str(cc),nargout=0);eng.setenv('PATH',str(cc)+';'+eng.getenv('PATH'),nargout=0)
        config=Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        if not config.is_file():raise RuntimeError('Supported MATLAB CUDA configuration missing')
        eng.mexcuda('-v','-R2018a','-f',str(config),'NVCCFLAGS=-std=c++17 --fmad=false -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80','-I'+str(root/'src/ext/HACApK'),'-I'+str(root/'src/core'),'-outdir',str(out),str(files[0]),nargout=0)
        if hashes!={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}:raise RuntimeError('Source changed during build')
        if subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()!=record['source_commit']:raise RuntimeError('Commit changed during build')
        record.update(accepted=True,matlab=eng.version(),binary_sha256=hashlib.sha256((out/'radia_cuda_hmatrix.mexw64').read_bytes()).hexdigest(),cuda_targets=['sm_75','sm_80','compute_80 PTX'],taskmanager_started=False,python_runtime_required=False)
    except BaseException as error:record['error']=f'{type(error).__name__}: {error}';raise
    finally:
        if eng is not None:eng.quit()
        (out/'build.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
