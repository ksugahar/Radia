"""Isolated resident HACApK algebra CUDA build, not a physical solver release."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,socket,subprocess,tempfile


def main():
    from vs2022_toolchain import discover_vs, compiler_bin as supported_compiler_bin
    parser=argparse.ArgumentParser();parser.add_argument('--require-clean-source',action='store_true');parser.add_argument('--mkl-root',type=Path);args=parser.parse_args()
    root=Path(__file__).resolve().parents[3]
    if socket.gethostname().lower()!='t7920':raise RuntimeError('CUDA build belongs to 100')
    dirty=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip())
    if args.require_clean_source and dirty:raise RuntimeError('Clean source required')
    files=[root/'packages/radia-gpu/native/radia_cuda_hmatrix.cu',root/'packages/radia-gpu/native/hmatrix_resident.cuh',root/'packages/radia-gpu/native/hmatrix_hacapk_pack.cuh',root/'packages/radia-gpu/native/hmatrix_entry_oracle.cuh',root/'packages/radia-gpu/native/hmatrix_aca_callback.cuh',root/'src/ext/HACApK/cHACApK_calc_entry_ij.h',root/'src/ext/HACApK/cHACApK_base.h',root/'packages/radia-gpu/native/analytic_matrix_build.cuh',root/'packages/radia-gpu/native/radia_field_exact.cuh',root/'src/core/rad_charge_quadrature.h',root/'src/core/rad_hdiv_vim.h',root/'src/core/rad_field_cloud_kernel.h',root/'src/core/rad_hdiv_field_math.inc',root/'src/core/rad_hdiv_field_kernel_api.inc',root/'src/core/rad_field_source_kernel.inc',Path(__file__).with_name('vs2022_toolchain.py'),Path(__file__),root/'Build.ps1']
    files += [root/'src/ext/HACApK/cHACApK_base.c',root/'src/ext/HACApK/cHACApK_aca.inc',root/'packages/radia-gpu/native/hmatrix_aca_host.c',root/'src/ext/HACApK/cHACApK_lib.c',root/'src/core/rad_hacapk_parallel.h'] + sorted((root/'src/ext/HACApK').glob('*.h'))
    files += [root/'packages/radia-gpu/native/hmatrix_material_solver.cuh',root/'packages/radia-gpu/native/hdiv_curved_gram.cuh']
    if args.mkl_root:
        files += [root/'packages/radia-gpu/native/hmatrix_full_assembly.cuh',root/'packages/radia-gpu/native/hmatrix_serial_host.c',root/'src/ext/HACApK/cHACApK_cpp_impl.c',root/'src/ext/HACApK/cHACApK_pca_cluster.c']
    files=list(dict.fromkeys(files))
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    base=root/'build-hmatrix-cuda';base.mkdir(exist_ok=True);out=Path(tempfile.mkdtemp(prefix='job-',dir=base));print(out,flush=True)
    record=dict(schema='radia.hmatrix-cuda-build.v1',accepted=False,source_dirty=dirty,source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),source_sha256=hashes,a100_used=False,scope='Resident numeric HACApK leaf application only; assembly/physical solver not accepted')
    eng=None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        vs=discover_vs();cc=supported_compiler_bin(vs);os.environ['NVCC_CCBIN']=str(cc)
        objects=[]
        full=bool(args.mkl_root);extra=[]
        if full:
            if not (args.mkl_root/'include/mkl.h').is_file() or not (args.mkl_root/'lib/mkl_rt.lib').is_file():raise RuntimeError('Explicit supported mkl-devel root required')
            os.environ['MKL_NUM_THREADS']='1';os.environ['MKL_THREADING_LAYER']='SEQUENTIAL'
            record['mkl_import_library_sha256']=hashlib.sha256((args.mkl_root/'lib/mkl_rt.lib').read_bytes()).hexdigest()
            extra=['-DRADIA_HMATRIX_FULL_HOST','-L'+str(args.mkl_root/'lib'),'-lmkl_rt']
        names=['cHACApK_base','cHACApK_cpp_impl','cHACApK_pca_cluster','cHACApK_lib','hmatrix_serial_host'] if full else ['hmatrix_aca_host','cHACApK_lib']
        for name in names:
            obj=out/(name+'.obj');batch=out/(name+'.cmd');objects.append(str(obj))
            source=root/('packages/radia-gpu/native' if name.startswith('hmatrix_') else 'src/ext/HACApK')/(name+'.c')
            include=f'/I"{args.mkl_root/"include"}"' if full else ''
            batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\n'+f'cl /nologo /c /TC /MD /O2 /Gy /DRAD_HACAPK_BATCH_ENTRIES {include} /I"{root/"src/core"}" /I"{root/"src/ext/HACApK"}" "{source}" /Fo:"{obj}"\n',encoding='utf-8')
            compiled=subprocess.run(['cmd','/c',str(batch)],capture_output=True)
            (out/(name+'-compile.log')).write_bytes(compiled.stdout+compiled.stderr);compiled.check_returncode()
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
        pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        eng.setenv('NVCC_CCBIN',str(cc),nargout=0);eng.setenv('PATH',str(cc)+';'+eng.getenv('PATH'),nargout=0)
        config=Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        if not config.is_file():raise RuntimeError('Supported MATLAB CUDA configuration missing')
        eng.mexcuda('-v','-R2018a','-f',str(config),'NVCCFLAGS=-std=c++17 --fmad=false -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80','-I'+str(root/'src/ext/HACApK'),'-I'+str(root/'src/core'),'-outdir',str(out),str(files[0]),*objects,*extra,'-lcublas','LINKFLAGS=$LINKFLAGS /OPT:REF',nargout=0)
        if hashes!={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}:raise RuntimeError('Source changed during build')
        if subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()!=record['source_commit']:raise RuntimeError('Commit changed during build')
        record.update(accepted=True,matlab=eng.version(),binary_sha256=hashlib.sha256((out/'radia_cuda_hmatrix.mexw64').read_bytes()).hexdigest(),cuda_targets=['sm_75','sm_80','compute_80 PTX'],taskmanager_started=False,python_runtime_required=False)
    except BaseException as error:record['error']=f'{type(error).__name__}: {error}';raise
    finally:
        if eng is not None:eng.quit()
        (out/'build.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
