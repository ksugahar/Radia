"""Isolated resident HACApK algebra CUDA build, not a physical solver release."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,socket,subprocess,tempfile
import re


def project_header_inputs(files, root):
    """Hash quoted project includes transitively; external headers belong to F4."""
    result=list(dict.fromkeys(files));seen=set(result)
    search=[root/'packages/radia-simulink/native',root/'src/core/gpu/hdiv',root/'src/core',root/'src/ext/HACApK']
    for path in result:
        if path.suffix not in ('.h','.hpp','.cuh','.cu','.c','.cpp'):
            continue
        for name in re.findall(r'^\s*#\s*include\s*"([^"]+)"',path.read_text(encoding='utf-8'),re.M):
            for base in [path.parent,*search]:
                candidate=(base/name).resolve()
                if candidate.is_file() and candidate.is_relative_to(root.resolve()):
                    if candidate not in seen:seen.add(candidate);result.append(candidate)
                    break
    return result


def main():
    from check_mex_error_boundaries import require_boundaries
    require_boundaries()
    from vs2022_toolchain import discover_vs, compiler_bin
    from matlab_cuda_toolchain import configure_engine, require_build_admission
    parser=argparse.ArgumentParser();parser.add_argument('--require-clean-source',action='store_true');parser.add_argument('--mkl-root',type=Path);parser.add_argument('--build-admission',type=Path);parser.add_argument('--build-admission-sha256');parser.add_argument('--entry-batch-size',type=int,choices=[256,1024],default=256);parser.add_argument('--fault-injection',action='store_true');args=parser.parse_args()
    if args.fault_injection and not args.mkl_root:parser.error('Fault injection requires the full serial C build')
    root=Path(__file__).resolve().parents[3]
    admission=require_build_admission(root,'hmatrix',args.build_admission,args.build_admission_sha256)
    dirty=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip())
    if args.require_clean_source and dirty:raise RuntimeError('Clean source required')
    files=[root/'packages/radia-simulink/native/radia_cuda_hmatrix.cu',root/'src/core/gpu/hdiv/hmatrix_resident.cuh',root/'packages/radia-simulink/native/hmatrix_hacapk_pack.cuh',root/'packages/radia-simulink/native/hmatrix_entry_oracle.cuh',root/'packages/radia-simulink/native/hmatrix_aca_callback.cuh',root/'packages/radia-simulink/native/hmatrix_leaf_bounds.h',root/'src/ext/HACApK/cHACApK_calc_entry_ij.h',root/'src/ext/HACApK/cHACApK_base.h',root/'src/core/gpu/hdiv/analytic_matrix_build.cuh',root/'src/core/gpu/hdiv/radia_field_exact.cuh',root/'src/core/rad_charge_quadrature.h',root/'src/core/rad_hdiv_vim.h',root/'src/core/rad_field_cloud_kernel.h',root/'src/core/rad_hdiv_field_math.inc',root/'src/core/rad_hdiv_field_kernel_api.inc',root/'src/core/rad_field_source_kernel.inc',Path(__file__).with_name('matlab_cuda_toolchain.py'),Path(__file__).with_name('vs2022_toolchain.py'),Path(__file__),root/'Build.ps1']
    files += [root/'src/ext/HACApK/cHACApK_base.c',root/'src/ext/HACApK/cHACApK_aca.inc',root/'packages/radia-simulink/native/hmatrix_aca_host.c',root/'src/ext/HACApK/cHACApK_lib.c',root/'src/core/rad_hacapk_parallel.h'] + sorted((root/'src/ext/HACApK').glob('*.h'))
    files += [root/'packages/radia-simulink/native/hmatrix_curved_entry_cache.cuh',root/'packages/radia-simulink/native/hmatrix_curved_entry_oracle.cuh',root/'src/core/gpu/hdiv/hdiv_assembly_profile.cuh',root/'src/core/gpu/hdiv/hdiv_curved_moments.cuh',root/'src/core/gpu/hdiv/hmatrix_material_solver.cuh',root/'src/core/gpu/hdiv/hdiv_curved_gram.cuh']
    if args.mkl_root:
        files += [root/'packages/radia-simulink/native/hmatrix_full_assembly.cuh',root/'packages/radia-simulink/native/hmatrix_serial_host.c',root/'packages/radia-simulink/native/hmatrix_c_memory.c',root/'packages/radia-simulink/native/hmatrix_c_memory.h',root/'src/ext/HACApK/cHACApK_cpp_impl.c',root/'src/ext/HACApK/cHACApK_pca_cluster.c']
    files += [root/'src/core/gpu/hdiv/hmatrix_host_fault.hpp',root/'src/core/gpu/hdiv/hmatrix_leaf_coverage.hpp']
    files.append(Path(__file__).with_name('check_mex_error_boundaries.py'))
    if args.fault_injection:
        files += [root/'validation_test/hmatrix_cuda/check_c_memory_transaction.c',root/'validation_test/hmatrix_cuda/check_host_fault_allocator.cpp']
    files += [Path(__file__).with_name('native_source_paths.py')]
    files=project_header_inputs(files,root)
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    base=root/'build-hmatrix-cuda';base.mkdir(exist_ok=True);out=Path(tempfile.mkdtemp(prefix='job-',dir=base));print(out,flush=True)
    record=dict(schema='radia.hmatrix-cuda-build.v1',accepted=False,source_dirty=dirty,source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),source_sha256=hashes,a100_used=False,scope='Resident numeric HACApK leaf application only; assembly/physical solver not accepted')
    record['build_admission']=admission
    record['entry_batch_size']=args.entry_batch_size
    record['fault_injection']=args.fault_injection
    record.pop('a100_used')
    eng=None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        cc=compiler_bin(discover_vs());os.environ['NVCC_CCBIN']=str(cc)
        objects=[]
        full=bool(args.mkl_root);extra=[]
        if full:
            if not (args.mkl_root/'include/mkl.h').is_file() or not (args.mkl_root/'lib/mkl_rt.lib').is_file():raise RuntimeError('Explicit supported mkl-devel root required')
            os.environ['MKL_NUM_THREADS']='1';os.environ['MKL_THREADING_LAYER']='SEQUENTIAL'
            record['mkl_import_library_sha256']=hashlib.sha256((args.mkl_root/'lib/mkl_rt.lib').read_bytes()).hexdigest()
            extra=['-DRADIA_HMATRIX_FULL_HOST','-DRAD_HACAPK_DENSE_LEAF_SINK','-DRAD_HACAPK_FRAME_ADMISSION','-L'+str(args.mkl_root/'lib'),'-lmkl_rt']
            if args.fault_injection:extra.append('-DRADIA_HMATRIX_FAULT_INJECTION')
        from native_source_paths import stage_logical_sources
        logical_root=out/'logical-source'
        logical_sources=stage_logical_sources(files,root,logical_root,'radia')
        (out/'logical-source-paths.json').write_text(json.dumps(logical_sources,sort_keys=True,indent=2))
        record['logical_source_paths_sha256']=hashlib.sha256((out/'logical-source-paths.json').read_bytes()).hexdigest()
        names=['cHACApK_base','cHACApK_cpp_impl','cHACApK_pca_cluster','cHACApK_lib','hmatrix_serial_host','hmatrix_c_memory'] if full else ['hmatrix_aca_host','cHACApK_lib']
        for name in names:
            obj=out/(name+'.obj');batch=out/(name+'.cmd');objects.append(str(obj))
            source=root/('packages/radia-simulink/native' if name.startswith('hmatrix_') else 'src/ext/HACApK')/(name+'.c')
            include=(f'/I"{args.mkl_root/"include"}" /DRAD_HACAPK_DENSE_LEAF_SINK /DRAD_HACAPK_FRAME_ADMISSION') if full else ''
            if full and name.startswith('cHACApK_'):
                include += f' /DRAD_HACAPK_C_MEMORY /I"{logical_root/"packages/radia-simulink/native"}" /FI"hmatrix_c_memory.h"'
            batch.write_text('@echo off\ncall "'+str(discover_vs()/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\n'+f'cl /nologo /c /TC /MD /O2 /Gy /DRADIA_BLAS_MKL=1 /DRAD_HACAPK_BATCH_ENTRIES /DRAD_HACAPK_BATCH_SIZE={args.entry_batch_size} {include} /I"{logical_root/"src/core"}" /I"{logical_root/"src/ext/HACApK"}" "{logical_root/source.relative_to(root)}" /Fo:"{obj}"\n',encoding='utf-8')
            compiled=subprocess.run(['cmd','/c',str(batch)],capture_output=True,cwd=root)
            (out/(name+'-compile.log')).write_bytes(compiled.stdout+compiled.stderr);compiled.check_returncode()
        if args.fault_injection:
            record['host_fault_contracts']={}
            for name,language in [('check_c_memory_transaction','c'),('check_host_fault_allocator','cpp')]:
                source=root/'validation_test/hmatrix_cuda'/(name+'.'+language)
                exe=out/(name+'.exe');batch=out/(name+'.cmd')
                linked=[str(out/'hmatrix_c_memory.obj')] if language=='c' else []
                batch.write_text('@echo off\ncall "'+str(discover_vs()/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\n'+
                    f'cl /nologo /MD /Od /EHsc /std:c++17 /I"{root/"packages/radia-simulink/native"}" /I"{root/"src/ext/HACApK"}" /I"{root/"src/core"}" /I"{args.mkl_root/"include"}" "{source}" '+
                    ' '.join('"'+p+'"' for p in linked)+f' /Fe:"{exe}" /Fo:"{out/(name+".obj")}"\n',encoding='utf-8')
                built=subprocess.run(['cmd','/c',str(batch)],capture_output=True,cwd=root)
                (out/(name+'-compile.log')).write_bytes(built.stdout+built.stderr);built.check_returncode()
                tested=subprocess.run([str(exe)],capture_output=True,timeout=30)
                (out/(name+'-run.log')).write_bytes(tested.stdout+tested.stderr);tested.check_returncode()
                record['host_fault_contracts'][name]={'accepted':True,'sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'exit_code':tested.returncode}
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
        pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        eng.setenv('NVCC_CCBIN',str(cc),nargout=0);eng.setenv('PATH',str(cc)+';'+eng.getenv('PATH'),nargout=0)
        config, toolchain_receipt = configure_engine(eng)
        record['toolchain']=toolchain_receipt
        (out/'toolchain.json').write_text(json.dumps(toolchain_receipt, indent=2)+'\n', encoding='utf-8')
        if not config.is_file():raise RuntimeError('Supported MATLAB CUDA configuration missing')
        eng.mexcuda('-v','-R2018a','-f',str(config),'NVCCFLAGS=-std=c++17 --fmad=false --resource-usage -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80','-DRAD_HACAPK_BATCH_SIZE='+str(args.entry_batch_size),'-I'+str(root/'src/ext/HACApK'),'-I'+str(root/'src/core'),'-outdir',str(out),str(files[0]),*objects,*extra,'-lcublas','-lcusolver','LINKFLAGS=$LINKFLAGS /OPT:REF',nargout=0)
        record.update(accepted=True,matlab=eng.version(),binary_sha256=hashlib.sha256((out/'radia_cuda_hmatrix.mexw64').read_bytes()).hexdigest(),cuda_targets=['sm_75','sm_80','compute_80 PTX'],taskmanager_started=False,python_runtime_required=False)
    except BaseException as error:record['error']=f'{type(error).__name__}: {error}';raise
    finally:
        if eng is not None:eng.quit()
        (out/'build.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
