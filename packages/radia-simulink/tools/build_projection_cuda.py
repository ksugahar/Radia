"""Build the experimental NGSolve projection mexcuda bridge in an owned job."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,socket,subprocess,sys,tempfile
import importlib.util
import importlib.metadata

def create_build_directory():
    """CI scratch belongs to the runner; local callers own their temp settings."""
    if os.environ.get('GITHUB_ACTIONS') == 'true':
        value = os.environ.get('RUNNER_TEMP')
        if not value:
            raise RuntimeError('CI projection build requires RUNNER_TEMP')
        base = Path(value)
        if not base.is_absolute() or not base.is_dir():
            raise RuntimeError('RUNNER_TEMP must be an existing absolute directory')
    else:
        base = Path(tempfile.gettempdir())
    return Path(tempfile.mkdtemp(prefix='projection-cuda-', dir=base))

def projection_dependencies():
    spec=importlib.util.find_spec('netgen')
    if spec is None or not spec.origin:
        raise RuntimeError('The selected interpreter requires Netgen')
    netgen=Path(spec.origin).resolve().parent
    prefix=Path(sys.base_prefix)
    files=[netgen/'lib'/name for name in ('libngsolve.lib','ngcore.lib','nglib.lib')]
    files.append(prefix/'libs/python312.lib')
    for directory in (netgen/'include',prefix/'include'):
        if not directory.is_dir():raise RuntimeError('Native dependency headers missing: '+str(directory))
        files.extend(p for p in directory.rglob('*') if p.is_file() and p.suffix in ('.h','.hpp','.hxx'))
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    identity=dict(interpreter=sys.executable,python_base_prefix=str(prefix),netgen_root=str(netgen),
        versions={name:importlib.metadata.version(name) for name in ('ngsolve','netgen-mesher')},
        native_dependency_sha256=hashes)
    return prefix,netgen,identity

def main():
    from check_mex_error_boundaries import require_boundaries
    require_boundaries()
    from vs2022_toolchain import discover_vs, compiler_bin
    from matlab_cuda_toolchain import configure_engine, require_build_admission
    parser=argparse.ArgumentParser();parser.add_argument('--require-clean-source',action='store_true');parser.add_argument('--build-admission',type=Path);parser.add_argument('--build-admission-sha256');parser.add_argument('--lifecycle-diagnostic',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parents[3]
    dirty=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip())
    if args.require_clean_source and dirty:raise RuntimeError('Clean source required')
    admission=require_build_admission(root,'projection',args.build_admission,args.build_admission_sha256)
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
    source=root/'packages/radia-simulink/native/radia_cuda_projection.cu'
    out=create_build_directory()
    receipt=dict(accepted=False,a100_used=False,scope='Experimental native mapped-target replay, not full RadFld migration',source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),taskmanager_started=False)
    receipt.pop('a100_used')
    receipt['build_admission']=admission
    receipt['lifecycle_diagnostic']=args.lifecycle_diagnostic
    receipt['fault_injection']=False  # This builder exposes no allocation-fault mode.
    sources=[source,source.with_name('radia_projection_host.cpp'),source.with_name('radia_projection_host.cuh'),source.with_name('hdiv_rt0_host.cuh'),source.with_name('hdiv_element_setup.cuh'),root/'src/core/gpu/hdiv/hdiv_material_staging.hpp',Path(__file__).with_name('matlab_cuda_toolchain.py'),Path(__file__).with_name('vs2022_toolchain.py'),Path(__file__),root/'Build.ps1']
    sources.append(source.with_name('mex_error_boundary.hpp'))
    sources.append(Path(__file__).with_name('check_mex_error_boundaries.py'))
    receipt['source_dirty']=dirty
    receipt['source_sha256']={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    print(str(out),flush=True);eng=None
    try:
        vs=discover_vs();cc=compiler_bin(vs)
        os.environ['NVCC_CCBIN']=str(cc);os.environ['PATH']=str(cc)+os.pathsep+os.environ['PATH']
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash')
        pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        prefix,netgen,dependencies=projection_dependencies()
        receipt['native_dependencies']=dependencies
        config, toolchain_receipt = configure_engine(eng)
        receipt['toolchain']=toolchain_receipt
        (out/'toolchain.json').write_text(json.dumps(toolchain_receipt, indent=2)+'\n', encoding='utf-8')
        obj=out/'radia_projection_host.obj'
        includes=[netgen/'include'/d for d in ('','include','general','linalg','meshing','gprim')]+[prefix/'include',root/'packages/radia-simulink/native',Path(eng.matlabroot())/'extern/include']
        cmd=['cl','/nologo','/c','/O2','/MD','/EHsc','/std:c++20','/DWIN32','/DNOMINMAX','/DMATLAB_MEXCMD_RELEASE=R2018a','/DUSE_MEX_CMD','/DMATLAB_MEX_FILE','/arch:AVX2',f'/Fo{obj}',*[f'/I{p}' for p in includes],str(root/'packages/radia-simulink/native/radia_projection_host.cpp')]
        if args.lifecycle_diagnostic:cmd.insert(1,'/DRADIA_PROJECTION_LIFECYCLE_TEST')
        batch=out/'compile.cmd'
        batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\n'+subprocess.list2cmdline(cmd)+'\nexit /b %errorlevel%\n',encoding='utf-8')
        subprocess.run(['cmd','/c',str(batch)],check=True)
        eng.mexcuda('-v','-f',str(config),'-R2018a','-DWIN32','-DNOMINMAX','NVCCFLAGS=-std=c++20 --fmad=false -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80',*[f'-I{netgen / "include" / d}' for d in ('','include','general','linalg','meshing','gprim')],f'-I{prefix / "include"}',f'-I{root / "packages/radia-simulink/native"}',f'-L{netgen / "lib"}',f'-L{prefix / "libs"}','-llibngsolve','-lnglib','-lngcore','-lpython312','-outdir',str(out),str(source),str(obj),nargout=0)
        if receipt['source_sha256']!={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}:raise RuntimeError('Source changed during build')
        if projection_dependencies()[2]!=dependencies:raise RuntimeError('Native dependencies changed during build')
        binary=out/'radia_cuda_projection.mexw64' 
        receipt.update(accepted=True,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),matlab=eng.version())
    except BaseException as e:receipt['error']=str(e);raise
    finally:
        if eng is not None:eng.quit()
        (out/'build.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()

