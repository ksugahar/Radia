"""Build the MEX-free pybind CSR adapter using the admitted bundled toolchain.

Called by Build.ps1 -CudaMexOnly -CsrPythonBindings. Copies no third-party DLL.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import sysconfig
import os
import shutil



def snapshot_headers(source, destination):
    """Compile exact SDK bytes from a neutral relative include directory."""
    source, destination = Path(source), Path(destination)
    entries = list(source.rglob('*'))
    if any(p.is_symlink() or getattr(p, 'is_junction', lambda: False)() for p in entries):
        raise RuntimeError('Linked dependency headers are not admitted')
    hashes = {p.relative_to(source).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in entries if p.is_file()}
    shutil.copytree(source, destination)
    if any(hashlib.sha256((destination/name).read_bytes()).hexdigest() != digest
           for name, digest in hashes.items()):
        raise RuntimeError('Dependency header snapshot mismatch')
    return hashes

def build_python_binding(engine, package, out, config, targets, toolchain):
    import ngsolve
    import netgen
    from vs2022_toolchain import discover_vs, compiler_bin
    root=Path(engine.matlabroot()).resolve()
    native=package.parents[1]/'src/core/gpu'
    component='csr'
    stem='csr_pcg_python'
    module='_csr_cuda'
    include=Path(netgen.__file__).parent/'include'
    common=include/'pybind11/detail/common.h'
    if not (include/'pybind11/pybind11.h').is_file() or not common.is_file():
        raise RuntimeError('NGSolve/Netgen pybind11 development headers required')
    common_bytes=common.read_bytes()
    # Match the installed NGSolve build headers, not an independently upgraded pybind.
    import re
    version_fields=re.findall(r'^#define PYBIND11_VERSION_(MAJOR|MINOR|PATCH)\s+(\d+)',common_bytes.decode(),re.M)
    pybind_version=dict(version_fields)
    if set(pybind_version)!={'MAJOR','MINOR','PATCH'}:
        raise RuntimeError('Unrecognized NGSolve pybind11 header version')
    engine.mexcuda('-v','-c','-f',str(config),targets,'-outdir',str(out),
                   str(native/(stem+'_cuda.cu')),nargout=0)
    obj=out/(stem+'_cuda.obj')
    if not obj.is_file():
        raise RuntimeError('mexcuda did not produce the standalone domain object')
    original_headers=out/'dependency-headers-original/netgen'
    header_root=out/'dependency-headers/netgen'
    header_hashes=snapshot_headers(include,original_headers)
    from native_source_paths import stage_logical_sources
    logical_headers=stage_logical_sources([p for p in original_headers.rglob('*') if p.is_file()],original_headers,header_root,'ngsolve/include')
    if any(v['original_sha256']!=header_hashes[k] for k,v in logical_headers.items()):raise RuntimeError('Header snapshot changed')
    (out/'logical-header-paths.json').write_text(json.dumps(logical_headers,sort_keys=True,indent=2))
    common=original_headers/'pybind11/detail/common.h'
    if common.read_bytes()!=common_bytes:raise RuntimeError('NGSolve pybind headers changed during staging')
    (out/'dependency-header-snapshot.json').write_text(json.dumps(header_hashes,sort_keys=True,indent=2))
    relative=lambda p: os.path.relpath(p,out)
    vs=discover_vs();compiler=compiler_bin(vs)
    binary=out/(module+sysconfig.get_config_var('EXT_SUFFIX'))
    host_obj=out/(stem+'.obj')
    python_lib=Path(sys.base_prefix)/'libs'/f'python{sys.version_info.major}{sys.version_info.minor}.lib'
    command=[str(compiler/'cl.exe'),'/nologo','/LD','/EHsc','/MD','/O2','/fp:precise','/std:c++20',
             '/utf-8',f'/I{relative(header_root)}',f'/I{sysconfig.get_path("include")}',
             f'/Fo{host_obj}',relative(native/(stem+'.cpp')),str(obj),'/link',f'/OUT:{binary}',
             f'/IMPLIB:{out / (stem+".lib")}',str(python_lib),
             str(root/'sys/cuda/win64/cuda/lib/x64/cudart.lib')]
    if component=='csr':
        command.insert(command.index('/link'),'/DRADIA_CSR_PYTHON_MODULE=_csr_cuda')
    batch=out/('compile-'+component+'-python.cmd')
    batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+
                     '" >nul\nif errorlevel 1 exit /b 1\n'+subprocess.list2cmdline(command)+
                     '\nexit /b %errorlevel%\n',encoding='utf-8')
    with (out/'python-native-build.log').open('w') as log:
        subprocess.run(['cmd','/c',str(batch)],stdout=log,stderr=subprocess.STDOUT,check=True,cwd=out)
    dependencies=subprocess.check_output([str(compiler/'dumpbin.exe'),'/DEPENDENTS',str(binary)],text=True)
    (out/'python-native-dependencies.txt').write_text(dependencies,encoding='utf-8')
    result=dict(schema='radia.'+component+'-python-build.v1',binary=binary.name,
                binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                python_version=sys.version,extension_suffix=sysconfig.get_config_var('EXT_SUFFIX'),
                pybind11_version=pybind_version,pybind_header_sha256=hashlib.sha256(common.read_bytes()).hexdigest(),ngsolve_version=ngsolve.__version__,
                netgen_module=str(netgen.__file__),toolchain=toolchain,
                cuda_object_sha256=hashlib.sha256(obj.read_bytes()).hexdigest(),
                dependencies=dependencies,matlab_engine_user_dependency=False,
                logical_header_paths_sha256=hashlib.sha256((out/'logical-header-paths.json').read_bytes()).hexdigest(),
                dependency_header_snapshot_sha256=hashlib.sha256((out/'dependency-header-snapshot.json').read_bytes()).hexdigest(),
                runtime_import_checked=False,redistribution_accepted=False,fault_injection=False,device_fault_injection=False)
    (out/'python-build.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result
