"""Staged standalone AMG pybind module; invoked through Build.ps1 only."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import sysconfig

def stage_binding_sources(package, out, include):
    """Own SDK copies and preserve source bodies with neutral diagnostic names."""
    from build_csr_python import snapshot_headers
    from native_source_paths import stage_logical_sources
    package, out, include = Path(package), Path(out), Path(include)
    original = out / 'dependency-headers-original/netgen'
    headers = out / 'dependency-headers/netgen'
    hashes = snapshot_headers(include, original)
    logical_headers = stage_logical_sources(
        [p for p in original.rglob('*') if p.is_file()], original, headers,
        'ngsolve/include')
    if any(value['original_sha256'] != hashes[name]
           for name, value in logical_headers.items()):
        raise RuntimeError('Resident AMG dependency header snapshot changed')
    native = package.parents[1] / 'src/core/gpu'
    logical_native = out / 'logical-native'
    native_paths = stage_logical_sources(
        [native/'resident_amg_python.cpp', native/'resident_amg_host.hpp'],
        native, logical_native, 'radia/native')
    records = {'dependency-header-snapshot': hashes,
               'logical-header-paths': logical_headers,
               'logical-native-paths': native_paths}
    provenance = {}
    for name, record in records.items():
        path = out / (name + '.json')
        path.write_text(json.dumps(record, sort_keys=True, indent=2))
        provenance[name.replace('-', '_') + '_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    return headers, logical_native/'resident_amg_python.cpp', provenance

def build(package,out,admission):
    import matlab.engine
    import netgen
    from matlab_cuda_toolchain import configure_engine
    from vs2022_toolchain import discover_vs,compiler_bin
    out.mkdir(parents=True,exist_ok=True)
    native=package.parents[1]/'src/core/gpu'
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
        headers, host_source, header_provenance = stage_binding_sources(package, out, include)
        command=[str(compiler/'cl.exe'),'/nologo','/LD','/EHsc','/MD','/O2','/fp:precise','/std:c++20','/utf-8','/DRADIA_AMG_PYTHON_MODULE=_resident_amg',f'/I{headers.relative_to(out)}',f'/I{sysconfig.get_path("include")}',f'/Fo{out/"resident_amg_python.obj"}',str(host_source.relative_to(out)),str(obj),'/link',f'/OUT:{binary}',f'/IMPLIB:{out/"resident_amg.lib"}',str(Path(sys.base_prefix)/'libs'/f'python{sys.version_info.major}{sys.version_info.minor}.lib')]+[str(root/f'sys/cuda/win64/cuda/lib/x64/{name}.lib') for name in ['cudart','cublas']]
        batch=out/'compile-resident-amg.cmd';batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\nif errorlevel 1 exit /b 1\n'+subprocess.list2cmdline(command)+'\nexit /b %errorlevel%\n')
        with (out/'resident-amg-host-build.log').open('w') as log:subprocess.run(['cmd','/c',str(batch)],stdout=log,stderr=subprocess.STDOUT,check=True,cwd=out)
        sources=[native/f'resident_amg_{suffix}' for suffix in ['cuda.cu','python.cpp','host.hpp']]
        sources += [Path(__file__), Path(__file__).with_name('build_cuda_mex.py'),
                    Path(__file__).with_name('matlab_cuda_toolchain.py'),
                    Path(__file__).with_name('vs2022_toolchain.py'),
                    Path(__file__).with_name('build_csr_python.py'),
                    Path(__file__).with_name('native_source_paths.py'), repo/'Build.ps1']
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
        record.update(header_provenance)
        (out/'resident-amg-build.json').write_text(json.dumps(record,indent=2)+'\n')
    finally:engine.quit()
