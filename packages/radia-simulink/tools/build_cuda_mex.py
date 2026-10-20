"""Build the isolated CUDA-only NGSolve assembled-operator bridge on an explicitly admitted build host."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import socket
import subprocess
import sys
import matlab.engine

def main():
    from matlab_cuda_toolchain import configure_engine,require_build_admission
    from check_mex_error_boundaries import require_boundaries
    require_boundaries()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--build-admission',type=Path);parser.add_argument('--build-admission-sha256')
    args = parser.parse_args()
    package = Path(__file__).resolve().parents[1]
    root = package.parents[1]
    admission=require_build_admission(root,'csr-ic',args.build_admission,args.build_admission_sha256)
    sources = [package / 'native/radia_cuda_csr_pcg.cu', package / 'native/radia_cuda_ic.cu']
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    eng = matlab.engine.start_matlab('-nodesktop -nosplash')
    try:
        matlab_root = Path(eng.matlabroot())
        config,toolchain=configure_engine(eng)
        targets = 'NVCCFLAGS=-gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80'
        for source in sources:
            libraries = ['-lcusparse','-lcublas'] if source.stem == 'radia_cuda_ic' else []
            eng.mexcuda('-v', '-f', str(config), targets, *libraries, '-outdir', str(out), str(source), nargout=0)
        source = sources[0]
        binary = out / 'radia_cuda_csr_pcg.mexw64'
        assert binary.is_file()
        device = json.loads(eng.eval("jsonencode(struct('name',gpuDevice().Name,'compute_capability',gpuDevice().ComputeCapability,'total_memory',gpuDevice().TotalMemory,'available_memory',gpuDevice().AvailableMemory,'toolkit',gpuDevice().ToolkitVersion,'driver',gpuDevice().GraphicsDriverVersion))"))
        record = dict(accepted=True,build_admission=admission,toolchain=toolchain,schema='radia.cuda-mex-build.v1', hostname=socket.gethostname(),
                      source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                      binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                      source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
                      source_dirty=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip()),
                      matlab=eng.version(), release=eng.eval("version('-release')"),
                      compiler_configuration=str(config), gpu=device, python=sys.executable,
                      python_version=platform.python_version(), operation_route='official MATLAB Engine for Python; task-owned session',
                      no_cpu_solver_fallback=True, cuda_targets=['sm_75','sm_80','compute_80 PTX'],
                      source_sha256_by_file={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources+[package/'native/mex_error_boundary.hpp',Path(__file__).with_name('check_mex_error_boundaries.py'),Path(__file__).with_name('matlab_cuda_toolchain.py'),Path(__file__).with_name('vs2022_toolchain.py'),Path(__file__),root/'Build.ps1']},
                      binary_sha256_by_file={p.with_suffix('.mexw64').name:hashlib.sha256((out/p.with_suffix('.mexw64').name).read_bytes()).hexdigest() for p in sources})
        (out/'build.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        # Keep operational identity in the local receipt, not public CI output.
        print(json.dumps({key:value for key,value in record.items()
                          if key not in ('hostname','python','compiler_configuration')},indent=2))
    finally:
        eng.quit()

if __name__ == '__main__': main()
