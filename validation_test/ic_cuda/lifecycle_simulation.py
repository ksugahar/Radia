"""Single-GPU test-only mismatch regression; never installs a runtime."""
import argparse
import hashlib
import json
from pathlib import Path
import matlab.engine


SOURCE = Path(__file__).resolve().parents[2] / 'packages/radia-simulink/native/radia_cuda_ic.cu'


def build(engine, out):
    """Compile the test-only IC gateway on the build host; creates no CUDA context."""
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    config = Path(engine.matlabroot()) / 'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
    if not config.is_file():
        raise RuntimeError('Supported mexcuda MSVC configuration missing')
    engine.mexcuda('-f', str(config), '-DRADIA_IC_LIFECYCLE_TEST',
                  'NVCCFLAGS=-gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80',
                  '-lcusparse', '-lcublas', '-outdir', str(out), str(SOURCE), nargout=0)
    receipt = {'binary_sha256': hashlib.sha256((out / 'radia_cuda_ic.mexw64').read_bytes()).hexdigest(),
               'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest()}
    (out / 'build.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


def check(engine, binary_dir, out):
    """Run the lifecycle regression on a prebuilt gateway of this exact source; never compiles."""
    binary_dir, out = Path(binary_dir).resolve(), Path(out).resolve()
    receipt = json.loads((binary_dir / 'build.json').read_text())
    binary = binary_dir / 'radia_cuda_ic.mexw64'
    if hashlib.sha256(SOURCE.read_bytes()).hexdigest() != receipt['source_sha256']:
        raise RuntimeError('Prebuilt IC lifecycle gateway was compiled from a different source')
    if hashlib.sha256(binary.read_bytes()).hexdigest() != receipt['binary_sha256']:
        raise RuntimeError('Prebuilt IC lifecycle gateway does not match its build receipt')
    out.mkdir(parents=True, exist_ok=False)
    engine.addpath(str(Path(__file__).parent), nargout=0)
    engine.addpath(str(binary_dir), '-begin', nargout=0)
    result = dict(engine.validate_ic_lifecycle(True, nargout=1))
    result.update(receipt)
    result['limitation'] = 'Simulates selected-device mismatch; physical multi-GPU switching remains unverified'
    (out / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--build-only', action='store_true', help='compile the gateway into --out (build host)')
    mode.add_argument('--binary-dir', type=Path, help='run on a prebuilt gateway (GPU host)')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    engine = matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
    try:
        result = build(engine, args.out) if args.build_only else check(engine, args.binary_dir, args.out)
        print(json.dumps(result, indent=2))
    finally:
        engine.quit()
