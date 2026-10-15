"""Single-GPU test-only mismatch regression; never installs a runtime."""
import argparse
import hashlib
import json
from pathlib import Path
import matlab.engine


def check(engine, out):
    root = Path(__file__).resolve().parents[2]
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = root / 'packages/radia-gpu/native/radia_cuda_ic.cu'
    config = Path(engine.matlabroot()) / 'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
    if not config.is_file():
        raise RuntimeError('Supported mexcuda MSVC configuration missing')
    engine.mexcuda('-f', str(config), '-DRADIA_IC_LIFECYCLE_TEST',
                  'NVCCFLAGS=-gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80',
                  '-lcusparse', '-lcublas', '-outdir', str(out), str(source), nargout=0)
    engine.addpath(str(Path(__file__).parent), nargout=0)
    engine.addpath(str(out), '-begin', nargout=0)
    result = dict(engine.validate_ic_lifecycle(True, nargout=1))
    result['binary_sha256'] = hashlib.sha256((out / 'radia_cuda_ic.mexw64').read_bytes()).hexdigest()
    result['source_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
    result['limitation'] = 'Simulates selected-device mismatch; physical multi-GPU switching remains unverified'
    (out / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    engine = matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
    try:
        print(json.dumps(check(engine, args.out), indent=2))
    finally:
        engine.quit()
