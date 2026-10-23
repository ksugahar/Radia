"""Start and close only this probe's MATLAB Engine; install nothing."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matlab-root', required=True)
    parser.add_argument('--gpu', required=True)
    parser.add_argument('--probe-dir', required=True)
    parser.add_argument('--result', type=Path, required=True)
    args = parser.parse_args()
    if args.result.exists():
        raise FileExistsError(args.result)
    report = {'schema':'radia.runner-engine-readiness.v1', 'passed':False,
              'shutdown':'not_started', 'python':sys.executable, 'python_version':sys.version}
    engine = None
    try:
        import matlab.engine
        report['engine_module'] = matlab.engine.__file__
        report['engine_version'] = importlib.metadata.version('matlabengine')
        # Always a NEW owned process. Never attach to a person's MATLAB.
        engine = matlab.engine.start_matlab('-nodesktop -nosplash')
        engine.addpath(str(Path(args.probe_dir).resolve()), nargout=0)
        payload = engine.probe_matlab_readiness(args.matlab_root,args.gpu,nargout=1)
        report['matlab'] = json.loads(payload)
        if not report['matlab']['passed']:
            raise RuntimeError('Engine MATLAB/license/GPU checks did not pass')
        report['passed'] = True
    except Exception as error:
        report['error'] = f'{type(error).__name__}: {error}'
    finally:
        if engine is not None:
            try:
                engine.quit()
                report['shutdown'] = 'complete'
            except Exception as error:
                report.update(passed=False, shutdown_error=f'{type(error).__name__}: {error}')
        with args.result.open('x',encoding='utf-8') as stream:
            json.dump(report,stream,indent=2)
        print(json.dumps(report))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
