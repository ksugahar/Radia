"""Validate installed IH bytes; Python is independent CPU validation only."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'packages/radia-gpu'))
from verify_wheel import verify, verify_installed_identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheel', type=Path, required=True)
    parser.add_argument('--installed-python', type=Path, required=True)
    parser.add_argument('--toolkit-tools', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    record = {'accepted': False, 'scope': 'Installed IH numerical and lifecycle subgate',
              'full_wheel_accepted': False, 'experimental_validation': False}
    engine = None
    try:
        wheel = verify(args.wheel)
        installed = json.loads(subprocess.check_output(
            [str(args.installed_python), '-I', '-c',
             'import json,radia_gpu; print(json.dumps(radia_gpu.verify()))'], text=True))
        verify_installed_identity(installed, wheel)
        payload = Path(installed['matlab_path']).resolve()
        binary = payload / 'radia_cuda_ih.mexw64'
        receipt = json.loads((payload.parent / 'ih-build.json').read_text(encoding='utf-8'))
        digest = hashlib.sha256(binary.read_bytes()).hexdigest()
        if digest != receipt['binary_sha256'] or receipt['source_commit'] != wheel['source_commit']:
            raise ValueError('Installed IH receipt mismatch')
        if not (args.toolkit_tools / 'model_check/model_check.p').is_file():
            raise ValueError('Official Simulink Toolkit tools are required')
        import matlab.engine
        engine = matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
        engine.addpath(str(ROOT / 'validation_test/ih_cuda'), nargout=0)
        result = engine.checkTemperatureBH(str(payload), str(args.out), str(payload), nargout=1)
        if not result['accepted'] or Path(engine.which('radia_cuda_ih')).resolve() != binary:
            raise ValueError('Installed IH numerical/lifecycle acceptance failed')
        engine.addpath(engine.genpath(str(args.toolkit_tools.resolve())), nargout=0)
        models = list(args.out.glob('*.slx'))
        if len(models) != 1:
            raise ValueError('Expected one IH model')
        model = models[0].resolve()
        engine.load_system(str(model), nargout=0)
        report = engine.model_check(str(model), nargout=1)
        if 'status: healthy' not in str(report):
            raise ValueError('Official Toolkit model check failed')
        engine.save_system(model.stem, str(model), nargout=0)
        engine.close_system(model.stem, 0, nargout=0)
        engine.load_system(str(model), nargout=0)
        if 'status: healthy' not in str(engine.model_check(str(model), nargout=1)):
            raise ValueError('Reopened model check failed')
        engine.close_system(model.stem, 0, nargout=0)
        subprocess.run([sys.executable, str(ROOT / 'validation_test/ih_cuda/check_temperature_bh_cpu.py'),
                        str(args.out)], check=True)
        cpu = json.loads((args.out / 'cpu-result.json').read_text(encoding='utf-8'))
        if max(cpu['relative_errors'].values()) > 1e-8 or cpu['true_relative_residual'] > 1e-6:
            raise ValueError('Release CPU parity gate failed')
        record.update(accepted=True, source_commit=wheel['source_commit'],
                      wheel_sha256=wheel['wheel_sha256'], manifest_sha256=wheel['manifest_sha256'],
                      binary_sha256=digest, cpu=cpu, toolkit_accepted=True,
                      matlab_version=engine.version(), loaded_mex=str(binary))
    except Exception as error:
        record['error'] = str(error)
        raise
    finally:
        if engine is not None:
            engine.quit()
        (args.out / 'installed-result.json').write_text(
            json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
