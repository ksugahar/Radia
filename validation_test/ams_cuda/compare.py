"""Interleave two AMS MEX builds on identical saved acceptance inputs."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import matlab.engine
from scipy.io import loadmat
from run import jsonable


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    builds = [args.before.resolve(), args.after.resolve()]
    record = {'schema': 'radia.ams-cuda-comparison.v1', 'accepted': False,
              'order': ['before', 'after', 'after', 'before'] * args.rounds,
              'builds': [], 'cases': []}
    eng = None
    try:
        for directory in builds:
            meta = json.loads((directory / 'build.json').read_text())
            mex = directory / 'radia_cuda_ams.mexw64'
            digest = hashlib.sha256(mex.read_bytes()).hexdigest()
            hashes = {key.replace('\\', '/'): value for key, value in meta['sha256'].items()}
            if digest != hashes['build-ams-cuda/radia_cuda_ams.mexw64']:
                raise RuntimeError('MEX does not match recorded build hash')
            record['builds'].append({'source_commit': meta['source_commit'],
                                     'source_dirty': meta['source_dirty'], 'mex_sha256': digest})
        eng = matlab.engine.start_matlab('-nodesktop -nosplash')
        eng.addpath(str(Path(__file__).resolve().parent), nargout=0)
        for source in args.inputs:
            source = source.resolve()
            snapshot = loadmat(source, simplify_cells=True)['snapshot']
            entry = {'input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                     'dofs': int(snapshot['A']['rows']), 'cycle': int(snapshot['cycle_type']),
                     'runs': []}
            times = [[], []]
            inclusive = [[], []]
            for block, variant in enumerate([0, 1, 1, 0] * args.rounds):
                target = args.out / f'case{len(record["cases"])}-block{block}.mat'
                eng.eval('clear radia_cuda_ams;', nargout=0)
                if block:
                    eng.rmpath(str(builds[previous]), nargout=0)
                eng.benchmark(str(source), str(target), str(builds[variant]), nargout=0)
                previous = variant
                data = jsonable(loadmat(target, simplify_cells=True))
                data.pop('solution'); data.pop('inclusiveSolution')
                samples = data['times']
                times[variant].extend(samples)
                inclusive[variant].append(data['inclusiveSeconds'])
                if (max(data['probeErrors']) > 1e-8 or
                        data['info']['true_relative_residual'] > 1e-9 or
                        not all(data['checks'].values())):
                    raise RuntimeError('Interleaved AMS parity or gateway check failed')
                entry['runs'].append({'variant': ['before', 'after'][variant], **data})
                print(f'case {len(record["cases"])} block {block}: '
                      f'{["before", "after"][variant]} {statistics.median(samples):.6f}s', flush=True)
            eng.eval('clear radia_cuda_ams;', nargout=0)
            eng.rmpath(str(builds[previous]), nargout=0)
            medians = list(map(statistics.median, times))
            entry.update(before_warm_median_seconds=medians[0],
                         after_warm_median_seconds=medians[1],
                         warm_speedup=medians[0] / medians[1],
                         warm_improvement_percent=100 * (1 - medians[1] / medians[0]),
                         warm_timing_improved=medians[1] < medians[0],
                         before_inclusive_median_seconds=statistics.median(inclusive[0]),
                         after_inclusive_median_seconds=statistics.median(inclusive[1]))
            record['cases'].append(entry)
            (args.out / 'result.json').write_text(json.dumps(record, indent=2) + '\n')
        record['accepted'] = True
    except Exception as exc:
        record['error'] = repr(exc)
        raise
    finally:
        if eng is not None:
            eng.quit()
        (args.out / 'result.json').write_text(json.dumps(record, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--rounds', type=int, default=3)
    args = parser.parse_args()
    if args.rounds < 1:
        parser.error('--rounds must be positive')
    main(args)
