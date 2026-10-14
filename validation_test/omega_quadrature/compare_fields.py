"""Compare saved quadrature re-solves, not certify convergence from two points."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def validate_identity(payload):
    def sha(value):
        return (isinstance(value, str) and len(value) == 64
                and all(c in '0123456789abcdefABCDEF' for c in value))

    if payload.get('schema') != 'radia.validation.omega-quadrature.v1':
        raise ValueError('missing or unsupported evidence schema')
    if not sha(payload.get('mesh_sha256')):
        raise ValueError('missing mesh hash')
    identity = payload.get('implementation', {})
    for key in ('native', 'python_files'):
        hashes = identity.get(key)
        if (not isinstance(hashes, dict) or not hashes
                or not all(isinstance(name, str) and name and sha(value) for name, value in hashes.items())):
            raise ValueError(f'missing or invalid {key} hashes')
    for key in ('solver_sha256', 'factory_sha256', 'runner_sha256', 'diagnostics_sha256'):
        if not sha(identity.get(key)):
            raise ValueError(f'missing {key}')
    for key in ('version', 'ngsolve', 'module', 'solver_file'):
        if not isinstance(identity.get(key), str) or not identity[key].strip():
            raise ValueError(f'missing {key}')
    if payload.get('runtime', {}).get('mode') not in ('wheel', 'research'):
        raise ValueError('missing runtime mode')
    case = payload.get('case', {})
    if case.get('case') != 'esrf6' or not sha(case.get('coil_source_sha256')):
        raise ValueError('missing ESRF source identity')
    for key in ('mu_r', 'radius'):
        value = case.get(key)
        if type(value) not in (float, int) or not np.isfinite(value) or value <= 0:
            raise ValueError(f'invalid physical control: {key}')
    for key in ('kelvin_center', 'physical_center'):
        vector = np.asarray(case.get(key), dtype=float)
        if vector.shape != (3,) or not np.isfinite(vector).all():
            raise ValueError(f'invalid physical control: {key}')


def compare_fields(low, high):
    for payload in (low, high):
        validate_identity(payload)
        for key in ('completed', 'source_unchanged', 'mesh_unchanged'):
            if payload.get(key) is not True:
                raise ValueError(f'unverified input: {key}')
        if len(payload['rows']) != 1:
            raise ValueError('expected one solve per input')
    for key in ('implementation', 'mesh_sha256', 'case', 'runtime'):
        if low[key] != high[key]:
            raise ValueError(f'incompatible {key}')
    for key in ('source_order', 'threads', 'evaluation_order', 'algebraic_only'):
        if low['controls'][key] != high['controls'][key]:
            raise ValueError(f'incompatible control: {key}')
    a, b = low['rows'][0], high['rows'][0]
    if a['order'] != b['order'] or a['ndof'] != b['ndof'] or a['bonus'] >= b['bonus']:
        raise ValueError('expected increasing assembly bonus at fixed FE order and DOFs')
    observations = [row['field_observations'] for row in (a, b)]
    for key in ('centres_m', 'samples_m', 'observable'):
        if observations[0][key] != observations[1][key]:
            raise ValueError(f'incompatible observation {key}')
    values, raw_samples = [], []
    for observation in observations:
        centres = np.asarray(observation['centres_m'], dtype=float)
        points = np.asarray(observation['samples_m'], dtype=float)
        samples = np.asarray(observation['B_samples_T'], dtype=float)
        average = np.asarray(observation['B_average_T'], dtype=float)
        if (centres.ndim != 2 or centres.shape[1] != 3 or len(centres) == 0
                or points.shape != (8*len(centres), 3) or samples.shape != points.shape
                or average.shape != centres.shape):
            raise ValueError('invalid field sample shape')
        if not all(np.isfinite(array).all() for array in (centres, points, samples, average)):
            raise ValueError('nonfinite field sample')
        computed = samples.reshape(-1, 8, 3).mean(axis=1)
        if not np.array_equal(computed, average):
            raise ValueError('saved average does not match samples')
        values.append(computed)
        raw_samples.append(samples)
    delta = values[1] - values[0]
    reference = float(np.linalg.norm(values[1]))
    difference = float(np.linalg.norm(delta))
    raw_delta = raw_samples[1] - raw_samples[0]
    raw_difference = float(np.linalg.norm(raw_delta))
    if not np.isfinite([reference, difference, raw_difference]).all():
        raise ValueError('field norm overflow')
    return {
        'schema': 'radia.validation.omega-quadrature-field-delta.v1',
        'acceptance': 'HOLD: two assembly rules are not a convergence certificate',
        'order': a['order'], 'bonuses': [a['bonus'], b['bonus']],
        'centre_count': len(delta), 'reference_rms_T': reference / np.sqrt(len(delta)),
        'difference_rms_T': difference / np.sqrt(len(delta)),
        'difference_relative_rms': difference / reference if reference else None,
        'difference_max_vector_T': float(np.max(np.linalg.norm(delta, axis=1))),
        'sample_difference_rms_T': raw_difference / np.sqrt(len(raw_delta)),
        'sample_difference_max_vector_T': float(np.max(np.linalg.norm(raw_delta, axis=1))),
        'mesh_sha256': low['mesh_sha256'],
        'native': low['implementation']['native'],
        'scope': 'linear fixed-source/fixed-order assembly-rule sensitivity only',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('low', type=Path)
    parser.add_argument('high', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    for path in (args.low, args.high):
        if args.output.resolve() == path.resolve() or (args.output.exists() and args.output.samefile(path)):
            parser.error('output must not overwrite evidence')
    raw = [path.read_bytes() for path in (args.low, args.high)]
    result = compare_fields(*(json.loads(data) for data in raw))
    result['input_sha256'] = [hashlib.sha256(data).hexdigest() for data in raw]
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
