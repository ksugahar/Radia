"""Re-executable P1 charge-diagonal campaign; completion is not numerical adoption.

CPU: actual production ChargeGram entries, not the new diagonal arithmetic.
GPU: actual production _native loader and _P1ChargeDiagonal, no loader bypass.
All numerical criteria remain pending user approval. No implicit backend switch.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import importlib.metadata as metadata
import json
import platform
import subprocess
import sys
import time

SCHEMA = 'radia.hdiv-p1-charge-diagonal.validation.v1'
CASES = ('base', 'oblique', 'translated', 'rotated', 'scaled', 'permuted_reference')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean_source(root, expected, git):
    if subprocess.check_output([str(git), 'rev-parse', 'HEAD'], cwd=root,
                               text=True, timeout=30).strip() != expected:
        raise ValueError('Source commit mismatch')
    if subprocess.check_output([str(git), 'status', '--porcelain'], cwd=root,
                               text=True, timeout=30).strip():
        raise ValueError('Source checkout is dirty')


def import_source(root):
    if 'radia' in sys.modules:
        raise ValueError('A fresh process is required; do not mix Radia providers')
    sys.path.insert(0, str(root / 'src'))


def transformed_inputs(base, label, np):
    data = {key: value.copy() if isinstance(value, np.ndarray) else value
            for key, value in base.items()}
    transform, shift = np.eye(3), np.zeros(3)
    if label == 'oblique':
        transform = np.array([[1.2, .15, .1], [0., .8, .2], [.1, 0., 1.1]])
    elif label == 'translated':
        shift = np.array([1.25, -.75, .5])
    elif label == 'rotated':
        transform = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    elif label == 'scaled':
        transform = 2. * np.eye(3)
    elif label not in ('base', 'permuted_reference'):
        raise ValueError('Unknown geometry case')
    for name, width in [('cell_verts', 4), ('face_verts', 3)]:
        geometry = base[name].reshape(-1, width, 3) @ transform.T + shift
        if label == 'permuted_reference':
            geometry[:, [1, 2]] = geometry[:, [2, 1]]
        data[name] = np.ascontiguousarray(geometry).ravel()
    return data


def physical_payload(inputs, np):
    """Map the existing Radia charge rules, without FE/orientation reassembly.

    Weights have physical measure only. The production diagonal entry applies
    the target monomial. CPU ChargeGram independently folds it into its weights.
    """
    vertices, points, weights, offsets = [], [], [], [0]
    for i, kind in enumerate(inputs['charge_kind']):
        surface = int(kind) == 1
        table = (inputs['face_verts'].reshape(-1, 3, 3) if surface else
                 inputs['cell_verts'].reshape(-1, 4, 3))
        v = table[int(inputs['charge_host'][i])]
        ref = (inputs['ref_tri_pts'].reshape(-1, 2) if surface else
               inputs['ref_tet_pts'].reshape(-1, 3))
        rule = inputs['ref_tri_w'] if surface else inputs['ref_tet_w']
        edges = v[1:] - v[0]
        measure = (np.linalg.norm(np.cross(edges[0], edges[1])) if surface else
                   abs(np.linalg.det(edges)))
        padded = np.zeros(12)
        padded[:v.size] = v.ravel()
        vertices.extend(padded.tolist())
        points.extend((v[0] + ref @ edges).ravel().tolist())
        weights.extend((rule * measure).tolist())
        offsets.append(len(weights))
    return dict(kinds=inputs['charge_kind'].tolist(),
                exponents=inputs['charge_expo'].ravel().tolist(),
                vertices=vertices, offsets=offsets, points=points, weights=weights)


def cpu_campaign(args, snapshot):
    import numpy as np
    import ngsolve as ng
    import threadpoolctl
    from radia.vim._vim import _prepare_gpu_charge_inputs, _rp
    if Path(_rp.__file__).resolve() != args.cpu_native.resolve() or sha(_rp.__file__) != args.cpu_native_sha256:
        raise ValueError('Actual CPU provider differs from the pinned independent oracle')
    receipt = json.loads(Path(str(args.cpu_native) + '.build.json').read_bytes())
    if receipt['source_commit'] != args.cpu_source_sha or receipt['binary_sha256'] != args.cpu_native_sha256 or receipt['source_dirty']:
        raise ValueError('CPU native receipt mismatch')
    ng.SetNumThreads(1)
    cases = []
    with threadpoolctl.threadpool_limits(limits=1), ng.TaskManager():
        fes = ng.HDiv(ng.Mesh(str(snapshot)), order=1)
        B, base = _prepare_gpu_charge_inputs(fes, backend='gpu')
        if fes.ndof != 12 or B.shape != (13, 12):
            raise ValueError('This campaign requires the declared one-TET P1 fixture')
        for label in CASES:
            inputs = transformed_inputs(base, label, np)
            gram = _rp._ChargeGramHMatrix(**inputs)
            payload = physical_payload(inputs, np)
            diagonal = np.asarray([gram.entry(i, i) for i in range(13)])
            scale = np.asarray(gram.charge_sigma())
            if (diagonal.shape != (13,) or scale.shape != (13,) or
                    not np.isfinite(diagonal).all() or not np.isfinite(scale).all() or
                    (diagonal <= 0).any() or (scale <= 0).any()):
                raise ValueError('Invalid CPU oracle self diagonal/scale')
            cases.append(dict(case=label, **payload,
                              cpu_diagonal=diagonal.tolist(), cpu_scale=scale.tolist()))
    if any(name.endswith('_hdiv_cuda') for name in sys.modules):
        raise ValueError('CPU oracle unexpectedly imported a GPU module')
    return dict(cases=cases, cpu_source_sha=args.cpu_source_sha,
                cpu_native_sha256=args.cpu_native_sha256,
                cpu_receipt_sha256=sha(str(args.cpu_native) + '.build.json'),
                cpu_native_actual_path=str(Path(_rp.__file__).resolve()),
                ngsolve=ng.__version__, ngsolve_module_path=ng.__file__,
                ngsolve_distribution_version=metadata.version('ngsolve'),
                fe_dofs=12, charge_map_shape=[13, 12], threads=1,
                shared_components=['production analytic potentials and first moments'],
                independence='Independent CPU ChargeGram constructor/entry pipeline; not an independent formulation oracle')


def gpu_campaign(args, admission, oracle):
    import numpy as np
    from radia.gpu import _native
    native = _native('_hdiv_cuda')  # Keep strict same-source CPU/GPU provenance guards.
    binary = Path(native.__file__)
    receipt = json.loads(Path(str(binary) + '.build.json').read_bytes())
    if sha(binary) != admission['native_sha256'] or receipt['source_commit'] != admission['source_sha']:
        raise ValueError('GPU artifact/source mismatch')
    cases = []
    for case in oracle['cases']:
        arrays = [np.asarray(case[key], dtype=dtype) for key, dtype in
                  [('kinds', np.int32), ('exponents', np.int32), ('vertices', np.float64),
                   ('offsets', np.uint64), ('points', np.float64), ('weights', np.float64)]]
        owner = native._P1ChargeDiagonal(*arrays, order=1, backend='gpu')
        try:
            actual, scale = owner.diagnostic_diagonal(), owner.diagnostic_scale()
            expected, expected_scale = np.asarray(case['cpu_diagonal']), np.asarray(case['cpu_scale'])
            if (actual.shape != (13,) or scale.shape != (13,) or
                    not np.isfinite(actual).all() or not np.isfinite(scale).all() or
                    (actual <= 0).any() or (scale <= 0).any()):
                raise ValueError('Invalid P1 diagonal result')
            cases.append(dict(case=case['case'], actual=actual.tolist(), scale=scale.tolist(),
                              relative_error=(np.abs(actual-expected)/expected).tolist(),
                              scale_relative_error=(np.abs(scale-expected_scale)/expected_scale).tolist(),
                              statistics=owner.stats()))
        finally:
            owner.close()
        closed = owner.stats()
        if not closed['closed'] or closed['owned_device_bytes'] != 0:
            raise ValueError('P1 owner did not close')
        cases[-1]['closed_statistics'] = closed
    return dict(cases=cases, native_actual_path=str(binary), native_sha256=sha(binary),
                native_source_sha=receipt['source_commit'], native_receipt_sha256=sha(str(binary)+'.build.json'))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--backend', choices=('cpu', 'gpu'), required=True)
    for key in ('mesh', 'contract', 'receipt', 'receipt-provenance', 'out', 'git'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--receipt-sha256', required=True)
    p.add_argument('--source-root', type=Path)
    p.add_argument('--cpu-source-sha')
    p.add_argument('--cpu-native', type=Path)
    p.add_argument('--cpu-native-sha256')
    p.add_argument('--prior-oracle', type=Path)
    p.add_argument('--prior-oracle-sha256')
    p.add_argument('--oracle', type=Path)
    p.add_argument('--oracle-sha256')
    p.add_argument('--admission', type=Path)
    p.add_argument('--admission-sha256')
    p.add_argument('--nvidia-smi', type=Path)
    args = p.parse_args(argv)
    if args.out.exists():
        raise FileExistsError('Preserve existing result JSON; choose a new output')
    result = dict(schema=SCHEMA, backend=args.backend, execution_complete=False,
                  criterion_status='pending-user-approval', numerically_accepted=False,
                  scope='P1 self diagonal only; no full Gram/solver/scale/performance qualification',
                  entry_sha256=sha(__file__), python=sys.executable, platform=platform.platform())
    started = time.monotonic()
    try:
        driver_root = Path(__file__).resolve().parents[2]
        driver_sha = subprocess.check_output([str(args.git), 'rev-parse', 'HEAD'],
                                             cwd=driver_root, text=True, timeout=30).strip()
        clean_source(driver_root, driver_sha, args.git)
        result['driver_source_sha'] = driver_sha
        admission = oracle = prior = None
        if args.backend == 'cpu':
            if not all((args.source_root, args.cpu_source_sha, args.cpu_native, args.cpu_native_sha256)):
                raise ValueError('Explicit independent CPU source/native pins required')
            clean_source(args.source_root, args.cpu_source_sha, args.git)
            if sha(args.cpu_native) != args.cpu_native_sha256:
                raise ValueError('CPU native hash mismatch')
            if args.prior_oracle:
                if not args.prior_oracle_sha256 or sha(args.prior_oracle) != args.prior_oracle_sha256:
                    raise ValueError('Prior frozen oracle digest mismatch')
                prior = json.loads(args.prior_oracle.read_bytes())
                if (prior['oracle_source'] != args.cpu_source_sha or not prior['executed'] or
                        tuple(c['case'] for c in prior['cases']) != CASES):
                    raise ValueError('Prior oracle identity/case coverage mismatch')
        else:
            if not all((args.oracle, args.oracle_sha256, args.admission, args.admission_sha256, args.nvidia_smi)):
                raise ValueError('Pinned oracle, manager admission and device probe required')
            if sha(args.oracle) != args.oracle_sha256 or sha(args.admission) != args.admission_sha256:
                raise ValueError('Oracle/admission digest mismatch')
            oracle, admission = json.loads(args.oracle.read_bytes()), json.loads(args.admission.read_bytes())
            if oracle['schema'] != SCHEMA or oracle['backend'] != 'cpu' or not oracle['execution_complete'] or tuple(c['case'] for c in oracle['cases']) != CASES:
                raise ValueError('Incomplete independent CPU oracle')
            if admission['gpu_execution_approved'] is not True or admission['oracle_sha256'] != args.oracle_sha256:
                raise ValueError('GPU execution is not approved for this oracle')
            if admission['source_sha'] != driver_sha:
                raise ValueError('Driver and GPU candidate source identities differ')
            deadline = datetime.datetime.fromisoformat(admission['deadline_utc'])
            if deadline.tzinfo is None or deadline <= datetime.datetime.now(datetime.timezone.utc):
                raise ValueError('GPU admission deadline invalid or expired')
            uuid = subprocess.check_output([str(args.nvidia_smi), '--query-gpu=uuid', '--format=csv,noheader'], text=True, timeout=10).strip()
            if uuid != admission['device_uuid']:
                raise ValueError('GPU UUID mismatch')
            result.update(oracle_sha256=args.oracle_sha256, admission_sha256=args.admission_sha256, device_uuid=uuid)
            if args.source_root:
                clean_source(args.source_root, admission['source_sha'], args.git)
        if args.source_root:
            import_source(args.source_root)
        from radia.hdiv_gpu import _checked_mesh_snapshot
        provenance = json.loads(args.receipt_provenance.read_bytes())
        with _checked_mesh_snapshot(args.mesh, args.contract, args.receipt, args.receipt_sha256, provenance) as (snapshot, mesh_hash, contract_hash):
            result.update(mesh_sha256=mesh_hash, contract_sha256=contract_hash,
                          receipt_sha256=sha(args.receipt), receipt_provenance_sha256=sha(args.receipt_provenance))
            if args.backend == 'cpu':
                result.update(cpu_campaign(args, snapshot))
                if prior is not None:
                    result['prior_frozen_reference'] = dict(sha256=args.prior_oracle_sha256,
                        source_sha=prior['oracle_source'], exact_values_reproduced=
                        [(a['case']==b['case'] and a['cpu_diagonal']==b['cpu_diagonal'] and a['cpu_scale']==b['cpu_scale'])
                         for a, b in zip(result['cases'], prior['cases'])])
            else:
                if mesh_hash != oracle['mesh_sha256'] or contract_hash != oracle['contract_sha256']:
                    raise ValueError('GPU input differs from the frozen CPU oracle')
                result.update(gpu_campaign(args, admission, oracle))
        clean_source(driver_root, driver_sha, args.git)
        if args.source_root:
            clean_source(args.source_root, args.cpu_source_sha if args.backend == 'cpu' else admission['source_sha'], args.git)
        result['execution_complete'] = True
    except Exception as error:
        result['error'] = repr(error)
    result.update(elapsed_seconds=time.monotonic()-started,
                  finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return 0 if result['execution_complete'] else 1


if __name__ == '__main__':
    sys.exit(main())
