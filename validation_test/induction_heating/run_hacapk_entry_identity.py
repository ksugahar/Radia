"""Identity and product gates for scalar P1 on-demand HACApK construction.

Baseline mode saves pre-refactor dense tables; compare mode reads those tables.
Large tables are scratch inputs, never publication artifacts. The JSON contains
only self-authored geometry identities, numerical gates, and neutral runtime data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import ngsolve as ng
import numpy as np

from radia import _radia_pybind as native
from run_loop_work_ring import ring_mesh


DEGREE = 7
SINGULAR = 6
EPS = 1e-10


def nodes_for(points, triangles, curved=False):
    corners = points[triangles]
    nodes = np.concatenate([corners, .5*(corners+np.roll(corners, -1, axis=1))], axis=1)
    if curved:
        # Self-authored smooth perturbation of the three mid-edge nodes.
        nodes[:, 3:, 2] += 2e-5
    return np.ascontiguousarray(nodes)


def geometry(level, curved=False):
    _, points, triangles = ring_mesh(level)
    return points, triangles, nodes_for(points, triangles, curved)


def digest(*arrays):
    return hashlib.sha256(b''.join(a.tobytes() for a in arrays)).hexdigest()


def dense(points, triangles, nodes):
    return native._AssembleSLDL_Galerkin(points, triangles, nodes, DEGREE, SINGULAR, 1)


def entry_gate(provider, expected, indices, actual_matrices=None):
    """Two sums of identical local contributions: gamma_(2*m+4) bound."""
    maximum = 0.
    count = 0
    for i, j in indices:
        sl, dl, asl, adl, m = provider.Entry(int(i), int(j))
        if actual_matrices is not None:
            sl, dl = (matrix[i, j] for matrix in actual_matrices)
        k = 2*m+4  # two accumulation paths and the dense final reduction
        gamma = k*np.finfo(float).eps/(1-k*np.finfo(float).eps)
        for name, actual, absolute_sum, reference in zip(('SL', 'DL'), (sl, dl), (asl, adl), expected):
            difference = abs(actual-reference[i, j])
            bound = gamma*absolute_sum
            if not np.isfinite(difference) or difference > bound:
                raise AssertionError(f'{name}({i},{j}) difference {difference} exceeds contribution bound {bound}')
            maximum = max(maximum, difference/bound if bound else 0.)
        count += 1
    return dict(entries=count, maximum_roundoff_bound_fraction=float(maximum))


def product_gate(points, triangles, nodes, matrices, threads, cache_bytes=8*1024*1024):
    ng.SetNumThreads(threads)
    with ng.TaskManager():
        sl, dl, provider = native._CreateP1HACApKGeometry(points, triangles, nodes, DEGREE, SINGULAR, cache_bytes)
        reference_handles = [native.HACApKBEMManager(points, m) for m in matrices]
        for handle in [sl, dl, *reference_handles]:
            assert handle.BuildHMatrix(aca_eps=EPS, leaf_size=64, eta=2., max_rank=-1, print_level=0)
        rng = np.random.default_rng(271828)
        probes = [np.ones(len(points)), points[:, 0], points[:, 2], rng.normal(size=len(points))]
        actions, records = [], []
        for name, handle, reference in zip(('SL', 'DL'), (sl, dl), reference_handles):
            for transpose in (False, True):
                for x in probes:
                    x = np.ascontiguousarray(x/np.linalg.norm(x))
                    a = handle.MatVec(x, transpose=transpose)
                    b = reference.MatVec(x, transpose=transpose)
                    relative = np.linalg.norm(a-b)/max(np.linalg.norm(b), np.finfo(float).tiny)
                    # Declared action guard, not an induced-norm theorem for ACA.
                    assert np.isfinite(relative) and relative <= 10*EPS, (name, relative)
                    actions.append(a)
                    records.append(float(relative))
        stats = dict(SL=sl.GetStats(), DL=dl.GetStats(), provider=provider.GetStats())
        assert stats['SL']['construction_route'] == 'p1-entry-on-demand'
        assert stats['SL']['source_dense_bytes'] == 0
    return actions, dict(threads=threads, max_action_difference=max(records), stats=stats)


def run(args):
    ng.SetNumThreads(1)
    records = []
    args.tables.mkdir(parents=True, exist_ok=True)
    for level in args.levels:
        for curved in ([False, True] if level == min(args.levels) else [False]):
            points, triangles, nodes = geometry(level, curved)
            label = f'{level}'+('-curved' if curved else '')
            table = args.tables/(label+'.npz')
            started = time.perf_counter()
            ng.SetNumThreads(1)
            with ng.TaskManager():
                matrices = dense(points, triangles, nodes)
            elapsed = time.perf_counter()-started
            record = dict(level=level, faces=len(triangles), vertices=len(points), curved=curved,
                          geometry_sha256=digest(points, triangles, nodes), dense_seconds=elapsed)
            if args.mode == 'baseline':
                np.savez(table, SL=matrices[0], DL=matrices[1], geometry_sha256=record['geometry_sha256'])
            else:
                with np.load(table) as old:
                    assert old['geometry_sha256'].item() == record['geometry_sha256']
                    old_matrices = [old['SL'], old['DL']]
                with ng.TaskManager():
                    _, _, provider = native._CreateP1HACApKGeometry(points, triangles, nodes, DEGREE, SINGULAR)
                    if len(triangles) <= 2048:
                        indices = [(i, j) for i in range(len(points)) for j in range(len(points))]
                    else:
                        rng = np.random.default_rng(314159)
                        indices = [(int(i), int(j)) for i, j in rng.integers(len(points), size=(512, 2))]
                        indices += [(i, i) for i in np.linspace(0, len(points)-1, 64, dtype=int)]
                    record['old_dense_entry_gate'] = entry_gate(provider, old_matrices, indices)
                    record['new_dense_entry_gate'] = entry_gate(provider, matrices, indices)
                    # Dense-refactor table identity on ALL entries, including large levels.
                    if not all(np.array_equal(a, b) for a, b in zip(matrices, old_matrices)):
                        # Use each entry's actual contribution bound, not a global fitted atol.
                        changed = set()
                        for a, b in zip(matrices, old_matrices):
                            changed.update(zip(*np.nonzero(a != b)))
                        record['dense_refactor_changed_entry_gate'] = entry_gate(provider, old_matrices, changed, actual_matrices=matrices)
                    record['dense_refactor_bitwise_equal'] = all(np.array_equal(a, b) for a, b in zip(matrices, old_matrices))
                    if level == min(args.levels):
                        _, _, uncached = native._CreateP1HACApKGeometry(points, triangles, nodes, DEGREE, SINGULAR, 0)
                        sample = indices[::max(1, len(indices)//256)]
                        for i, j in sample:
                            assert uncached.Entry(int(i), int(j)) == provider.Entry(int(i), int(j))
                        assert uncached.GetStats()['cache_bytes'] == 0
                        record['cache_off_entries'] = len(sample)
                        padded = np.vstack([points, [[.1, .1, .1]]])
                        _, _, with_interior = native._CreateP1HACApKGeometry(padded, triangles, nodes, DEGREE, SINGULAR, 0)
                        assert with_interior.Entry(len(points), 0) == (0., 0., 0., 0., 0)
                        assert with_interior.Entry(0, len(points)) == (0., 0., 0., 0., 0)
                        invalid = nodes.copy(); invalid[0, 0, 0] += .001
                        try:
                            native._CreateP1HACApKGeometry(points, triangles, invalid, DEGREE, SINGULAR)
                        except ValueError:
                            pass
                        else:
                            raise AssertionError('Mismatched geometry accepted')
                        try:
                            provider.Entry(-1, 0)
                        except IndexError:
                            pass
                        else:
                            raise AssertionError('Invalid entry accepted')
                if not curved and level <= args.product_max_level:
                    action1, first = product_gate(points, triangles, nodes, matrices, 1)
                    action4, second = product_gate(points, triangles, nodes, matrices, 4)
                    thread_error = max(np.linalg.norm(a-b)/max(np.linalg.norm(a), np.finfo(float).tiny)
                                       for a, b in zip(action1, action4))
                    assert thread_error <= 1e-12, thread_error
                    record.update(products=[first, second], max_thread_difference=float(thread_error))
                    if level == min(args.levels):
                        action0, uncached_record = product_gate(points, triangles, nodes, matrices, 1, 0)
                        assert uncached_record['stats']['provider']['cache_bytes'] == 0
                        cache_error = max(np.linalg.norm(a-b)/max(np.linalg.norm(a), np.finfo(float).tiny)
                                          for a, b in zip(action1, action0))
                        assert cache_error <= 1e-12
                        record.update(cache_off_product=uncached_record, max_cache_difference=float(cache_error))
            records.append(record)
            print(json.dumps(record), flush=True)
    return dict(schema='radia.hacapk-p1-entry-identity.v1', mode=args.mode,
                os=platform.system(), ngsolve_version=ng.__version__, numpy_version=np.__version__,
                quadrature=dict(regular_degree=DEGREE, singular_order=SINGULAR),
                entry_bound='gamma_(2*m+4) times sum of absolute local pair contributions',
                product_action_guard=10*EPS, thread_action_guard=1e-12, records=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('baseline', 'compare'), required=True)
    parser.add_argument('--tables', type=Path, required=True)
    parser.add_argument('--levels', type=int, nargs='+', default=[32, 64, 112, 152])
    parser.add_argument('--product-max-level', type=int, default=64)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
