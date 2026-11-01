"""Matched finite-volume PEEC / independent volume A-V FEM principle case.

An independently authored planar return path surrounds a window, with a T fin
on the upper bar and a 0.5 mm feed gap in the lower bar. The ideal port crosses
only that gap; its field is omitted by BOTH routes. This is a partial-potential
configuration, not a complete coil/source model or an adopted validation.

FEM uses no PEEC matrix or fitted material. It discretizes -Delta A = mu J
with vector H1 and div J = 0 with scalar H1 inside the conductor, where
J=-j omega sigma (A+grad Vt). This is the differential realization of the
volume-Neumann integral for this declared OPEN gap-source model. In particular,
curl-curl with a free gradient gauge and internal Dirichlet V is not used:
it can cancel the applied V with a gradient A and produce a spurious current.
The volume potential of open current need not be divergence-free; this case
does NOT claim a closed-current Coulomb-gauge Maxwell solution. A complete
source wire would have to be included in BOTH routes for that claim.
The postprocessed Cartesian vector Poisson potential is assembled separately.
Exterior, h/p and quadrature convergence remain explicit obligations.
"""
import argparse
import math
import numbers
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np


MU0 = 4e-7*np.pi
SIGMA = 5.8e7
# Millimetres. Include EVERY measurement plane in both volume decompositions.
KNOTS = ([-3., -2., -1., -.8, -.25, .25, .8, 1., 2., 3.],
         [-3., -2., -1., -.2, 0., .2, .5, 1., 1.2, 1.5], [-.25, .25])
REGIONS = ('body', 'root', 'fin_shaft', 'tip', 'outside_window')


def region_name(x, y):
    if not -.8 <= x < .8:
        return 'outside_window'
    return 'body' if y < -.2 else 'root' if y < .2 else 'fin_shaft' if y < 1.2 else 'tip'


def metal_at(x, y):
    return ((-1. < y < 0.) or (0. < y < 1.5 and abs(x) < 1.)
            or (-3. < y < -1. and abs(x) > 2.)
            or (-3. < y < -2. and abs(x) > .25))


def grid(level):
    levels = (level,)*3 if isinstance(level, int) and not isinstance(level, bool) else level
    if (not isinstance(levels, (tuple, list)) or len(levels) != 3
            or any(isinstance(n, bool) or not isinstance(n, int) or n < 1 for n in levels)):
        raise ValueError('Three positive integer axis scales or one isotropic scale required')
    edges = [np.r_[*(np.linspace(a, b, n+1)[:-1] for a, b in zip(k[:-1], k[1:])), k[-1]]*1e-3 for k, n in zip(KNOTS, levels)]
    centres = [(e[1:]+e[:-1])/2 for e in edges]
    x, y, z = np.meshgrid(*centres, indexing='ij')
    mask = ((y < 0) & (y > -.001) | (y > 0) & (abs(x) < .001)
            | (y < -.001) & (abs(x) > .002)
            | (y < -.002) & (abs(x) > .00025))
    source, sink = [], []
    for i, j, k in np.argwhere(mask):
        if centres[1][j] >= -.002:
            continue
        if np.isclose(edges[0][i+1], -.00025, rtol=0, atol=1e-15):
            source.append((int(i), int(j), int(k), 0, 1))
        if np.isclose(edges[0][i], .00025, rtol=0, atol=1e-15):
            sink.append((int(i), int(j), int(k), 0, -1))
    if not source or not sink:
        raise RuntimeError('Physical feed faces absent')
    return edges, mask, source, sink


# User-owned physical limit (direct user decision 2026-11-01, record SHA256
# ba27c3ce7fa01f91664b2d22470944da47f93d707b91c9fc6b155bc14ecf56d6),
# superseding the earlier 0.01 and 0.03 limits for these four regional
# quantities only.  Residual, parity and surrogate criteria are separate.
PHYSICAL_TARGETS = dict(R=.05, partial_L=.05, root_heat=.05, tip_heat=.05)


def physical_verdict(candidate, reference, reference_resolved):
    """Compare regional quantities with the matched independent volume FEM.

    The four quantities are real and finite (allocated partial L may be
    signed).  Relative difference is abs(candidate-reference)/abs(reference).
    A zero reference or an unresolved FEM reference (convergence or
    uncertainty not established) gives no verdict; within-limit raw values
    alone never adopt.  reference_resolved must be a bool.
    """
    if not isinstance(reference_resolved, bool):
        raise TypeError(f'reference_resolved must be bool, got {type(reference_resolved).__name__}')

    def real(name, side, value):
        if isinstance(value, bool) or not isinstance(value, numbers.Real):
            raise TypeError(f'{name}: {side} must be a real number, got {value!r}')
        value = float(value)
        if not math.isfinite(value):
            raise ValueError(f'{name}: {side} must be finite, got {value!r}')
        return value

    differences = {}
    for name, limit in PHYSICAL_TARGETS.items():
        value, ref = real(name, 'candidate', candidate[name]), real(name, 'reference', reference[name])
        if ref == 0:
            raise ValueError(f'{name}: zero reference is unresolved; no relative verdict')
        relative = abs(value - ref) / abs(ref)
        if not math.isfinite(relative):
            raise ValueError(f'{name}: nonfinite relative difference {relative!r}')
        differences[name] = dict(relative=relative, limit=limit)
    within = all(d['relative'] <= d['limit'] for d in differences.values())
    if not reference_resolved:
        status = 'NOT_ADOPTED: reference FEM convergence/uncertainty unresolved'
    else:
        status = 'PASS' if within else 'FAIL'
    return dict(status=status, within_raw_limits=within, differences=differences)


def specification():
    return dict(schema='radia.ih.fin-return-principle.v1', units='SI', sigma_S_per_m=SIGMA,
        mu_r=1., phasor='exp(+j omega t); complex peak; I=1 A',
        knots_mm=KNOTS, feed_gap_mm=.5,
        return_path='lower bar and both end bars, same copper as upper bar and fin',
        source='left lower-bar x=-0.25 mm gap face', sink='right lower-bar x=+0.25 mm gap face',
        port_model='ideal current/voltage port across gap; no source-wire magnetic field in either route',
        region_frame=dict(origin_m=[0., 0., 0.], fin_axis=[0., 1., 0.], span_axis=[1., 0., 0.],
                          root_half_width_m=.0002, tip_start_m=.0012, axial_bounds_m=[-.0008, .0008]),
        region_R='2*mean Joule loss / abs(I_peak)^2',
        region_L='Re integral_region conj(J) dot A_volume-Neumann / abs(I_peak)^2; half reciprocal mutual pairing',
        exterior='PEEC free-space Neumann; FEM finite all-Dirichlet exterior, must converge with box size',
        physical_targets=PHYSICAL_TARGETS,
        physical_verdict='NOT_EVALUATED: no converged reference yet',
        shared=['geometry', 'sigma', 'frequency', 'ports', 'regions', 'phasor/normalization'],
        independent=['current spaces', 'meshes', 'assembly', 'linear solver', 'volume-potential integration'])


def peec_case(level, frequency, order, checkpoint_folder=None):
    from radia.volume_peec import build_volume_peec
    from radia.volume_fin import solve_volume_fin, measure_volume_fin
    from radia.fin_observables import FinMeasurementRegions
    edges, mask, source, sink = grid(level)
    model = build_volume_peec(edges, mask, SIGMA, source_faces=source, sink_faces=sink)
    start = time.perf_counter()
    result = solve_volume_fin(model, frequency_hz=frequency, drive_current_peak_A=1., quadrature_order=order)
    solve_s = time.perf_counter()-start
    if checkpoint_folder is not None:
        # Preserve the completed circuit BEFORE the potentially costly region
        # pairing. This is evidence, not an implicit solver restart/fallback.
        import h5py
        write(checkpoint_folder/'solve-result.json', dict(status='CIRCUIT_COMPLETE_MEASUREMENT_PENDING',
            topology_sha256=result.topology_sha256, frequency_hz=frequency,
            cells=len(model.cell_indices), branches=len(model.branches),
            R_ohm=result.resistance_ohm, L_H=result.partial_inductance_H,
            KCL=result.kcl_relative_residual, KVL=result.kvl_relative_residual, solve_s=solve_s))
        with h5py.File(checkpoint_folder/'solve.h5', 'x') as h5:
            h5.attrs['schema'] = 'radia.volume-fin.solve-checkpoint.v1'
            h5.attrs['units'] = 'SI'
            h5.attrs['matrix_layout'] = 'row-major logical indexing'
            h5.attrs['index_base'] = 0
            h5.attrs['missing_cell_sentinel'] = -1
            h5.attrs['topology_sha256'] = result.topology_sha256
            h5.attrs['frequency_hz'] = result.frequency_hz
            h5.attrs['drive_current_peak_A'] = result.drive_current_peak_A
            h5.attrs['specification_json'] = json.dumps(specification())
            h5.attrs['runner_source_sha'] = subprocess.check_output(
                ['git', '-C', str(Path(__file__).resolve().parents[2]), 'rev-parse', 'HEAD'], text=True).strip()
            h5.attrs['runner_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            artifacts = {}
            for name in ('radia.peec_matrices', 'radia._radia_pybind', 'radia.volume_fin', 'radia.volume_peec'):
                path = Path(sys.modules[name].__file__).resolve()
                item = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                if path.suffix.lower() == '.pyd':
                    receipt = path.with_name(path.name+'.build.json')
                    item.update(receipt=str(receipt), receipt_sha256=hashlib.sha256(receipt.read_bytes()).hexdigest())
                artifacts[name] = item
            h5.attrs['actual_kernel_artifacts_json'] = json.dumps(artifacts)
            group = h5.create_group('topology')
            for name in model.__dataclass_fields__:
                value = getattr(model, name)
                if isinstance(value, np.ndarray):
                    group.create_dataset(name, data=value)
                else:
                    group.attrs[name] = value
            h5.create_dataset('branch_current_peak_A', data=result.branch_current_peak_A)
            h5.create_dataset('inductance_H', data=result.inductance_H)
        path = checkpoint_folder/'solve.h5'
        write(checkpoint_folder/'solve-checkpoint.json', dict(schema='radia.volume-fin.solve-checkpoint-receipt.v1',
            path=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), bytes=path.stat().st_size,
            status='CLOSED_HDF5_HASHED', topology_sha256=result.topology_sha256))
    measured = measure_volume_fin(model, result, FinMeasurementRegions(**specification()['region_frame']),
        quadrature_order=order, convention=specification()['region_L']+'; '+specification()['port_model'])
    return dict(level=level, cells=len(model.cell_indices), branches=len(model.branches),
        R_ohm=result.resistance_ohm, L_H=result.partial_inductance_H,
        regions=measured, solve_s=solve_s, total_s=time.perf_counter()-start,
        KCL=result.kcl_relative_residual, KVL=result.kvl_relative_residual)


def make_mesh(folder, maxh, padding):
    import ngsolve as ng
    from netgen.occ import Box, Glue, OCCGeometry, Pnt
    shapes = {}
    for xa, xb in zip(KNOTS[0][:-1], KNOTS[0][1:]):
        for ya, yb in zip(KNOTS[1][:-1], KNOTS[1][1:]):
            x, y = (xa+xb)/2, (ya+yb)/2
            if not metal_at(x, y):
                continue
            name = region_name(x, y)
            box = Box(Pnt(xa*1e-3, ya*1e-3, -.00025), Pnt(xb*1e-3, yb*1e-3, .00025))
            shapes[name] = box if name not in shapes else shapes[name]+box
    for name, shape in shapes.items():
        shape.mat('copper_'+name)
        shape.maxh = maxh
        for face in shape.faces:
            face.name = 'interface'
            x, y, z = map(float, face.center)
            if y < -.002 and abs(x+.00025) < 1e-12:
                face.name = 'source'
            if y < -.002 and abs(x-.00025) < 1e-12:
                face.name = 'sink'
    copper = None
    for shape in shapes.values():
        copper = shape if copper is None else copper+shape
    bounds = ((-.003-padding, -.003-padding, -.00025-padding),
              (.003+padding, .0015+padding, .00025+padding))
    outer_box = Box(Pnt(*bounds[0]), Pnt(*bounds[1]))
    air = outer_box-copper
    air.mat('air'); air.maxh = maxh*4
    for face in air.faces:
        centre = tuple(map(float, face.center))
        face.name = 'outer' if any(abs(centre[a]-bounds[b][a]) < 1e-12 for a in range(3) for b in (0, 1)) else 'interface'
    assembly = Glue([air, *shapes.values()])
    # Glue merges each copper/air shared face. Name the FINAL topology: naming
    # only copper faces beforehand loses the port to the air-side interface.
    for face in assembly.faces:
        x, y, z = map(float, face.center)
        if y < -.002 and abs(x+.00025) < 1e-12:
            face.name = 'source'
        elif y < -.002 and abs(x-.00025) < 1e-12:
            face.name = 'sink'
    geometry = OCCGeometry(assembly)
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=maxh*4))
    path = folder/'mesh.vol'
    mesh.ngmesh.Save(str(path))
    labels = dict(schema='radia.vol-label-contract.v1', application='fin-return-principle-v1', strict_labels=True,
        required=dict(materials=['air', *('copper_'+r for r in REGIONS)], boundaries=['outer', 'source', 'sink', 'interface']),
        allowed=dict(materials=['air', *('copper_'+r for r in REGIONS)], boundaries=['outer', 'source', 'sink', 'interface']))
    write(folder/'labels.json', labels)
    args = [sys.executable, '-B', '-m', 'cubit_mesh_export.check', str(path), '--contract', str(folder/'labels.json'),
            '--report-json', str(folder/'check-vol.json'), '--format', 'json']
    with (folder/'check-vol.log').open('wb') as log:
        rc = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc or not json.loads((folder/'check-vol.json').read_text())['passed']:
        raise RuntimeError('Mandatory check-vol failed before FE space/solve')
    # Solver consumes the checked, round-tripped file rather than in-memory CAD.
    mesh = ng.Mesh(str(path))
    # Position, normal, area and adjacent material are checked independently
    # from the raw .vol in the campaign receipt as well as these FE integrals.
    n = ng.specialcf.normal(3)
    port_geometry = {}
    for name, expected_x in (('source', -.00025), ('sink', .00025)):
        boundary = mesh.Boundaries(name)
        area = float(ng.Integrate(1., mesh, definedon=boundary))
        plane_error = float(ng.Integrate((ng.x-expected_x)**2, mesh, definedon=boundary))
        normal_error = float(ng.Integrate(1-n[0]**2, mesh, definedon=boundary))
        if abs(area-5e-7) > 1e-18 or plane_error > 1e-30 or abs(normal_error) > 1e-18:
            raise RuntimeError(f'Physical port mismatch: {name} {area=} {plane_error=} {normal_error=}')
        port_geometry[name] = dict(area_m2=area, plane_x_m=expected_x,
                                   integrated_plane_error=plane_error, normal_error=normal_error)
    write(folder/'port-integrals.json', port_geometry)
    return mesh


def checked_solve(form, load, gf, lift, free):
    import ngsolve as ng
    form.Assemble(); load.Assemble()
    residual = load.vec.CreateVector(); residual.data = load.vec-form.mat*lift
    gf.vec.data = lift + form.mat.Inverse(freedofs=free, inverse='sparsecholesky')*residual
    check = load.vec.CreateVector(); check.data = load.vec-form.mat*gf.vec
    ids = np.flatnonzero(np.asarray(free, dtype=bool))
    r = np.linalg.norm(check.FV().NumPy()[ids])/max(np.linalg.norm(residual.FV().NumPy()[ids]), 1e-300)
    if not np.isfinite(r) or r > 1e-6:
        raise RuntimeError(f'Independent FE true residual actual={r} expected<=1e-6; no fallback')
    return float(r)


def fem_case(mesh, frequency, order, integration_order):
    import ngsolve as ng
    omega = 2*np.pi*frequency
    conductor = mesh.Materials('copper_.*')
    a_space = ng.VectorH1(mesh, order=order, complex=True, dirichlet='outer')
    v_space = ng.H1(mesh, order=order+1, complex=True, definedon=conductor, dirichlet='source|sink')
    space = a_space*v_space
    (a, v), (b, w) = space.TnT()
    form = ng.BilinearForm(space, symmetric=True)
    form += ng.InnerProduct(ng.grad(a), ng.grad(b))/MU0*ng.dx
    # Air carries no current in the volume integral either. The vector Laplace
    # operator is coercive with the outer Dirichlet condition, so no fictitious
    # air conductivity or gradient-gauge regularization is required.
    form += 1j*omega*SIGMA*a*b*ng.dx(definedon=conductor)
    form += 1j*omega*SIGMA*(ng.grad(v)*b+(a+ng.grad(v))*ng.grad(w))*ng.dx(definedon=conductor)
    gf = ng.GridFunction(space)
    gf.components[1].Set(ng.CF(1.+0j), definedon=mesh.Boundaries('source'))
    lift = gf.vec.CreateVector(); lift.data = gf.vec
    load = ng.LinearForm(space)
    start = time.perf_counter()
    residual = checked_solve(form, load, gf, lift, space.FreeDofs())
    av_s = time.perf_counter()-start
    ga, gv = gf.components
    jraw = -1j*omega*SIGMA*(ga+ng.grad(gv))
    indicator = ng.GridFunction(v_space)
    indicator.Set(ng.CF(1.+0j), definedon=mesh.Boundaries('source'))
    current = -complex(ng.Integrate(jraw*ng.grad(indicator), mesh, definedon=conductor, order=integration_order))
    if not np.isfinite(current) or abs(current) == 0:
        raise RuntimeError(f'Invalid terminal current {current}')
    j = jraw/current
    # Independent measurement assembly: -Delta A = mu J. No PEEC kernels,
    # matrices or current basis. Same FE degree permits an algebraic identity
    # check against the coupled A solution, separately from discretization.
    h1 = ng.VectorH1(mesh, order=order, complex=True, dirichlet='outer')
    ap, bp = h1.TnT()
    potential = ng.BilinearForm(h1, symmetric=True)
    potential += ng.InnerProduct(ng.grad(ap), ng.grad(bp))*ng.dx
    forcing = ng.LinearForm(h1); forcing += MU0*j*bp*ng.dx(definedon=conductor)
    ga_neumann = ng.GridFunction(h1)
    zero = ga_neumann.vec.CreateVector(); zero[:] = 0
    t = time.perf_counter()
    poisson_residual = checked_solve(potential, forcing, ga_neumann, zero, h1.FreeDofs())
    heat, partial = {}, {}
    for name in REGIONS:
        region = mesh.Materials('copper_'+name)
        heat[name] = float(ng.Integrate((j*ng.Conj(j)).real/(2*SIGMA), mesh, definedon=region, order=integration_order))
        partial[name] = float(ng.Integrate((ng.Conj(j)*ga_neumann).real, mesh, definedon=region, order=integration_order))
    port_z = 1j*omega/current
    volume_z = 2*sum(heat.values())+1j*omega*sum(partial.values())
    return dict(ndof=space.ndof, ndof_A=a_space.ndof, ndof_V=v_space.ndof, ne=mesh.ne, nv=mesh.nv,
        order_A=order, order_V=order+1, quadrature=integration_order,
        R_ohm=2*sum(heat.values()), L_H=sum(partial.values()),
        mean_power_W=heat, equivalent_resistance_ohm={k: 2*v for k, v in heat.items()},
        partial_L_H=partial, terminal_impedance_ohm=[port_z.real, port_z.imag],
        port_work_relative_difference=float(abs(port_z-volume_z)/abs(port_z)),
        port_current_unscaled_A=[current.real, current.imag],
        true_relative_residual=residual, poisson_true_relative_residual=poisson_residual,
        av_s=av_s, poisson_measure_s=time.perf_counter()-t, total_s=time.perf_counter()-start,
        air_sigma_ratio=0., A_space='vector H1, volume-potential differential equation',
        region_L_scope='finite-box Dirichlet volume potential; exterior convergence pending')


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')


def native_trace(path):
    """Opt-in diagnostics; trace overhead makes these times non-benchmarks."""
    if sys.getprofile() is not None:
        raise RuntimeError('Existing profiler active; no implicit replacement')
    stream = path.open('x', encoding='utf-8', newline='\n')
    started = time.perf_counter()
    entries = 0
    methods = {'build', 'solve_branch_currents', 'solve_z_matrix'}
    def emit(event, name, **fields):
        stream.write(json.dumps(dict(event=event, name=name, elapsed_s=time.perf_counter()-started, **fields))+'\n')
        stream.flush()
        os.fsync(stream.fileno())
    def profile(frame, event, arg):
        nonlocal entries
        if event in ('call', 'return') and frame.f_code.co_name in ('solve_volume_fin', 'measure_volume_fin'):
            emit(event, frame.f_code.co_name, source_file=frame.f_code.co_filename)
        if event not in ('c_call', 'c_return', 'c_exception'):
            return
        name = getattr(arg, '__name__', '')
        if name == 'inductance_entry' and event == 'c_call':
            entries += 1
            if entries % 100000 == 0:
                emit('entry_count', name, calls=entries)
        elif name in methods:
            # pybind's c_profile argument is the underlying PyCapsule method,
            # not the user-visible bound method. Identify the actual receiver
            # at the pinned Python call site rather than guessing its capsule.
            if (frame.f_code.co_name == 'solve_volume_fin'
                    and Path(frame.f_code.co_filename).name == 'volume_fin.py'):
                receiver = frame.f_locals.get('builder' if name == 'build' else 'solver')
                owner = type(receiver).__name__
                if owner in ('PEECBuilder', 'MNASolver'):
                    emit(event, owner+'.'+name,
                         profiler_receiver_type=type(getattr(arg, '__self__', None)).__name__)
    emit('start', 'native-trace', overhead='diagnostic only, not a timing benchmark')
    sys.setprofile(profile)
    def finish():
        sys.setprofile(None)
        emit('finish', 'native-trace', entry_calls=entries)
        stream.close()
    return finish


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--route', required=True, choices=['peec', 'fem'])
    parser.add_argument('--frequency', required=True, type=float)
    parser.add_argument('--level', type=int)
    parser.add_argument('--axis-levels', type=int, nargs=3)
    parser.add_argument('--volume-order', type=int, default=64)
    parser.add_argument('--maxh-mm', type=float, default=.5)
    parser.add_argument('--padding-mm', type=float, default=3.)
    parser.add_argument('--fe-order', type=int, default=2)
    parser.add_argument('--trace-native', action='store_true')
    parser.add_argument('--checkpoint-solve', action='store_true')
    args = parser.parse_args()
    if args.route == 'peec' and (args.level is None) == (args.axis_levels is None):
        raise ValueError('Choose exactly one explicit isotropic --level or --axis-levels')
    if args.output.exists():
        raise FileExistsError(args.output)
    if args.route != 'peec' and (args.trace_native or args.checkpoint_solve):
        raise ValueError('Native tracing and circuit checkpoints are PEEC-only')
    if args.checkpoint_solve:
        import h5py  # fail at startup if the explicitly requested dependency is absent
    import ngsolve as ng
    report = dict(status='INCOMPLETE', specification=specification(), argv=sys.argv,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  route=args.route, frequency_hz=args.frequency)
    write(args.output, report)
    finish_trace = native_trace(args.output.with_name('native-trace.jsonl')) if args.trace_native else None
    if args.route == 'peec':
        try:
            with ng.TaskManager():
                report['result'] = peec_case(args.axis_levels if args.axis_levels is not None else args.level,
                    args.frequency, args.volume_order, args.output.parent if args.checkpoint_solve else None)
        finally:
            if finish_trace is not None:
                finish_trace()
    else:
        with ng.TaskManager():
            mesh = make_mesh(args.output.parent, args.maxh_mm*1e-3, args.padding_mm*1e-3)
            report['mesh_sha256'] = hashlib.sha256((args.output.parent/'mesh.vol').read_bytes()).hexdigest()
            report['result'] = fem_case(mesh, args.frequency, args.fe_order, 2*args.fe_order+10)
    report['status'] = 'DIAGNOSTIC_COMPLETE_NOT_ADOPTED'
    write(args.output, report)
    if args.route == 'fem' and report['result']['port_work_relative_difference'] > 1e-6:
        report['status'] = 'FAILED_ORIGINAL_SYSTEM_PORT_WORK'
        write(args.output, report)
        raise RuntimeError('FEM port/volume work mismatch exceeds existing 1e-6 consistency limit; no physical verdict')


if __name__ == '__main__':
    main()
