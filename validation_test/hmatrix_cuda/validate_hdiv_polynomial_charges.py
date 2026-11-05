"""Charge operators of the typed and flat HDiv routes against the polynomial cube reference; JSON record.

On a mesh of the unit cube a polynomial magnetization has a closed reference energy (polynomial_box_charges). The record
keeps four quantities apart, per route and explicit rule profile:

- ``polynomial``: Q_h = u^T B^T G B u against the reference, the independent operator check on polynomial directions;
- ``operator``: difference to the finest typed profile on the whole space, in the Frobenius norm and in the energy norm
  (largest generalized eigenvalue of the difference against the mass); two routes, no independent reference;
- ``solution``: a linear solve (M/chi + N) x = f by a dense factorization, as coefficient and as field (mass norm)
  differences, with its residual, so that solver error cannot stand in for operator error;
- ``mass``: the projection error of each polynomial mode, which bounds how far the mode lies outside the space.

The rational pyramid charge has no polynomial reference. Its energies on one pyramid are compared with the reduction
of rational_pyramid_reference, and on a mesh the record holds the mutual convergence only.
Limits are not applied here: the record reports actual values. Timings are diagnostic.
"""
from pathlib import Path
import argparse
import json
import os
import platform
import time

import numpy as np
import scipy.linalg as sl

from hdiv_mixed_meshes import build_mesh, mixed_grid, oriented, pyramid_cube, tet_cube, transition_cube
from polynomial_box_charges import check_reference, magnetization_modes, mode_energies
from rational_pyramid_reference import energies, separated_integral

# (outer, pyramid, graded near, plain near of cell targets, inner switch, far levels). The first row is the proposed
# profile; the others keep the graded rule, the near inner and one far rule, so the series converges to a common reference.
TYPED_PROFILES = {1: [(5, 6, 10, 8, .25, ((4, 1.1), (3, 1.7))), (5, 6, 10, 0, 0., ()), (8, 8, 12, 0, 0., ()), (10, 10, 16, 0, 0., ()), (12, 10, 20, 0, 0., ()), (16, 14, 30, 0, 0., ())],
                  2: [(6, 14, 10, 8, .25, ((5, 1.1), (4, 1.6), (3, 3.4))), (6, 14, 10, 0, 0., ()), (8, 14, 12, 0, 0., ()), (6, 30, 16, 0, 0., ()), (10, 30, 24, 0, 0., ()), (16, 40, 30, 0, 0., ())]}
FLAT_PROFILES = [dict(), dict(intorder=20, far_quad=8), dict(intorder=30, far_quad=10), dict(intorder=40, far_quad=12), dict(intorder=60, far_quad=14),
                 dict(intorder=40, far_quad=12, eps=1e-14), dict(intorder=60, far_quad=14, eps=1e-14)]
MESHES = dict(tets=tet_cube, pyramids=pyramid_cube, transition=transition_cube, grid=lambda: mixed_grid((2, 2, 2)))


def typed_route(mesh, order, profile):
    import ngsolve as ng
    from radia.vim import _vim as V
    rules = V._mixed_gram_rules
    V._mixed_gram_rules = lambda p: profile
    try:
        with ng.TaskManager():
            fes = ng.HDiv(mesh, order=order)
            B, G, M, _ = V._build_charge_gram_mixed(fes, build_hmatrix=False)
            B = B.toarray()
            charges = range(B.shape[0])
            dense = np.array([[G.entry(i, j) for j in charges] for i in charges])
    finally:
        V._mixed_gram_rules = rules
    return B.T @ dense @ B, M.toarray()


def flat_route(mesh, order, options):
    import ngsolve as ng
    from radia.vim import _vim as V
    if options:
        options = dict(dict(eps=1e-11, ho_far_factor=1e30), **options)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        B, G, M = V.build_charge_gram(fes, **options)[:3]
        B = B.toarray()
        action = np.column_stack([G.matvec(B[:, j]) for j in range(B.shape[1])])
    return B.T @ action, (M.toarray() if hasattr(M, "toarray") else np.asarray(M))


def modes_on_mesh(mesh, order, modes):
    """Mass projection of each mode with a long rule, and its L2 distance to the mode."""
    import ngsolve as ng
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        u, v = fes.TnT()
        mass = ng.BilinearForm(fes)
        mass += u*v*ng.dx(bonus_intorder=12)
        mass.Assemble()
        rows, columns, values = mass.mat.COO()
        dense = np.zeros((fes.ndof, fes.ndof))
        dense[np.asarray(rows), np.asarray(columns)] = np.asarray(values)
        factor = sl.cho_factor(dense)
        coordinates = (ng.x - .5, ng.y - .5, ng.z - .5)
        out = np.zeros((fes.ndof, len(modes)))
        distance = []
        for k, (component, exponents) in enumerate(modes):
            scalar = ng.CF(1.)
            for coordinate, exponent in zip(coordinates, exponents):
                for _ in range(exponent):
                    scalar = scalar*coordinate
            field = ng.CF(tuple(scalar if c == component else ng.CF(0.) for c in range(3)))
            load = ng.LinearForm(fes)
            load += ng.InnerProduct(field, v)*ng.dx(bonus_intorder=12)
            load.Assemble()
            out[:, k] = sl.cho_solve(factor, load.vec.FV().NumPy().copy())
            gf = ng.GridFunction(fes)
            gf.vec.FV().NumPy()[:] = out[:, k]
            difference = gf - field
            distance.append(float(np.sqrt(ng.Integrate(ng.InnerProduct(difference, difference), mesh, order=2*order + 12))))
    return out, distance


def polynomial_errors(N, U, reference):
    actual = U.T @ N @ U
    difference = actual - reference
    diagonal = np.abs(np.diag(difference))/np.maximum(np.abs(np.diag(reference)), 1e-300)
    return dict(relative_frobenius=float(np.linalg.norm(difference)/np.linalg.norm(reference)),
                largest_entry=float(np.abs(difference).max()), reference_largest_entry=float(np.abs(reference).max()),
                largest_diagonal_relative=float(diagonal.max()),
                uniform_factor_error=[float(abs(actual[k, k] - reference[k, k])) for k in uniform_modes(U.shape[1])])


def uniform_modes(count):
    """The constant modes come first per component of magnetization_modes."""
    return [0, count//3, 2*count//3]


def solve(N, M, chi, field_load):
    matrix = M/chi + N
    x = sl.solve(matrix, field_load, assume_a="sym")
    return x, float(np.linalg.norm(matrix @ x - field_load)/np.linalg.norm(field_load))


def difference_record(a, b, M):
    return dict(relative_frobenius=float(np.linalg.norm(a - b)/np.linalg.norm(b)),
                energy_norm=float(np.abs(sl.eigh(.5*(a - b + (a - b).T), M, eigvals_only=True)).max()),
                asymmetry=float(np.linalg.norm(a - a.T)/np.linalg.norm(a)))


def rational_record():
    """Gram entries of one pyramid whose reference coordinates are physical, per profile, against the reduction."""
    import ngsolve as ng
    from radia.vim import _vim as V
    points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1)]
    mesh = build_mesh(points, [oriented(points, [0, 1, 2, 3, 4])])
    if [tuple(mesh[v].point) for v in mesh[ng.ElementId(ng.VOL, 0)].vertices] != points:
        raise ValueError("The element no longer keeps the reference vertex order")
    reference = energies()
    finer = energies(heights=60, panel=18)
    rows = []
    rules = V._mixed_gram_rules
    try:
        for profile in TYPED_PROFILES[2]:
            V._mixed_gram_rules = lambda p: profile
            with ng.TaskManager():
                _, G, _, _ = V._build_charge_gram_mixed(ng.HDiv(mesh, order=2), build_hmatrix=False)
            actual = dict(constant=4*np.pi*G.entry(0, 0), rational=4*np.pi*G.entry(4, 4), mixed=4*np.pi*G.entry(0, 4))
            rows.append(dict(profile=profile, **{key: abs(actual[key] - finer[key])/finer[key] for key in finer}))
            print("one pyramid", profile, rows[-1], flush=True)
    finally:
        V._mixed_gram_rules = rules
    return dict(reference={key: float(value) for key, value in finer.items()},
                reduction_change=max(abs(reference[key] - finer[key])/finer[key] for key in finer),
                separated_at_one_half=separated_integral(.5), rows=rows)


def run_mesh(name, order, reference, modes, routes, chi):
    import ngsolve as ng
    mesh = MESHES[name]()
    U, distance = modes_on_mesh(mesh, order, modes)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        load = ng.LinearForm(fes)
        load += ng.InnerProduct(ng.CF((0, 0, 1)), fes.TestFunction())*ng.dx(bonus_intorder=12)
        load.Assemble()
    field_load = load.vec.FV().NumPy().copy()
    rows = []
    results = []
    for label, route, profile in routes:
        start = time.perf_counter()
        N, M = route(mesh, order, profile)
        seconds = time.perf_counter() - start
        results.append((label, profile, N, M, seconds))
    finest = max((r for r in results if r[0].startswith("typed")), key=lambda r: r[1][2])
    N_best, M_best = finest[2], finest[3]
    x_best, _ = solve(N_best, M_best, chi, field_load)
    for label, profile, N, M, seconds in results:
        x, residual = solve(N, M, chi, field_load)
        change = x - x_best
        rows.append(dict(route=label, profile=profile, ndof=int(N.shape[0]), build_seconds_diagnostic=seconds,
                         polynomial=polynomial_errors(N, U, reference),
                         operator=difference_record(N, N_best, M_best),
                         solution=dict(chi=chi, relative_residual=residual,
                                       coefficient=float(np.linalg.norm(change)/np.linalg.norm(x_best)),
                                       field=float(np.sqrt(change @ M_best @ change/(x_best @ M_best @ x_best)))),
                         mass_vs_finest=float(np.linalg.norm(M - M_best)/np.linalg.norm(M_best))))
        print(name, order, label, profile, rows[-1]["polynomial"]["relative_frobenius"], rows[-1]["operator"]["energy_norm"],
              rows[-1]["solution"]["field"], rows[-1]["solution"]["coefficient"], round(seconds, 1), flush=True)
    return dict(mesh=name, order=order, largest_mode_distance=float(max(distance)), rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--orders", type=int, nargs="+", choices=[1, 2], default=[1, 2])
    parser.add_argument("--meshes", nargs="+", choices=sorted(MESHES), default=sorted(MESHES))
    parser.add_argument("--chi", type=float, default=10.)
    parser.add_argument("--gauss", type=int, nargs=2, default=[24, 32])
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise FileExistsError(args.out)
    import ngsolve as ng
    result = dict(schema="radia.hdiv-polynomial-charges.v1", accepted=False, adoption="diagnostic", units="SI", gauss=args.gauss,
                  closed_forms=check_reference(args.gauss[0]), ngsolve_version=ng.__version__, machine=platform.machine(),
                  source_sha=os.environ.get("RADIA_SOURCE_SHA"), cases=[])
    try:
        if 2 in args.orders:
            result["rational"] = rational_record()
        for order in args.orders:
            modes = magnetization_modes(order)
            reference = mode_energies(modes, n=args.gauss[0])
            finer = mode_energies(modes, n=args.gauss[1])
            reference_change = float(np.abs(finer - reference).max()/np.abs(finer).max())
            for name in args.meshes:
                routes = [("typed", typed_route, profile) for profile in TYPED_PROFILES[order]]
                if name == "tets":
                    routes += [("flat", flat_route, options) for options in FLAT_PROFILES]
                case = run_mesh(name, order, finer, modes, routes, args.chi)
                case["reference_gauss_change"] = reference_change
                result["cases"].append(case)
                args.out.write_text(json.dumps(result, indent=1, allow_nan=False), encoding="utf-8")
        result["completed"] = True
    except BaseException as exc:
        result["error"] = repr(exc)
        raise
    finally:
        args.out.write_text(json.dumps(result, indent=1, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
