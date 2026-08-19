"""Section-optics design interface on a real magnet: the two gates that matter.

The fast unit tests check the design model against closed forms without a
field solve.  Two of its claims cannot be checked that way, because they
are about the bridge between the solver and the design model:

ROW GATE
    The design driver is handed exact linear response rows, and the chain
    fit is a linear least squares of linear field functionals, so the
    chain coefficients have an exact row representation

        theta = P . state + theta_incident

    Building ``P`` means composing the fit pseudo-inverse with the native
    configured-field functional API, and every part of that composition
    has a convention that can be got wrong -- the ``1/(4pi)`` the native
    rows already carry, the frame axes, the metric factor on the
    transverse rows.  This project has been bitten by the first of those
    three times.  The gate compares the predicted coefficients against the
    ones the fit actually produced; anything above 1e-8 relative means the
    interface is wrong and no design result built on it means anything.

COMPOSITION IDENTITY
    ``M = M_after . M_S . M_before`` is exact by construction, but only if
    the chain's element maps really do compose the way the design model
    assumes.  The unit test checks that on a synthetic slab; this checks
    it on a curved, fringe-dominated magnet where the frame metric and the
    per-element curvature are actually doing something.

Runs the TURBO-regime example at its coarse mesh so the whole check is a
few minutes rather than the hour the design-resolution mesh costs.  The
gates are mesh-independent; only their runtime is not.
"""
import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "src" / "radia").exists())
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "docs" / "section_optics"))
sys.path.insert(0, str(REPO / "validation_test" / "ffag_topopt"))

import ngsolve as ng  # noqa: E402
from netgen.occ import OCCGeometry  # noqa: E402
from turbo_magnet import (  # noqa: E402
    AMPERE_TURNS, RIGIDITY_T_M, build_upper_half_yoke, coil_filaments,
    orbit_entrance,
)
from validation_canonical_hcurl_ctype import (  # noqa: E402
    frame_axes, sample_frame_cloud,
)

import radia as rad  # noqa: E402
from radia import _radia_pybind as _native  # noqa: E402
from radia.accelerator_lie_topopt import (  # noqa: E402
    _fourth_order_lie_map_from_vector_potential_polynomials,
)
from radia.accelerator_magnet_topopt import (  # noqa: E402
    CoilBuilderHDivSource, PlanarDesignOrbit,
)
from radia.accelerator_section_optics import (  # noqa: E402
    section_composition_defect, snap_breaks_to_section,
)
from radia.beam_canonical_hcurl import (  # noqa: E402
    CanonicalHCurlChain, graded_breaks,
)
from radia.topology_optimization import (  # noqa: E402
    solve_hdiv_mmm_active_elements,
)
from radia.vim._vim import build_charge_gram  # noqa: E402

MU0 = 4.0e-7 * np.pi
MU_R = 1000.0
HALF_WIDTH, HALF_HEIGHT = 0.012, 0.003
ELEMENTS, ORDER_X, ORDER_S, GRADE = 32, 5, 2, 6.0
ROW_GATE = 1.0e-8
COMPOSITION_GATE = 1.0e-12


def chain_design_rows(chain, x, y, s):
    """The chain fit's own least-squares design matrix, and its metric factor.

    Mirrors ``CanonicalHCurlChain.fit_frame_samples`` exactly, including
    the factor ``g`` that multiplies the transverse rows, so the
    pseudo-inverse of what this returns IS the fit operator.
    """
    index, zeta = chain._locate(s)
    offsets = np.concatenate(([0], np.cumsum(
        [element.dimension for element in chain.elements])))
    count = x.size
    design = np.zeros((3 * count, offsets[-1]))
    metric = np.zeros(count)
    for e, element in enumerate(chain.elements):
        mask = index == e
        if not np.any(mask):
            continue
        xi = x[mask] / element.half_width_m
        eta = y[mask] / element.half_height_m
        gbx, gby, bs_columns = element.b_row_columns(xi, eta, zeta[mask])
        g, _ = element._metric(xi, zeta[mask])
        where = np.flatnonzero(mask)
        columns = slice(offsets[e], offsets[e + 1])
        design[where, columns] = gbx
        design[count + where, columns] = gby
        design[2 * count + where, columns] = bs_columns
        metric[where] = g
    return design @ chain._reduced, metric


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxh", type=float, default=0.018,
                        help="iron mesh size; the gates do not depend on it")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name(
                            "validation_section_optics_rows.json"))
    options = parser.parse_args(argv)

    started = time.perf_counter()
    rad.UtiDelAll()
    ng.SetNumThreads(int(options.threads))
    loops = coil_filaments()
    groups = [(np.stack((loop[:-1], loop[1:]), axis=1),
               AMPERE_TURNS / 2.0) for loop in loops]
    source = CoilBuilderHDivSource(segment_groups=tuple(groups))
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(
            build_upper_half_yoke()).GenerateMesh(maxh=float(options.maxh)))
        fes = ng.HDiv(mesh, order=1, discontinuous=True)
        _, gram, _m = build_charge_gram(
            fes, eps=1.0e-10, leafsize=256, eta=2.0, internal_interfaces=True,
            image_masks=[4], image_signs=[-1.0])
        rhs_vector = source.assemble_hdiv_rhs(fes)
    state = solve_hdiv_mmm_active_elements(
        charge_gram=gram, fes=fes, inv_chi=1.0 / (MU_R - 1.0), rhs=rhs_vector,
        response_matrix=np.zeros((1, fes.ndof)),
        active_elements=np.ones(mesh.ne, dtype=bool))[0]
    evaluator = gram.create_field_evaluator(
        np.ascontiguousarray(state, dtype=np.float64), 32, 0.05, 256,
        500000000, 1.0e-5, 16)

    def b_batch(points):
        values = np.ascontiguousarray(
            np.asarray(points, dtype=float).reshape(-1, 3))
        return MU0 * (source.h_field(values)
                      + np.asarray(evaluator.field(values, "auto"),
                                   dtype=float) / (4.0 * np.pi))

    print(f"mesh {mesh.ne} elements, ndof {fes.ndof}; solved in "
          f"{time.perf_counter()-started:.0f} s", flush=True)

    entry, heading, exit_x = orbit_entrance()
    coil = rad.ObjCnt([rad.ObjFlmCur(loop.tolist(), AMPERE_TURNS / 2.0)
                       for loop in loops])
    positions, tangents, stations, curvature, _l, oop, _s = (
        _native.track_reference_orbit_native(
            evaluator, MU0 / (4.0 * np.pi), int(coil), False,
            RIGIDITY_T_M, entry, heading, exit_x, 5.0e-4, 1.2, 1.0e-6, 129))
    orbit = PlanarDesignOrbit(
        positions=positions, tangents=tangents,
        magnetic_rigidity=RIGIDITY_T_M, bend_axis=np.array([0.0, 0.0, 1.0]),
        path_length_stations=stations, signed_curvature_per_m=curvature)
    s_total = float(stations[-1])
    seg_mids = 0.5 * (stations[:-1] + stations[1:])
    monitor_s = np.linspace(0.0, s_total, 401)
    by_orbit = np.einsum(
        "ij,ij->i", b_batch(orbit.position_at(monitor_s)),
        np.asarray([frame_axes(orbit, sv)[1] for sv in monitor_s]))

    # The section is a geometric region.  Choose it once from the baseline
    # field, then pin it: an element break is snapped onto its start so the
    # interval is represented exactly by whole elements.
    provisional = graded_breaks(monitor_s,
                                np.abs(np.gradient(by_orbit, monitor_s)),
                                ELEMENTS, strength=GRADE)
    by_break = np.interp(provisional, monitor_s, by_orbit)
    body = np.abs(by_break) > 0.9 * np.abs(by_orbit).max()
    section_start = float(provisional[int(np.flatnonzero(body)[-1])])
    breaks, begin, snapped = snap_breaks_to_section(provisional, section_start)

    chain = CanonicalHCurlChain(
        breaks, HALF_WIDTH, HALF_HEIGHT, order_x=ORDER_X, order_s=ORDER_S,
        curvature_per_m=lambda s: -float(np.interp(s, seg_mids, curvature)))
    rng = np.random.default_rng(20260819)
    cloud = sample_frame_cloud(orbit, b_batch, (0.0, s_total), HALF_WIDTH,
                               HALF_HEIGHT, rng, ELEMENTS, 20, breaks=breaks)
    fit = chain.fit_frame_samples(*cloud)
    print(f"orbit {1e3*s_total:.2f} mm, planarity {oop:.1e} m; fit rel "
          f"{fit.relative_residual:.3e}; section from "
          f"{1e3*section_start:.2f} mm (element {begin}, break snapped "
          f"{1e6*snapped:.1f} um)", flush=True)

    # ---- gate 1: the observation rows against the fit --------------------
    x_cloud, y_cloud, s_cloud = (np.asarray(v, dtype=float) for v in cloud[:3])
    design, metric = chain_design_rows(chain, x_cloud, y_cloud, s_cloud)
    operator = np.linalg.pinv(design)
    count = x_cloud.size
    frames = np.asarray([frame_axes(orbit, sv) for sv in s_cloud])
    horizontal, vertical, tangent = frames[:, 0], frames[:, 1], frames[:, 2]
    weights = (operator[:, :count, None] * (metric[:, None] * horizontal)
               + operator[:, count:2 * count, None]
               * (metric[:, None] * vertical)
               + operator[:, 2 * count:, None] * tangent)
    points = orbit.position_at(s_cloud) + (x_cloud[:, None] * horizontal
                                           + y_cloud[:, None] * vertical)
    tick = time.perf_counter()
    rows = MU0 * np.asarray(gram.configured_field_functional_rows(
        np.ascontiguousarray(points), np.ascontiguousarray(weights)),
        dtype=float)
    incident = np.einsum("knc,nc->k", weights, MU0 * source.h_field(points))
    theta_fit, *_ = np.linalg.lstsq(chain._reduced, fit.coefficients,
                                    rcond=None)
    row_gate = float(np.max(np.abs(rows @ state + incident - theta_fit))
                     / np.max(np.abs(theta_fit)))
    print(f"\nROW GATE  {row_gate:.3e} relative "
          f"(limit {ROW_GATE:.0e}; rows built in "
          f"{time.perf_counter()-tick:.0f} s)", flush=True)

    # ---- gate 2: the composition identity on the real chain --------------
    ay, a_s, lengths, curvatures = chain.lie_element_spoly_arrays(degree=5)

    def transfer(first, last):
        return _fourth_order_lie_map_from_vector_potential_polynomials(
            ay[first:last], a_s[first:last], lengths[first:last],
            RIGIDITY_T_M,
            reference_curvature_per_m=curvatures[first:last],
            longitudinal_component="covariant",
            reference_orbit_tolerance=2.0e-2,
            parameter_jacobians=False).transfer.factorization.R

    before = transfer(0, begin)
    section = transfer(begin, chain.element_count)
    whole = transfer(0, chain.element_count)
    composition = section_composition_defect(before, section, whole)
    print(f"COMPOSITION DEFECT  {composition:.3e} "
          f"(limit {COMPOSITION_GATE:.0e})", flush=True)
    print(f"  section R_S[1,0] {section[1,0]:+.6f}, R_S[3,2] "
          f"{section[3,2]:+.6f}; whole R[1,0] {whole[1,0]:+.6f}")
    for name, matrix in (("before", before), ("section", section),
                         ("whole", whole)):
        print(f"  det {name:8s} horizontal "
              f"{np.linalg.det(matrix[:2,:2]):.9f}, vertical "
              f"{np.linalg.det(matrix[2:4,2:4]):.9f}")

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "radia_version": getattr(rad, "__version__", "unknown"),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "host": platform.node(),
        "mesh_elements": int(mesh.ne), "ndof": int(fes.ndof),
        "iron_maxh_m": float(options.maxh),
        "orbit_length_m": s_total, "planarity_m": float(oop),
        "fit_relative_residual": float(fit.relative_residual),
        "section_start_m": section_start, "section_element": int(begin),
        "break_snap_m": float(snapped),
        "row_gate_relative": row_gate, "row_gate_limit": ROW_GATE,
        "composition_defect": composition,
        "composition_limit": COMPOSITION_GATE,
        "section_R": section[:4, :4].tolist(),
        "whole_R": whole[:4, :4].tolist(),
        "wall_seconds": time.perf_counter() - started,
    }
    options.output.write_text(json.dumps(summary, indent=2) + "\n",
                              encoding="utf-8")
    print(f"\nwrote {options.output.name}; total "
          f"{time.perf_counter()-started:.0f} s")
    rad.UtiDelAll()
    failures = []
    if not (row_gate < ROW_GATE):
        failures.append(f"row gate {row_gate:.3e} >= {ROW_GATE:.0e}")
    if not (composition < COMPOSITION_GATE):
        failures.append(
            f"composition defect {composition:.3e} >= {COMPOSITION_GATE:.0e}")
    if failures:
        raise SystemExit("FAILED: " + "; ".join(failures))
    print("both gates passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
