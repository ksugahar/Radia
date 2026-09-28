"""TEAM 28 (3-D HCurl eddy bubble + VIM): source-field quadrature vs gap.

The reduced system (sampled resistance, tet-volume inductance) is assembled
once on the production current basis (``intorder=10``).  The coil drive
``int mode . A_ext`` and the Lorentz force ``int J x B_ext`` are then sampled
on bases built from the same mesh, parent space and ports with a higher
``intorder``; their response vectors must equal the production ones, so only
the source quadrature changes.  The disk is moved towards the coil through the
field evaluation offset (mesh fixed).  The coil is the CoilBuilder solid
winding pack evaluated by Radia.  Diagnostic evidence only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

_LOG = None
_T0 = time.perf_counter()
DISK_BOTTOM_M = 0.0108


def log(*parts):
    line = f"[{time.perf_counter() - _T0:8.1f}s] " + " ".join(str(p) for p in parts)
    if _LOG is not None:
        with open(_LOG, "a", encoding="utf-8") as stream:
            stream.write(line + "\n")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_mesh_and_ports(maxh_m):
    import netgen.occ as occ
    import ngsolve as ng

    disk = occ.Cylinder(occ.Pnt(0.0, 0.0, DISK_BOTTOM_M), occ.Z, 0.065, 0.003)
    disk.mat("Al")
    for face in disk.faces:
        face.name = "disk_air"
    mesh = ng.Mesh(occ.OCCGeometry(disk).GenerateMesh(maxh=maxh_m))
    fes = ng.HCurl(mesh, order=6, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += (ng.curl(u) * ng.curl(v) + 100.0 * u * v) * ng.dx
    metric = ng.BilinearForm(fes)
    metric += u * v * ng.dx
    training = (
        ng.CF((-ng.y, ng.x, 0.0)),
        (ng.x * ng.x + ng.y * ng.y) * ng.CF((-ng.y, ng.x, 0.0)),
        ng.z * ng.CF((-ng.y, ng.x, 0.0)),
    )
    with ng.TaskManager():
        stiffness.Assemble()
        metric.Assemble()
        ports = []
        for vector_potential in training:
            form = ng.LinearForm(fes)
            form += vector_potential * ng.curl(v) * ng.dx
            form.Assemble()
            ports.append(form.vec.FV().NumPy().copy())
    return mesh, fes, stiffness, metric, np.column_stack(ports)


def build_basis(ng, vim, mesh, fes, stiffness, metric, ports, intorder):
    with ng.TaskManager():
        return vim.NgsolveEddyBubbleHCurlBasis(
            mesh, fes, stiffness, metric, ports, steps=1,
            conductive_materials="Al", volume_materials="Al",
            intorder=int(intorder), parent_order=6, current_gram_rtol=1.0e-10)


def force_operator(current_basis, unit_b):
    return np.transpose(np.sum(current_basis.weights[None, :, None]
                               * np.cross(current_basis.modes, np.conj(unit_b)[None, :, :]),
                               axis=1), (1, 0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True,
                        help="checkout holding validation_test/maglev (drivers are imported from it)")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maxh", type=float, default=0.025)
    parser.add_argument("--source-orders", type=int, nargs="+", default=[10, 14, 18])
    parser.add_argument("--gaps-mm", type=float, nargs="+", default=[10.8, 5.0, 2.0, 1.0, 0.5])
    parser.add_argument("--ripple-gaps-mm", type=float, nargs=3, default=[1.0, 3.0, 0.05],
                        metavar=("START", "STOP", "STEP"))
    parser.add_argument("--arc-max-segment-length", type=float, default=0.002)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    global _LOG
    _LOG = args.output.with_suffix(".log")
    log("start", vars(args))

    maglev = args.repo / "validation_test" / "maglev"
    sys.path.insert(0, str(maglev))
    import ngsolve as ng
    import radia
    import radia.vim as vim
    import team28_hcurl_vim_force as t28
    import team28_coilbuilder_eddy_bubble as t28cb
    ng.SetNumThreads(args.threads)

    current = t28.REFERENCE_COIL_CURRENT_A
    s = 2.0j * np.pi * t28.FREQUENCY_HZ
    payload = {"schema": "radia.validation.team28-source-quadrature.v1",
               "host": platform.node(), "radia_version": getattr(radia, "__version__", None),
               "radia_module": radia.__file__, "ngsolve": ng.__version__,
               "script_sha256": sha256(__file__),
               "driver_sha256": {"team28_hcurl_vim_force.py": sha256(maglev / "team28_hcurl_vim_force.py"),
                                 "team28_coilbuilder_eddy_bubble.py":
                                     sha256(maglev / "team28_coilbuilder_eddy_bubble.py")},
               "controls": {**{k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
                            "production_intorder": 10, "coil_current_A": current,
                            "frequency_hz": t28.FREQUENCY_HZ, "sigma_S_per_m": t28.SIGMA_AL,
                            "disk_bottom_m": DISK_BOTTOM_M},
               "acceptance": "HOLD: diagnostic, not TEAM 28 acceptance", "rows": []}

    def save():
        args.output.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    mesh, fes, stiffness, metric, ports = build_mesh_and_ports(args.maxh)
    bases = {}
    for q in sorted(set([10] + list(args.source_orders))):
        started = time.perf_counter()
        bases[q] = build_basis(ng, vim, mesh, fes, stiffness, metric, ports, q)
        log("basis", q, bases[q].current_basis.n_samples, f"{time.perf_counter() - started:.1f}s")
    production = bases[10]
    reference_vectors = np.asarray(production.response_basis.vectors)
    consistency = {}
    for q, basis in bases.items():
        vectors = np.asarray(basis.response_basis.vectors)
        same_shape = vectors.shape == reference_vectors.shape
        relative = (float(np.linalg.norm(vectors - reference_vectors) / np.linalg.norm(reference_vectors))
                    if same_shape else None)
        consistency[str(q)] = {"rank": int(basis.rank), "samples": int(basis.current_basis.n_samples),
                               "response_vectors_relative_difference": relative}
        if relative is None or relative > 1e-10:
            raise RuntimeError(f"intorder {q} changed the reduced modes: {consistency[str(q)]}")
    payload["basis_consistency"] = consistency
    payload["mesh"] = {"elements": int(mesh.ne), "parent_ndof": int(fes.ndof), "rank": int(production.rank)}
    save()

    started = time.perf_counter()
    with ng.TaskManager():
        interaction = production.tet_volume_interaction(
            mesh, fes, degree=5, projection_quad=7, outer_quad=4,
            projection_tolerance=1.0e-10, materials="Al")
    system = production.assemble_vim(sigma=t28.SIGMA_AL, interaction=interaction)
    log("system", f"{time.perf_counter() - started:.1f}s")

    def force_at(gap_m, q):
        basis = bases[q]
        offset = float(gap_m) - DISK_BOTTOM_M
        a, b, _ = t28cb.coilbuilder_fields(
            basis.current_basis.points, coil_current_A=current,
            arc_max_segment_length_m=args.arc_max_segment_length, height_offset_m=offset)
        rhs = vim.ExternalVectorPotentialRHS(basis.current_basis, a / current)
        # Reduced harmonic solve (R + sL) c = -s P i, Faraday drive of coil current i.
        coefficients = np.linalg.solve(system.impedance(s), -s * rhs * current)
        operator = force_operator(basis.current_basis, b / current)
        return 0.5 * np.real(operator @ coefficients * np.conj(current))

    top = max(args.source_orders)
    # The production sampling at the nominal gap must reproduce the stored lane.
    stored = json.loads((maglev / "team28_coilbuilder_eddy_bubble_results.json").read_text(encoding="utf-8"))
    stored_force = np.asarray(stored["details"]["observables"]["reference_current_force_N"])
    nominal = force_at(DISK_BOTTOM_M, 10)
    payload["nominal_vs_stored_relative"] = float(
        np.linalg.norm(nominal - stored_force) / np.linalg.norm(stored_force))
    save()
    log("nominal vs stored", payload["nominal_vs_stored_relative"])
    # The stored lane is an older release; record the drift, gate only a gross mismatch.
    if args.maxh == 0.025 and payload["nominal_vs_stored_relative"] > 1e-2:
        raise RuntimeError("production sampling departs from the stored TEAM 28 force by >1%")
    for gap_mm in args.gaps_mm:
        row = {"gap_mm": gap_mm, "force_N": {}}
        for q in sorted(bases):
            started = time.perf_counter()
            row["force_N"][str(q)] = force_at(gap_mm * 1e-3, q).tolist()
            log("gap", gap_mm, "q", q, f"Fz={row['force_N'][str(q)][2]:.9e}",
                f"{time.perf_counter() - started:.1f}s")
        ftop = row["force_N"][str(top)][2]
        row["fz_relative_to_top"] = {q: abs(v[2] - ftop) / abs(ftop) for q, v in row["force_N"].items()}
        payload["rows"].append(row)
        save()

    start, stop, step = args.ripple_gaps_mm
    gaps = np.arange(start, stop + 0.5 * step, step)
    ripple = {"gaps_mm": gaps.tolist(), "fz_N": {}}
    for q in (10, top):
        ripple["fz_N"][str(q)] = [float(force_at(g * 1e-3, q)[2]) for g in gaps]
        log("ripple sweep", q)
    f10, ftop = np.asarray(ripple["fz_N"]["10"]), np.asarray(ripple["fz_N"][str(top)])
    difference = f10 - ftop
    # Ripple: residual of each curve after a smooth (quartic) fit in the gap.
    def residual(values):
        coefficients = np.polyfit(gaps, values, 4)
        return values - np.polyval(coefficients, gaps)
    ripple["production_minus_top_relative"] = (difference / np.abs(ftop)).tolist()
    ripple["quartic_residual_relative"] = {
        "10": float(np.max(np.abs(residual(f10))) / np.max(np.abs(f10))),
        str(top): float(np.max(np.abs(residual(ftop))) / np.max(np.abs(ftop)))}
    stiff10, stifftop = np.gradient(f10, gaps * 1e-3), np.gradient(ftop, gaps * 1e-3)
    ripple["stiffness_relative_difference_max"] = float(
        np.max(np.abs(stiff10 - stifftop)) / np.max(np.abs(stifftop)))
    payload["ripple"] = ripple
    payload["completed"] = True
    save()
    log("done", json.dumps({"rows": [{"gap_mm": r["gap_mm"], **r["fz_relative_to_top"]}
                                     for r in payload["rows"]],
                            "ripple": ripple["quartic_residual_relative"],
                            "stiffness_rel_diff": ripple["stiffness_relative_difference_max"]}))


if __name__ == "__main__":
    try:
        main()
    except BaseException:
        import traceback
        log("FAILED\n" + traceback.format_exc())
        raise
