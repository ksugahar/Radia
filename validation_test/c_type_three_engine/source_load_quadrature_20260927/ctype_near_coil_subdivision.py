"""Near-coil field: effect of subdividing the air source load near the coil.

Linear order-2 mixed total/reduced Omega on the C-type coarse Kelvin mesh.
The reduced (air) load mu0 H_s . grad v keeps the production tetrahedral rule
(absolute order 8) everywhere except on the near set -- air elements that
straddle the coil conductor surface, lie inside it, or sit within one element
size of it -- where the reference tetrahedron is split L times (Bey, 8^L
children) with the same order-8 rule on every child.  A whole-mesh order-20
load is the independent comparison.  B is sampled along three lines through
the +x straight leg of the coil and at the gap core.  Diagnostic evidence
only.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import ngsolve as ng
from ngsolve import TET, H1, BitArray, IntegrationRule, LinearForm, dx, grad

MU0 = 4.0e-7 * np.pi
_LOG = None
_T0 = time.perf_counter()


def log(*parts):
    line = f"[{time.perf_counter() - _T0:8.1f}s] " + " ".join(str(p) for p in parts)
    if _LOG is not None:
        with open(_LOG, "a", encoding="utf-8") as stream:
            stream.write(line + "\n")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def bey_children(tet):
    x0, x1, x2, x3 = tet
    m = lambda a, b: 0.5 * (a + b)
    x01, x02, x03, x12, x13, x23 = m(x0, x1), m(x0, x2), m(x0, x3), m(x1, x2), m(x1, x3), m(x2, x3)
    return [(x0, x01, x02, x03), (x01, x1, x12, x13), (x02, x12, x2, x23), (x03, x13, x23, x3),
            (x01, x02, x03, x13), (x01, x02, x12, x13), (x02, x03, x13, x23), (x02, x12, x13, x23)]


def subdivided_tet_rule(levels, order):
    """Composite reference-tetrahedron rule: 8^levels Bey children, one base rule each."""
    base = IntegrationRule(TET, int(order))
    xi = np.array([[p.point[0], p.point[1], p.point[2]] for p in base])
    w = np.array([p.weight for p in base])
    tets = [tuple(np.array(v, dtype=float) for v in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))]
    for _ in range(int(levels)):
        tets = [child for tet in tets for child in bey_children(tet)]
    points, weights = [], []
    for v0, v1, v2, v3 in tets:
        jac = np.stack([v1 - v0, v2 - v0, v3 - v0], axis=1)
        det = abs(np.linalg.det(jac))
        for p, wp in zip(xi @ jac.T + v0, w):
            points.append(tuple(float(c) for c in p))
            weights.append(float(wp * det))
    if abs(sum(weights) - 1.0 / 6.0) > 1e-13:
        raise RuntimeError("composite rule does not integrate the reference volume")
    return IntegrationRule(points, weights), len(tets)


def observation_lines():
    y = 0.13125
    xs = np.arange(0.0255, 0.1001, 0.001)
    zs = np.arange(0.0005, 0.1001, 0.001)
    return {
        "x_line_z0": np.array([[x, y, 0.0] for x in xs]),
        "x_line_z0.06": np.array([[x, y, 0.06] for x in xs]),
        "z_line_x0.0475": np.array([[0.0475, y, z] for z in zs]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging", type=Path, required=True,
                        help="directory with mesh/, run_three_engine.py, the reference JSON "
                             "and ctype_rhs_quadrature.py")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--levels", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--whole-order", type=int, default=20)
    parser.add_argument("--near-ratio", type=float, default=1.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--inverse", default="pardiso")
    args = parser.parse_args()
    global _LOG
    _LOG = args.output.with_suffix(".log")
    log("start", vars(args))
    ng.SetNumThreads(args.threads)
    ng.SetHeapSize(200 * 1000 * 1000)  # a level-2 composite rule has 8000 points per element

    import radia as rad
    from radia import kelvin_solver as ks
    from radia.kelvin_identify_ngsolve import detect_kelvin_offset

    staging = args.staging
    rte = load_module(staging / "run_three_engine.py", "_rte")
    helpers = load_module(staging / "ctype_rhs_quadrature.py", "_quad")
    reference = json.loads((staging / "c_type_20260908_linear_kelvin_after_coarse_worker_a.json")
                           .read_text(encoding="utf-8"))
    report_path = staging / "mesh" / "mesh_result.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    kelvin_vol = staging / "mesh" / "kelvin_domain.vol"
    if report["artifacts"]["kelvin_domain_vol_sha256"] != sha256(kelvin_vol):
        raise RuntimeError("kelvin_domain.vol hash differs from mesh_result.json")
    if reference["mesh_result_sha256"] != sha256(report_path):
        raise RuntimeError("reference was produced on a different mesh contract")

    rad.UtiDelAll()
    coil, coil_manifest = rte.build_coil()
    builder = helpers.build_builder()
    mesh = ng.Mesh(str(kelvin_vol))
    offset = np.asarray(detect_kelvin_offset(mesh), dtype=float)
    physical = np.asarray(report["kelvin_physical_center_m"], dtype=float)
    kelvin_center = tuple((physical + offset).tolist())
    kelvin_radius = float(report["kelvin_radius_m"])
    h_s = rad.RadiaField(coil, "h")
    order, bonus = 2, 4

    gap_points = rte.observation_points()
    gap_core = np.abs(gap_points[:, 0]) <= float(reference["gap_core_half_length_m"]) + 1e-14
    lines = observation_lines()
    signed, distance_manifest = helpers.coil_distance_model(builder, maxh=0.004)
    line_info = {}
    for name, pts in lines.items():
        d = signed(pts)
        line_info[name] = {"points_m": pts.tolist(), "signed_distance_to_conductor_m": d.tolist()}

    payload = {"schema": "radia.validation.ctype-near-coil-subdivision.v1",
               "platform_class": platform.system(), "radia_version": getattr(rad, "__version__", None),
               "radia_module": rad.__file__, "ngsolve": ng.__version__, "threads": args.threads,
               "script_sha256": sha256(__file__),
               "helpers_sha256": sha256(staging / "ctype_rhs_quadrature.py"),
               "three_engine_sha256": sha256(staging / "run_three_engine.py"),
               "mesh_result_sha256": sha256(report_path), "coil": coil_manifest,
               "controls": {"order": order, "bonus_intorder": bonus, "mu_r": 1000.0,
                            "inverse": args.inverse, "base_rule_order": 8,
                            "near_ratio": args.near_ratio, "levels": args.levels,
                            "whole_order": args.whole_order},
               "distance_model": distance_manifest, "lines": line_info,
               "acceptance": "HOLD: diagnostic, not three-engine acceptance",
               "loads": {}, "fields_T": {}}

    def save():
        args.output.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    # Near set on the air material.
    air_rows = helpers.material_element_numbers(mesh, "air")
    centroids, sizes, _, vertex_sets = helpers.element_geometry(mesh, "air")
    d_centroid = signed(centroids)
    vertex_d = np.array([signed(v) for v in vertex_sets])
    straddles = (vertex_d.min(axis=1) < 0) & (vertex_d.max(axis=1) > 0)
    near_mask = straddles | (d_centroid < 0) | (d_centroid / sizes < args.near_ratio)
    near = BitArray(mesh.ne)
    near.Clear()
    far = BitArray(mesh.ne)
    far.Clear()
    for nr, is_near in zip(air_rows, near_mask):
        (near if is_near else far).Set(int(nr))
    payload["near_set"] = {"air_elements": int(len(air_rows)), "near_elements": int(near_mask.sum()),
                           "straddling": int(straddles.sum()),
                           "inside": int(((d_centroid < 0) & ~straddles).sum())}
    save()
    log("near set", payload["near_set"])

    kelvin_trace = None
    with ng.TaskManager():
        kelvin_trace = ks.project_source_interface_potential(
            mesh, h_s, "kelvin_int", order=order, relative_tolerance=0.05)
        source = ks.project_source_total_hodge(mesh, h_s, ("iron",), order=order,
                                               bonus_intorder=bonus)
    fes_reduced = ng.Periodic(H1(mesh, order=order, definedon=mesh.Materials("air|kelvin")))
    test = fes_reduced.TestFunction()
    air = mesh.Materials("air")

    def assemble(measure, label):
        started = time.perf_counter()
        form = LinearForm(fes_reduced)
        form += MU0 * h_s * grad(test) * measure
        with ng.TaskManager():
            form.Assemble()
        log("load", label, f"{time.perf_counter() - started:.1f}s")
        return form.vec.FV().NumPy().copy()

    def solve(label, vector):
        started = time.perf_counter()
        with ng.TaskManager():
            result = ks.solve_magnetostatic_mixed_total_reduced_omega_kelvin(
                mesh, h_s, source["potential"], kelvin_radius, kelvin_center,
                mu_r_by_material={"iron": 1000.0}, reduced_materials=("air",),
                total_materials=("iron", "kelvin"), interface_boundary="iron_air_interface",
                order=order, dirichlet_bbbnd="GND", bonus_intorder=bonus, inverse=args.inverse,
                kelvin_mats=("kelvin",), kelvin_match_exact=True,
                kelvin_interface_boundary="kelvin_int",
                kelvin_source_potential=kelvin_trace["potential"],
                total_source_h=source["harmonic_field"], total_source_materials=("iron",),
                source_rhs_reduced=vector)
        if fes_reduced.ndof != result["fes"].components[0].ndof:
            raise RuntimeError("standalone reduced space does not match the solver's")
        fields = {"gap": rte.median_plane_projection(
            gap_points, rte.evaluate_cf(result["B_cf"], mesh, gap_points))}
        for name, pts in lines.items():
            fields[name] = rte.evaluate_cf(result["B_cf"], mesh, pts)
        if not all(np.isfinite(v).all() for v in fields.values()):
            raise RuntimeError(f"{label}: non-finite field")
        payload["fields_T"][label] = {k: v.tolist() for k, v in fields.items()}
        payload["loads"].setdefault(label, {})["solve_s"] = time.perf_counter() - started
        save()
        log("solved", label)
        return fields

    base = {TET: IntegrationRule(TET, 8)}
    far_vector = assemble(dx(definedon=air, definedonelements=far, intrules=base), "far q8")
    near_q8 = assemble(dx(definedon=air, definedonelements=near, intrules=base), "near q8")
    production = assemble(dx(definedon=air, intrules=base), "whole q8")
    split_check = float(np.linalg.norm(far_vector + near_q8 - production) / np.linalg.norm(production))
    payload["near_far_split_check"] = split_check
    if split_check > 1e-12:
        raise RuntimeError(f"near + far does not reproduce the whole load: {split_check}")
    vectors = {"production_q8": production}
    for level in sorted(set(args.levels)):
        rule_l, children = subdivided_tet_rule(level, 8)
        near_l = assemble(dx(definedon=air, definedonelements=near, intrules={TET: rule_l}),
                          f"near L{level}")
        vectors[f"near_split_L{level}"] = far_vector + near_l
        payload["loads"][f"near_split_L{level}"] = {"children": children,
                                                    "points_per_element": len(rule_l)}
    vectors[f"whole_q{args.whole_order}"] = assemble(
        dx(definedon=air, intrules={TET: IntegrationRule(TET, int(args.whole_order))}),
        f"whole q{args.whole_order}")
    for label, vector in vectors.items():
        payload["loads"].setdefault(label, {})["relative_change_vs_production"] = float(
            np.linalg.norm(vector - production) / np.linalg.norm(production))

    fields = {label: solve(label, vector) for label, vector in vectors.items()}
    finest = f"near_split_L{max(args.levels)}"
    summary = {}
    for label, f in fields.items():
        row = {}
        for name in ("gap",) + tuple(lines):
            ref = np.asarray(fields[finest][name])
            val = np.asarray(f[name])
            if name == "gap":
                ref, val = ref[gap_core], val[gap_core]
            diff = np.linalg.norm(val - ref, axis=1)
            scale = np.linalg.norm(ref, axis=1)
            row[name] = {"relative_rms": float(np.sqrt(np.mean(diff ** 2)) / np.sqrt(np.mean(scale ** 2))),
                         "max_pointwise_relative": float(np.max(diff / np.maximum(scale, 1e-30))),
                         "max_abs_T": float(diff.max())}
        summary[label] = row
    payload["vs_finest_split"] = {"reference": finest, "rows": summary}
    payload["completed"] = True
    save()
    log("done summary", json.dumps(summary))


if __name__ == "__main__":
    try:
        main()
    except BaseException:
        import traceback
        log("FAILED\n" + traceback.format_exc())
        raise
