"""Source-load quadrature on the C-type mixed total/reduced Omega route.

Load-only interventions on the linear order-2 C-type coarse mesh: the
stiffness, interface and Kelvin terms keep the production bonus while only a
coil-source load changes its quadrature.

* air:   the reduced-region volume load mu0 H_s . grad v, reassembled with an
         explicit tetrahedral rule of absolute order q and passed through
         ``source_rhs_reduced``;
* hodge: the iron total-Hodge projection of H_s at projection bonus b;
* surface_flux: the air load moved to the enclosing faces (div H_s = 0).

Each variant is compared with the production solve and with the stored
HDiv-MMM / reduced-A fields of the same mesh.  A local map records the
element-wise air/iron source-integral quadrature error against the signed
distance from the element to the coil conductor, in units of element size.
The report is diagnostic evidence, not three-engine acceptance.
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
from ngsolve import (TET, H1, IntegrationRule, LinearForm, VectorL2, dx, grad)

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


def build_builder():
    from radia.coil_builder import CoilBuilder
    radius, straight_x, straight_y = 0.0225, 0.050, 0.0625
    centre = np.array([0.0, 0.13125, 0.0])
    start = centre + np.array([0.5 * straight_x + radius, -0.5 * straight_y, 0.0])
    return (CoilBuilder(-2000.0).set_start(start).set_cross_section(0.035, 0.105)
            .add_straight(straight_y).add_arc(radius, 90.0)
            .add_straight(straight_x).add_arc(radius, 90.0)
            .add_straight(straight_y).add_arc(radius, 90.0)
            .add_straight(straight_x).add_arc(radius, 90.0))


def rule(q):
    return {TET: IntegrationRule(TET, int(q))}


def material_element_numbers(mesh, material):
    names = mesh.GetMaterials()
    return np.array([el.nr for el in mesh.Elements(ng.VOL) if names[el.index] == material])


def element_source_integrals(mesh, cf, material, q):
    """Per-element integral of a vector CF over one material with a fixed rule.

    Rows follow ``material_element_numbers``.  The order-0 VectorL2 layout is
    checked at run time against a constant field, not assumed.
    """
    space = VectorL2(mesh, order=0)
    measure = dx(definedon=mesh.Materials(material), intrules=rule(q))
    form = LinearForm(cf * space.TestFunction() * measure)
    form.Assemble()
    values = form.vec.FV().NumPy().copy()
    n = mesh.ne
    if space.ndof != 3 * n:
        raise RuntimeError("unexpected VectorL2 order-0 size")
    rows = material_element_numbers(mesh, material)
    return np.stack([values[:n], values[n:2 * n], values[2 * n:]], axis=1)[rows]


def check_vector_l2_layout(mesh, material):
    values = element_source_integrals(mesh, ng.CoefficientFunction((1.0, 2.0, 3.0)), material, 2)
    if not (np.allclose(values[:, 1], 2 * values[:, 0], rtol=1e-12, atol=0)
            and np.allclose(values[:, 2], 3 * values[:, 0], rtol=1e-12, atol=0)
            and np.all(values[:, 0] > 0)):
        raise RuntimeError("VectorL2 order-0 layout is not component-blocked")


def coil_distance_model(builder, maxh):
    """Signed distance to the conductor from a fine surface sample (negative inside)."""
    from netgen.occ import OCCGeometry
    from scipy.spatial import cKDTree
    shape = builder.to_occ()
    coil_mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))
    points = np.array([v.point for v in coil_mesh.vertices])
    boundary, centroids = set(), []
    for element in coil_mesh.Elements(ng.BND):
        ids = [v.nr for v in element.vertices]
        boundary.update(ids)
        centroids.append(points[ids].mean(axis=0))
    surface = np.vstack([points[sorted(boundary)], np.array(centroids)])
    tree = cKDTree(surface)

    def signed(xyz):
        distance, _ = tree.query(xyz)
        inside = np.array([coil_mesh(*map(float, p)).nr >= 0 for p in xyz])
        return np.where(inside, -distance, distance)

    return signed, {"surface_samples": int(len(surface)), "maxh_m": maxh,
                    "coil_mesh_elements": int(coil_mesh.ne)}


def element_geometry(mesh, material):
    wanted = set(material_element_numbers(mesh, material).tolist())
    centroids, sizes, volumes, vertex_sets = [], [], [], []
    for element in mesh.Elements(ng.VOL):
        if element.nr not in wanted:
            continue
        xyz = np.array([mesh[v].point for v in element.vertices])
        centroids.append(xyz.mean(axis=0))
        edges = [np.linalg.norm(xyz[i] - xyz[j]) for i in range(4) for j in range(i + 1, 4)]
        sizes.append(max(edges))
        volumes.append(abs(np.linalg.det(np.stack([xyz[1] - xyz[0], xyz[2] - xyz[0],
                                                   xyz[3] - xyz[0]]))) / 6.0)
        vertex_sets.append(xyz)
    return (np.array(centroids), np.array(sizes), np.array(volumes), vertex_sets)


def local_map(mesh, h_s, material, orders, reference_order, signed_distance):
    check_vector_l2_layout(mesh, material)
    with ng.TaskManager():
        reference = element_source_integrals(mesh, h_s, material, reference_order)
    log("local map reference", material, f"q={reference_order}")
    centroids, sizes, volumes, vertex_sets = element_geometry(mesh, material)
    if len(centroids) != len(reference):
        raise RuntimeError(f"{material}: element count mismatch")
    d_centroid = signed_distance(centroids)
    vertex_d = np.array([signed_distance(v) for v in vertex_sets])
    straddles = (vertex_d.min(axis=1) < 0) & (vertex_d.max(axis=1) > 0)
    ratio = d_centroid / sizes
    scale = np.linalg.norm(reference, axis=1)
    global_scale = float(np.sqrt(np.sum(scale ** 2)))
    bins = [("straddles_conductor", straddles),
            ("inside_conductor", (~straddles) & (d_centroid < 0)),
            ("d/h<0.5", (~straddles) & (d_centroid >= 0) & (ratio < 0.5)),
            ("0.5<=d/h<1", (~straddles) & (ratio >= 0.5) & (ratio < 1)),
            ("1<=d/h<2", (~straddles) & (ratio >= 1) & (ratio < 2)),
            ("2<=d/h<4", (~straddles) & (ratio >= 2) & (ratio < 4)),
            ("d/h>=4", (~straddles) & (ratio >= 4))]
    rows = {}
    for q in orders:
        with ng.TaskManager():
            values = element_source_integrals(mesh, h_s, material, q)
        log("local map", material, f"q={q}")
        error = np.linalg.norm(values - reference, axis=1)
        total = float(np.sqrt(np.sum(error ** 2)))
        per_bin = {}
        for name, mask in bins:
            count = int(mask.sum())
            if count == 0:
                per_bin[name] = {"elements": 0}
                continue
            relative = error[mask] / np.maximum(scale[mask], 1e-300)
            per_bin[name] = {
                "elements": count,
                "median_relative_element_error": float(np.median(relative)),
                "max_relative_element_error": float(relative.max()),
                "share_of_global_error_squared": float(np.sum(error[mask] ** 2) / max(total ** 2, 1e-300)),
            }
        rows[str(q)] = {"global_relative_error": total / global_scale, "bins": per_bin}
    return {"material": material, "reference_order": reference_order, "elements": int(len(reference)),
            "element_size_m": {"min": float(sizes.min()), "median": float(np.median(sizes)),
                               "max": float(sizes.max())},
            "straddling_elements": int(straddles.sum()), "orders": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-dir", type=Path, required=True)
    parser.add_argument("--three-engine", type=Path, required=True,
                        help="run_three_engine.py of the checkout that produced the reference")
    parser.add_argument("--reference", type=Path, required=True,
                        help="stored linear three-engine result on the same mesh")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--bonus", type=int, default=4)
    parser.add_argument("--air-orders", type=int, nargs="+", default=[8, 12, 16, 20])
    parser.add_argument("--hodge-bonuses", type=int, nargs="+", default=[8, 12, 16])
    parser.add_argument("--map-orders", type=int, nargs="+", default=[6, 8, 12, 16])
    parser.add_argument("--map-reference-order", type=int, default=20)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--inverse", default="sparsecholesky",
                        help="direct solver for every solve; MKL-free sparsecholesky by default")
    args = parser.parse_args()
    global _LOG
    _LOG = args.output.with_suffix(".log")
    log("start", vars(args))
    ng.SetNumThreads(args.threads)

    import radia as rad
    from radia import kelvin_solver as ks
    from radia.kelvin_identify_ngsolve import detect_kelvin_offset

    rte = load_module(args.three_engine, "_rte")
    reference = json.loads(args.reference.read_text(encoding="utf-8"))
    mesh_report_path = args.mesh_dir / "mesh_result.json"
    kelvin_vol = args.mesh_dir / "kelvin_domain.vol"
    mesh_report = json.loads(mesh_report_path.read_text(encoding="utf-8"))
    if mesh_report["artifacts"]["kelvin_domain_vol_sha256"] != sha256(kelvin_vol):
        raise RuntimeError("kelvin_domain.vol hash differs from mesh_result.json")
    if reference["mesh_result_sha256"] != sha256(mesh_report_path):
        raise RuntimeError("reference result was produced on a different mesh contract")
    if reference["mode"] != "linear":
        raise RuntimeError("reference must be a linear run")
    points = rte.observation_points()
    if not np.allclose(points, np.asarray(reference["observation_points_m"]), atol=0, rtol=0):
        raise RuntimeError("observation points differ from the reference")
    gap_core = np.abs(points[:, 0]) <= float(reference["gap_core_half_length_m"]) + 1e-14

    rad.UtiDelAll()
    coil, coil_manifest = rte.build_coil()
    builder = build_builder()
    mesh = ng.Mesh(str(kelvin_vol))
    offset = np.asarray(detect_kelvin_offset(mesh), dtype=float)
    kelvin_center = tuple((np.asarray(mesh_report["kelvin_physical_center_m"]) + offset).tolist())
    kelvin_radius = float(mesh_report["kelvin_radius_m"])
    h_s = rad.RadiaField(coil, "h")
    order, bonus = int(args.order), int(args.bonus)

    stored = {name: rte.median_plane_projection(points, np.asarray(values))
              for name, values in reference["fields_T"].items()}
    payload = {"schema": "radia.validation.ctype-source-load-quadrature.v1",
               "platform_class": platform.system(), "python": sys.version,
               "radia_version": getattr(rad, "__version__", None), "radia_module": rad.__file__,
               "ngsolve": ng.__version__, "threads": args.threads,
               "kelvin_solver_sha256": sha256(ks.__file__),
               "script_sha256": sha256(__file__), "three_engine_sha256": sha256(args.three_engine),
               "reference": str(args.reference), "reference_sha256": sha256(args.reference),
               "mesh_result_sha256": sha256(mesh_report_path), "mesh_elements": mesh.ne,
               "controls": {"order": order, "bonus_intorder": bonus, "mu_r": args.mu_r,
                            "inverse": args.inverse,
                            "source_trace_tolerance": 0.05, "coil": coil_manifest},
               "acceptance": "HOLD: diagnostic, not three-engine acceptance",
               "variants": {}}

    def save():
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    def metrics(field):
        projected = rte.median_plane_projection(points, field)
        out = {}
        for name, values in stored.items():
            out[f"vs_{name}"] = rte.relative_rms(values[gap_core], projected[gap_core])
        return projected, out

    # Production route first: the harness must reproduce it before any variant counts.
    started = time.perf_counter()
    # The same call as run_three_engine.solve_omega (linear), plus the chosen inverse.
    from radia.static_electromagnet import solve_static_electromagnet_mixed_total_reduced_omega
    with ng.TaskManager():
        production_result = solve_static_electromagnet_mixed_total_reduced_omega(
            mesh, rad.RadiaField(coil, "h"), rte.MIXED_DOMAIN, kelvin_radius, kelvin_center,
            order=order, linear_mu_r_by_material={"iron": args.mu_r}, bh_table=None,
            source_trace_tolerance=0.05, source_potential_contract="total_hodge",
            bonus_intorder=bonus, kelvin_source_h=None, inverse=args.inverse)
    production_field = rte.evaluate_cf(production_result["B_cf"], mesh, points)
    production, production_metrics = metrics(production_field)
    payload["production"] = {
        "runtime_s": time.perf_counter() - started, "gap_core": production_metrics,
        "source_trace": production_result["static_electromagnet_contract"]["source_trace"]}
    del production_result
    save()
    log("production", production_metrics)

    with ng.TaskManager():
        kelvin_trace = ks.project_source_interface_potential(
            mesh, h_s, "kelvin_int", order=order, relative_tolerance=0.05)
    hodge_cache = {}

    def hodge(b):
        if b not in hodge_cache:
            with ng.TaskManager():
                hodge_cache[b] = ks.project_source_total_hodge(
                    mesh, h_s, ("iron",), order=order, bonus_intorder=b)
        return hodge_cache[b]

    def solve(label, *, hodge_bonus=bonus, air_vector=None, reduced_source_load="volume",
              detail=None):
        started = time.perf_counter()
        source = hodge(hodge_bonus)
        with ng.TaskManager():
            result = ks.solve_magnetostatic_mixed_total_reduced_omega_kelvin(
                mesh, h_s, source["potential"], kelvin_radius, kelvin_center,
                mu_r_by_material={"iron": args.mu_r}, reduced_materials=("air",),
                total_materials=("iron", "kelvin"), interface_boundary="iron_air_interface",
                order=order, dirichlet_bbbnd="GND", bonus_intorder=bonus, inverse=args.inverse,
                kelvin_mats=("kelvin",), kelvin_match_exact=True,
                kelvin_interface_boundary="kelvin_int",
                kelvin_source_potential=kelvin_trace["potential"],
                total_source_h=source["harmonic_field"], total_source_materials=("iron",),
                source_rhs_reduced=air_vector, reduced_source_load=reduced_source_load,
                return_system=air_vector is None and label == "harness_default")
        field = rte.evaluate_cf(result["B_cf"], mesh, points)
        projected, row = metrics(field)
        row["vs_production"] = rte.relative_rms(production[gap_core], projected[gap_core])
        row["runtime_s"] = time.perf_counter() - started
        row["hodge_bonus"] = hodge_bonus
        row["iron_relative_harmonic_norm"] = source["relative_harmonic_norm"]
        if detail:
            row.update(detail)
        payload["variants"][label] = row
        save()
        log(label, {k: row[k] for k in row if k.startswith("vs_")})
        return result

    harness = solve("harness_default")
    solver_reduced_ndof = harness["fes"].components[0].ndof
    del harness
    if abs(payload["variants"]["harness_default"]["vs_production"]) > 1e-10:
        raise RuntimeError("harness does not reproduce the production route")
    # A standalone space with the solver's constructor arguments; assembling on a
    # component extracted from the product space is not supported.  Dirichlet and
    # gauge flags do not renumber DOFs, and the air_default_vector solve below
    # must reproduce production before any air variant counts.
    fes_reduced = ng.Periodic(H1(mesh, order=order, definedon=mesh.Materials("air|kelvin")))
    if fes_reduced.ndof != solver_reduced_ndof:
        raise RuntimeError("standalone reduced space does not match the solver's")

    test = fes_reduced.TestFunction()

    def air_load(measure, label=""):
        started = time.perf_counter()
        form = LinearForm(fes_reduced)
        form += MU0 * h_s * grad(test) * measure
        with ng.TaskManager():
            form.Assemble()
        log("air load assembled", label, f"{time.perf_counter() - started:.1f}s")
        return form.vec.FV().NumPy().copy()

    default_vector = air_load(dx(definedon=mesh.Materials("air"), bonus_intorder=bonus), "default")
    calibration = {}
    for q in range(2 * order + bonus - 2, 2 * order + bonus + 3):
        explicit = air_load(dx(definedon=mesh.Materials("air"), intrules=rule(q)), f"calibration q={q}")
        calibration[str(q)] = float(np.linalg.norm(explicit - default_vector)
                                    / np.linalg.norm(default_vector))
    payload["air_load_rule_calibration"] = calibration
    default_order = min(calibration, key=calibration.get)
    payload["air_load_default_absolute_order"] = int(default_order)
    save()
    solve("air_default_vector", air_vector=default_vector)
    if payload["variants"]["air_default_vector"]["vs_production"] > 1e-10:
        raise RuntimeError("the externally assembled air load does not reproduce production")

    air_vectors = {}
    for q in sorted(set(args.air_orders)):
        air_vectors[q] = air_load(dx(definedon=mesh.Materials("air"), intrules=rule(q)), f"q={q}")
        change = float(np.linalg.norm(air_vectors[q] - default_vector) / np.linalg.norm(default_vector))
        solve(f"air_q{q}", air_vector=air_vectors[q], detail={"air_load_relative_change": change})
    for b in sorted(set(args.hodge_bonuses)):
        solve(f"hodge_b{b}", hodge_bonus=b)
    top_q, top_b = max(args.air_orders), max(args.hodge_bonuses)
    solve(f"air_q{top_q}_hodge_b{top_b}", hodge_bonus=top_b, air_vector=air_vectors[top_q])
    solve("air_surface_flux", reduced_source_load="surface_flux")

    signed, distance_manifest = coil_distance_model(builder, maxh=0.004)
    payload["distance_model"] = distance_manifest
    payload["local_maps"] = {}
    for material in ("air", "iron"):
        started = time.perf_counter()
        payload["local_maps"][material] = local_map(
            mesh, h_s, material, args.map_orders, args.map_reference_order, signed)
        payload["local_maps"][material]["runtime_s"] = time.perf_counter() - started
        save()
        log("local map", material)
    payload["completed"] = True
    save()


if __name__ == "__main__":
    try:
        main()
        log("done")
    except BaseException:
        import traceback
        log("FAILED\n" + traceback.format_exc())
        raise
