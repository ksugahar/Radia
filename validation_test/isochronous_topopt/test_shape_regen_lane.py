"""Density design -> watertight STL -> Cubit tet/hex -> HDiv re-evaluation.

The lane freezes both halves of the manufacturing-shape bridge: readable
Radia density-to-STL code and the headless ``cubit_stl_to_vol`` MCP route.
It deliberately compares a coarse tet reference, coarse Sculpt hex, and
fine Sculpt hex.  A mesh is not accepted merely because Cubit wrote files:
the Netgen boundary must equal the Gmsh topological skin, volume closure and
orientation must pass, two solvers must agree, and the adjoint load must
match an independently sampled demagnetizing field.

The direct/free-sideset exporter regression lives in
``validation_test/cubit/test_free_sideset_export.py``.  It guards the 2026
incident where a Sculpt volume with no ``surfaceelements`` looked solver
ready but removed every surface charge from the HDiv formulation.
"""
import json
import os
import time
from datetime import datetime, timezone
from math import cos, pi, sin
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
radia = pytest.importorskip("radia")
pytest.importorskip("trimesh")
pytest.importorskip("skimage")
pytest.importorskip("netCDF4")

from netgen.occ import Cylinder, HalfSpace, OCCGeometry, Pnt, Vec, Z  # noqa: E402
from ngsolve import HDiv, Integrate, Mesh, SetNumThreads, TaskManager  # noqa: E402

from radia.isochronous_topopt import (  # noqa: E402
    MU0, DensityAdjointVIM, HelmholtzFilter, density_to_s,
    field_functional_load, gradient_pair_points, optimize_density,
    orbit_arc_points, uniform_field_load, verify_design_iron_only)
from radia.topopt_cad import (  # noqa: E402
    iso_stl_from_grid, nodal_from_element_density, write_levelset_exodus)

RESULTS = Path(__file__).with_name("results_shape_regen_lane.json")


def _cubit_available():
    try:
        from radia_mcp.cubit import session as _cs  # noqa: F401
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(not _cubit_available(),
                                reason="radia_mcp cubit server not importable")


def make_sector_mesh(R1=0.05, R2=0.15, angle_deg=60.0, z0=0.02, thick=0.03,
                     maxh=0.02):
    ang = angle_deg * pi / 180.0
    ring = (Cylinder(Pnt(0, 0, z0), Z, r=R2, h=thick)
            - Cylinder(Pnt(0, 0, z0 - 0.01), Z, r=R1, h=thick + 0.02))
    hs1 = HalfSpace(Pnt(0, 0, 0), Vec(0, -1, 0))
    hs2 = HalfSpace(Pnt(0, 0, 0), Vec(-sin(ang), cos(ang), 0))
    return Mesh(OCCGeometry(ring * hs1 * hs2).GenerateMesh(maxh=maxh))


@pytest.fixture(scope="module")
def regen(tmp_path_factory):
    out = tmp_path_factory.mktemp("shape_regen")
    timings = {}
    SetNumThreads(4)
    chi_iron = 1000.0
    span = (pi / 12, pi / 2 - pi / 12)
    obj_pts, obj_radial = orbit_arc_points(0.115, 0.0, 7, span=span)
    pair_pts, pair_wts = gradient_pair_points(
        obj_pts, np.full(len(obj_pts), 1.0 / len(obj_pts)), delta=0.01,
        direction=obj_radial)

    def state_builder(fes):
        return uniform_field_load(fes, (0.0, 0.0, 1.0e5))

    def objective_builder(fes):
        return field_functional_load(fes, pair_pts, pair_wts, axis=2,
                                     scale=MU0, bonus_intorder=10)

    def constraint_builder(radius):
        def build(fes):
            cpts, _ = orbit_arc_points(radius, 0.0, 7, span=span)
            return field_functional_load(
                fes, cpts, np.full(len(cpts), 1.0 / len(cpts)), axis=2,
                scale=MU0, bonus_intorder=10)
        return build

    con_builders = [constraint_builder(r) for r in (0.08, 0.10)]
    started = time.perf_counter()
    with TaskManager():
        mesh = make_sector_mesh()
        fes = HDiv(mesh, order=1)
        prob = DensityAdjointVIM(fes, eps=1e-7)
        f_state = state_builder(fes)
        f_obj = objective_builder(fes)
        cons = [b(fes) for b in con_builders]
        filt = HelmholtzFilter(mesh, radius=0.012)
        lin0 = prob.linearize(
            density_to_s(filt.apply(np.full(prob.n_el, 0.5)), chi_iron),
            f_state, [f_obj] + cons)
        targets = [float(v) for v in lin0.values[1:]]
        result = optimize_density(prob, f_state, f_obj, cons, targets,
                                  chi_iron=chi_iron, volume_fraction=0.5,
                                  density_filter=filt, move_limit=0.1,
                                  max_iterations=30)
        verification = verify_design_iron_only(
            prob, result.density, state_builder,
            [objective_builder] + con_builders, chi_iron=chi_iron,
            density_filter=filt, gram_kwargs=dict(eps=1e-7))
    timings["density_design"] = time.perf_counter() - started

    started = time.perf_counter()
    rho_f = np.clip(filt.apply(result.density), 0.0, 1.0)
    nodal = nodal_from_element_density(mesh, rho_f)

    # Cubit meshes STL geometry TO THE FACETS, so the facet count is
    # the solver-mesh-size lever: decimate to ~1500 faces to keep the
    # re-evaluation solve inside the lane budget (an undecimated
    # resolution-96 surface forced 47k tets)
    stl = out / "design.stl"
    stl_info = iso_stl_from_grid(mesh, nodal, stl, level=0.5,
                                 resolution=96, smooth_iterations=3,
                                 target_faces=1500)
    import trimesh

    source_surface = trimesh.load_mesh(stl, process=True)
    coarse_surface = source_surface.simplify_quadric_decimation(face_count=300)
    coarse_surface.remove_unreferenced_vertices()
    coarse_surface.fix_normals()
    coarse_stl = out / "design_300.stl"
    coarse_surface.export(coarse_stl)
    coarse_stl_info = {
        "n_faces": int(len(coarse_surface.faces)),
        "n_vertices": int(len(coarse_surface.vertices)),
        "volume": float(abs(coarse_surface.volume)),
        "watertight": bool(coarse_surface.is_watertight),
    }
    assert coarse_stl_info["watertight"] is True
    timings["surface_generation"] = time.perf_counter() - started

    from radia_mcp.cubit.server import cubit_stl_to_vol

    mesh_specs = {
        "tet_reference": (coarse_stl, "tet", 0.01, 0.05),
        "hex_coarse": (coarse_stl, "hex", 0.01, 0.12),
        "hex_fine": (stl, "hex", 0.0, 0.03),
    }
    mesh_results = {}
    for name, (mesh_stl, scheme, size, closure_gate) in mesh_specs.items():
        started = time.perf_counter()
        mesh_results[name] = json.loads(cubit_stl_to_vol(
            stl_path=str(mesh_stl), scheme=scheme, size=size,
            closure_tolerance=closure_gate,
            out_vol=str(out / f"{name}.vol"),
            out_msh=str(out / f"{name}.msh")))
        timings[f"mesh_{name}"] = time.perf_counter() - started
    return SimpleNamespace(out=out, mesh=mesh, nodal=nodal,
                           chi_iron=chi_iron,
                           state_builder=state_builder,
                           objective_builder=objective_builder,
                           verification=verification,
                           stl=stl, stl_info=stl_info,
                           coarse_stl=coarse_stl,
                           coarse_stl_info=coarse_stl_info,
                           mesh_results=mesh_results,
                           timings=timings)


def test_grid_iso_stl_is_watertight_with_small_drift(regen):
    assert regen.stl_info["watertight"] is True
    # measured -0.18 % at 3 Taubin iterations (resolution 96); gate 2 %
    assert abs(regen.stl_info["smoothing_volume_drift"]) < 0.02
    assert abs(regen.stl_info["decimation_volume_drift"]) < 0.02


def test_cubit_iso_blocks_union_watertight(regen):
    """The official ATO block semantics on a boundary-touching design."""
    import trimesh
    from netCDF4 import Dataset
    from radia_mcp.cubit.server import _run_batch

    exo = regen.out / "design_lsd.exo"
    write_levelset_exodus(regen.mesh, regen.nodal, exo, level=0.5)
    iso_e = str(regen.out / "design_iso.e").replace(os.sep, "/")
    r = _run_batch(None, [
        "set dev on",
        f'import mesh "{str(exo).replace(os.sep, "/")}" nodal_var "LSD" '
        "no_geom",
        'create tri iso tet all nodal_var "LSD"',
        f'export mesh "{iso_e}" block all overwrite',
    ], timeout_s=900)
    assert r.get("status") == "ok", r

    ds = Dataset(iso_e.replace("/", os.sep), "r")
    try:
        xyz = np.stack([ds.variables["coordx"][:], ds.variables["coordy"][:],
                        ds.variables["coordz"][:]], axis=1)
        ids = list(np.asarray(ds.variables["eb_prop1"][:]))
        tri = {}
        for i, bid in enumerate(ids, start=1):
            v = ds.variables.get(f"connect{i}")
            et = getattr(v, "elem_type", "") if v is not None else ""
            if et.upper().startswith("TRI"):
                tri[int(bid)] = np.asarray(v[:], dtype=int) - 1
    finally:
        ds.close()
    assert len(tri) >= 2, f"expected 2 new TRI blocks, got {sorted(tri)}"
    new_ids = sorted(tri)[-2:]
    faces = np.vstack([tri[b] for b in new_ids])
    m = trimesh.Trimesh(vertices=xyz, faces=faces, process=True)
    m.merge_vertices()
    m.remove_unreferenced_vertices()
    m.fix_normals()
    assert m.is_watertight, (len(m.faces), new_ids)
    assert abs(m.volume) > 0.0


def test_stl_to_vol_gates(regen):
    for name, result in regen.mesh_results.items():
        assert result["status"] == "ok", (name, result)
        assert result["gates"]["no_inverted_elements"] is True
        assert result["gates"]["closure_ok"] is True
        assert result["gates"]["boundary_faces_ok"] is True, result
        assert result["vol_boundary_faces"] == result["msh_skin_faces"] > 0
        assert Path(result["vol"]).is_file()

    # Sculpt cell refinement must improve the geometry-volume closure.
    assert (regen.mesh_results["hex_fine"]["closure"]
            < regen.mesh_results["hex_coarse"]["closure"])


def test_functional_reevaluation_on_regenerated_body(regen):
    """Cross-check the regenerated body with two solvers and two meshes."""
    chi_eval = 100.0

    def evaluate(name):
        started = time.perf_counter()
        mesh = Mesh(str(regen.out / f"{name}.vol"))
        fes = HDiv(mesh, order=1)
        problem = DensityAdjointVIM(fes, eps=1e-7)
        state = regen.state_builder(fes)
        objective = regen.objective_builder(fes)
        s = np.full(problem.n_el, 1.0 / chi_eval)
        field_native, it_native = problem.solve(
            s, state, tol=1e-8, maxiter=5000, solver="native")
        field_ngsolve, it_ngsolve = problem.solve(
            s, state, tol=1e-8, maxiter=5000, solver="ngsolve-cg")
        j_native = float(ng.InnerProduct(objective.vec, field_native.vec))
        j_ngsolve = float(ng.InnerProduct(objective.vec, field_ngsolve.vec))
        obj_pts, obj_radial = orbit_arc_points(
            0.115, 0.0, 7, span=(pi / 12, pi / 2 - pi / 12))
        pair_pts, pair_wts = gradient_pair_points(
            obj_pts, np.full(len(obj_pts), 1.0 / len(obj_pts)),
            delta=0.01, direction=obj_radial)
        from radia.isochronous_topopt import demag_field_from_solution

        sampled = demag_field_from_solution(
            problem.demag, field_native, pair_pts)
        j_direct = float(MU0 * np.dot(pair_wts, sampled[:, 2]))
        return {
            "ne": int(mesh.ne),
            "ndof": int(fes.ndof),
            "native_iterations": int(it_native),
            "ngsolve_iterations": int(it_ngsolve),
            "J_native": j_native,
            "J_ngsolve": j_ngsolve,
            "J_direct_field": j_direct,
            "solver_parity_relative": (
                abs(j_native - j_ngsolve) / max(abs(j_ngsolve), 1e-300)),
            "reciprocity_relative": (
                abs(j_native - j_direct) / max(abs(j_direct), 1e-300)),
            "elapsed_s": time.perf_counter() - started,
            "problem": problem,
        }

    with TaskManager():
        evaluations = {
            name: evaluate(name)
            for name in ("tet_reference", "hex_coarse", "hex_fine")
        }
        fine_problem = evaluations["hex_fine"].pop("problem")
        for name in ("tet_reference", "hex_coarse"):
            evaluations[name].pop("problem")
        demag_factor_z = float(fine_problem.demag.DemagFactor(
            ng.CoefficientFunction((0.0, 0.0, 1.0))))

        started = time.perf_counter()
        stair_mesh = regen.verification.iron_mesh
        stair_fes = HDiv(stair_mesh, order=1)
        stair_problem = DensityAdjointVIM(stair_fes, eps=1e-7)
        stair_field, stair_iterations = stair_problem.solve(
            np.full(stair_problem.n_el, 1.0 / chi_eval),
            regen.state_builder(stair_fes), tol=1e-8,
            maxiter=5000, solver="native")
        J_stair = float(ng.InnerProduct(
            regen.objective_builder(stair_fes).vec, stair_field.vec))
        stair_elapsed = time.perf_counter() - started

    J_fine = evaluations["hex_fine"]["J_native"]
    J_coarse = evaluations["hex_coarse"]["J_native"]
    J_tet = evaluations["tet_reference"]["J_native"]
    fine_vs_stair = (J_fine - J_stair) / abs(J_stair)
    fine_vs_tet = abs(J_fine - J_tet) / abs(J_tet)
    coarse_vs_tet = abs(J_coarse - J_tet) / abs(J_tet)

    def compact_mesh(result):
        return {
            key: result[key]
            for key in (
                "scheme", "size", "stl_volume", "mesh_volume", "closure",
                "closure_tolerance", "total_negative", "min_jacobian_det",
                "vol_boundary_faces", "msh_skin_faces", "by_type",
                "dof_estimate", "gates")
        }

    stl_record = {
        key: value for key, value in regen.stl_info.items()
        if key != "path"
    }

    timing_candidates = dict(regen.timings)
    timing_candidates.update({
        f"physics_{name}": values["elapsed_s"]
        for name, values in evaluations.items()
    })
    timing_candidates["physics_staircase"] = stair_elapsed
    top_timings = dict(sorted(
        timing_candidates.items(), key=lambda item: item[1], reverse=True)[:4])
    record = dict(
        schema="radia.isochronous-topopt-shape-regen/v4",
        executed_at_utc=datetime.now(timezone.utc).isoformat(),
        radia_version=str(radia.__version__),
        ngsolve_version=str(ng.__version__),
        cubit_runtime="Coreform Cubit 2025.12 headless batch",
        timing_breakdown_s=top_timings,
        stl=stl_record,
        coarse_stl=regen.coarse_stl_info,
        meshes={name: compact_mesh(result)
                for name, result in regen.mesh_results.items()},
        evaluations=evaluations,
        ne_staircase=int(stair_mesh.ne),
        chi_eval=chi_eval,
        demag_factor_z_fine_hex=demag_factor_z,
        J_staircase=J_stair,
        staircase_iterations=int(stair_iterations),
        fine_hex_vs_staircase_relative=fine_vs_stair,
        fine_hex_vs_tet_relative=fine_vs_tet,
        coarse_hex_vs_tet_relative=coarse_vs_tet,
        boundary_face_incident=(
            "fixed: direct/free Sculpt sideset faces are exported to .vol; "
            "all three meshes now require vol_boundary_faces == msh_skin_faces"),
    )
    RESULTS.write_text(json.dumps(record, indent=1), encoding="utf-8")

    for values in evaluations.values():
        assert values["solver_parity_relative"] < 1e-8, values
        assert values["reciprocity_relative"] < 1e-6, values
        assert values["native_iterations"] < 5000, values
        assert values["ngsolve_iterations"] < 5000, values
    assert demag_factor_z > 0.2, demag_factor_z
    assert np.sign(J_fine) == np.sign(J_stair)
    assert abs(fine_vs_stair) < 0.10, fine_vs_stair
    assert fine_vs_tet < 0.05, fine_vs_tet
    assert fine_vs_tet < coarse_vs_tet, (fine_vs_tet, coarse_vs_tet)
