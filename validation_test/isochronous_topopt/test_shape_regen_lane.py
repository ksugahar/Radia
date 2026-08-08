"""Shape-regeneration lane: density design -> STL -> Cubit hex/tet -> re-eval.

Joins the DESIGN-side bridge (``radia.topopt_cad``) to the MESH-side cubit
MCP tool (``cubit_stl_to_vol``) on the REAL sector surrogate design of
``test_design_loop_lane.py``, and re-evaluates the objective functional on
the regenerated all-hex body.  Needs a Cubit license (batch, headless per
the driving policy); the all-hex re-evaluation is the production lane.

Measured 2026-08-08 (LAB, Cubit 2025.12), the golden anchors:

* grid iso STL at resolution 96 / 3 Taubin iterations / 1500-face
  decimation: watertight, smoothing drift -0.18 %, decimation drift
  -0.18 %;
* Cubit ``create tri iso`` on the same nodal field: the two new TRI
  blocks (fixed design-boundary caps + free iso surface) UNION to a
  watertight solid -- the official ATO block semantics confirmed on a
  design whose iron touches the domain boundary;
* ``cubit_stl_to_vol`` scheme=tet closure 0.4 %, scheme=hex (Sculpt +
  mesh-based geometry) closure 2.2 % at the smoke's sizes (hex is gated
  loosely: Sculpt respects its own cell size, not the facet density);
* the boundary-face incident this lane now locks: a bare Sculpt export
  carried ZERO ``.vol`` surface elements, so the VIM charge Gram
  silently lost all surface charge -- uniform-field J came out exactly
  demag-free (51x off, state CG 1 iteration, DemagFactor 0.000000).
  Fixed by binding mesh-based geometry before export; gated by
  ``boundary_faces_ok`` (tool), the ``vim.ChargeGram`` ZERO-boundary
  raise, and the DemagFactor lock below;
* functional re-evaluation runs on the HEX body: the 16k-tet
  facet-conforming mesh from the decimated marching-cubes STL stalls
  every native CG flavor (mass-Riesz / Jacobi / cluster-tree hit the
  20000-iteration cap even at tol=1e-6, chi=100), a recorded
  preconditioner-lane data point -- the 1.5k-hex Sculpt mesh solves the
  same functional in ~25 iterations (J_hex +0.7355 vs staircase +0.7197,
  delta +2.2 %).
"""
import json
import os
from math import cos, pi, sin
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia")
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
    return SimpleNamespace(out=out, mesh=mesh, nodal=nodal,
                           chi_iron=chi_iron,
                           state_builder=state_builder,
                           objective_builder=objective_builder,
                           verification=verification,
                           stl=stl, stl_info=stl_info)


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


@pytest.mark.parametrize("scheme,closure_gate", [("tet", 0.03),
                                                 ("hex", 0.08)])
def test_stl_to_vol_gates(regen, scheme, closure_gate):
    from radia_mcp.cubit.server import cubit_stl_to_vol

    r = json.loads(cubit_stl_to_vol(
        stl_path=str(regen.stl), scheme=scheme,
        closure_tolerance=closure_gate,
        out_vol=str(regen.out / f"regen_{scheme}.vol"),
        out_msh=str(regen.out / f"regen_{scheme}.msh")))
    assert r["status"] == "ok", r
    assert r["gates"]["no_inverted_elements"] is True
    assert r["gates"]["closure_ok"] is True
    # the boundary-face incident lock: the exported .vol must carry the
    # complete skin (a bare Sculpt free mesh exported ZERO surface
    # elements and the demag solve went silently demag-free)
    assert r["gates"]["boundary_faces_ok"] is True, r
    assert r["vol_boundary_faces"] == r["msh_skin_faces"] > 0
    assert Path(r["vol"]).is_file()


def test_functional_reevaluation_on_regenerated_body(regen):
    """The payoff: the smoothed manufacturing-shape functional.

    The comparison is GEOMETRY vs GEOMETRY (staircase iron-only extraction
    vs regenerated smooth body), so both sides are evaluated with the same
    machinery at chi=100.  The Cubit all-hex route is intentional: the same
    surface produced 16,102 tets whose facet-conforming slivers stall every
    native CG flavor (mass-Riesz / Jacobi / cluster-tree all reached the
    20,000-iteration fail-loud cap, even at tol=1e-6), while Sculpt
    produced 1,537 hexes that solve in ~25 iterations.  Tet generation
    remains a mesh gate; its solver scaling belongs to the preconditioner
    lane.

    The DemagFactor lock guards the 2026-08-08 boundary-face incident
    end-to-end: a bare-Sculpt ``.vol`` with zero surface elements made the
    charge Gram silently demag-free (DemagFactor 0.000000, J exactly the
    mass-only response, 51x off).  A healthy flat-sector body must show an
    order-1/2 demag factor.  Locked loosely: same sign, delta within
    +/-50 %.
    """
    vol = regen.out / "regen_hex.vol"
    assert vol.is_file(), "hex regeneration must run first"
    chi_eval = 100.0
    with TaskManager():
        # staircase reference: the exact-void iron-only extraction
        m1 = regen.verification.iron_mesh
        fes1 = HDiv(m1, order=1)
        prob1 = DensityAdjointVIM(fes1, eps=1e-7)
        lin1 = prob1.linearize(
            np.full(prob1.n_el, 1.0 / chi_eval),
            regen.state_builder(fes1), [regen.objective_builder(fes1)],
            tol=1e-10, solver="native")
        # regenerated smooth body
        m2 = Mesh(str(vol))
        fes2 = HDiv(m2, order=1)
        prob2 = DensityAdjointVIM(fes2, eps=1e-7)
        demag_factor_z = float(prob2.demag.DemagFactor(
            ng.CoefficientFunction((0.0, 0.0, 1.0))))
        lin2 = prob2.linearize(
            np.full(prob2.n_el, 1.0 / chi_eval),
            regen.state_builder(fes2), [regen.objective_builder(fes2)],
            tol=1e-10, solver="native")
    J_smooth = float(lin2.values[0])
    J_stair = float(lin1.values[0])
    delta = (J_smooth - J_stair) / abs(J_stair)
    record = dict(
        schema="radia.isochronous-topopt-shape-regen/v3",
        stl=regen.stl_info,
        reevaluation_mesh_family="hex",
        ne_staircase=int(m1.ne), ne_regen_hex=int(m2.ne),
        chi_eval=chi_eval,
        demag_factor_z_hex=demag_factor_z,
        J_smooth=J_smooth, J_staircase=J_stair, delta=delta,
        state_iterations=int(lin2.state_iterations),
        adjoint_iterations=[int(v) for v in lin2.adjoint_iterations],
        boundary_face_incident=(
            "a bare-Sculpt export wrote surfaceelements=0 and the demag "
            "solve was silently demag-free (J=+37.137, DemagFactor 0); "
            "fixed by mesh-based geometry before export, gated by "
            "boundary_faces_ok + the vim.ChargeGram ZERO-boundary raise + "
            "the demag_factor_z_hex lock here"),
        tet_solver_negative=(
            "16,102-tet native CG at chi=100 capped 20,000 iterations for "
            "mass-Riesz (tol 1e-10/1e-8/1e-6), Jacobi, and cluster-tree; "
            "tet export remains a gated mesh product, but production "
            "re-evaluation routes through the 1,537-hex Sculpt mesh"),
    )
    RESULTS.write_text(json.dumps(record, indent=1), encoding="utf-8")
    # broken-boundary exports read exactly 0 here; flat sector measured 0.477
    assert demag_factor_z > 0.2, demag_factor_z
    assert np.sign(J_smooth) == np.sign(J_stair)
    assert abs(delta) < 0.5, delta
