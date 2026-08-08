"""Hex-native design lane: the SLP loop runs end-to-end on Sculpt hexes.

The shape-regeneration lane covers density -> STL -> Sculpt -> re-evaluation
for a TET-designed density.  This lane locks the complementary half of the
Sculpt + HDiv-MMM pairing: the DESIGN loop itself on a Sculpt overlay-grid
hex mesh -- sector domain STL -> `cubit_stl_to_vol` hex -> `DensityAdjointVIM`
(RT1-hex charge Gram) -> `optimize_density` -> the hex path of
`iron_only_mesh` / `verify_design_iron_only`.  Needs a Cubit license
(headless batch per the driving policy); wall ~3-4 min on a development host.

Measured 2026-08-09 (LAB, Cubit 2025.12), the golden anchors:

* sector STL from the netgen surface export: watertight, V = 3.1415e-4
  (the analytic sector volume to 5 digits);
* Sculpt design mesh at size 0.012: 360 hexes, min quality 0.540, all
  `cubit_stl_to_vol` gates pass (closure 7.5 % at this coarse cell);
* design loop (chi=1000, 50 % volume budget, move 0.1): monotone,
  +3.3 % over 5 accepted iterates on the 2026-08-09 tree.  NOTE the
  accepted-iterate count and exact gain track the SLP driver
  (`topology_optimization.solve_lp_update`), which is under active
  development -- the tet reference lane measured 30 iterates / +16.1 %
  on an earlier tree state, and the SAME early-stop reproduces on tets
  there.  The gates below are therefore tree-robust (monotonicity and
  a wide gain band), NOT an iterate-count lock;
* hex iron-only extraction: 170 elements; matched-0/1 ersatz bands
  [+2.7 %, -4.7 %, -10.4 %] on the same tree state.
"""
import json
from math import cos, pi, sin
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia")
pytest.importorskip("trimesh")

from netgen.occ import Cylinder, HalfSpace, OCCGeometry, Pnt, Vec, Z  # noqa: E402
from ngsolve import HDiv, Mesh, SetNumThreads, TaskManager  # noqa: E402

from radia.isochronous_topopt import (  # noqa: E402
    MU0, DensityAdjointVIM, HelmholtzFilter, density_to_s,
    field_functional_load, gradient_pair_points, optimize_density,
    orbit_arc_points, uniform_field_load, verify_design_iron_only)

RESULTS = Path(__file__).with_name("results_hex_native_design_lane.json")


def _cubit_available():
    try:
        from radia_mcp.cubit import session as _cs  # noqa: F401
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(not _cubit_available(),
                                reason="radia_mcp cubit server not importable")


@pytest.fixture(scope="module")
def hexdesign(tmp_path_factory):
    from radia_mcp.cubit.server import cubit_stl_to_vol

    out = tmp_path_factory.mktemp("hex_native")
    SetNumThreads(4)
    ang = 60.0 * pi / 180.0
    ring = (Cylinder(Pnt(0, 0, 0.02), Z, r=0.15, h=0.03)
            - Cylinder(Pnt(0, 0, 0.01), Z, r=0.05, h=0.05))
    sector = ring * HalfSpace(Pnt(0, 0, 0), Vec(0, -1, 0)) \
        * HalfSpace(Pnt(0, 0, 0), Vec(-sin(ang), cos(ang), 0))
    with TaskManager():
        sec_mesh = Mesh(OCCGeometry(sector).GenerateMesh(maxh=0.01))
    sec_stl = out / "sector.stl"
    sec_mesh.ngmesh.Export(str(sec_stl), "STL Format")

    mesh_report = json.loads(cubit_stl_to_vol(
        stl_path=str(sec_stl), scheme="hex", size=0.012,
        closure_tolerance=0.15,
        out_vol=str(out / "design_hex.vol"),
        out_msh=str(out / "design_hex.msh")))

    chi_iron = 1000.0
    span = (pi / 12, pi / 2 - pi / 12)
    obj_pts, obj_radial = orbit_arc_points(0.115, 0.0, 7, span=span)
    pair_pts, pair_wts = gradient_pair_points(
        obj_pts, np.full(7, 1.0 / 7), delta=0.01, direction=obj_radial)

    def state_builder(fes):
        return uniform_field_load(fes, (0.0, 0.0, 1.0e5))

    def objective_builder(fes):
        return field_functional_load(fes, pair_pts, pair_wts, axis=2,
                                     scale=MU0, bonus_intorder=10)

    def constraint_builder(radius):
        def build(fes):
            cpts, _ = orbit_arc_points(radius, 0.0, 7, span=span)
            return field_functional_load(
                fes, cpts, np.full(7, 1.0 / 7), axis=2, scale=MU0,
                bonus_intorder=10)
        return build

    con_builders = [constraint_builder(r) for r in (0.08, 0.10)]
    with TaskManager():
        mesh = Mesh(str(out / "design_hex.vol"))
        fes = HDiv(mesh, order=1)
        prob = DensityAdjointVIM(fes, eps=1e-7)
        f_state = state_builder(fes)
        cons = [b(fes) for b in con_builders]
        filt = HelmholtzFilter(mesh, radius=0.012)
        lin0 = prob.linearize(
            density_to_s(filt.apply(np.full(prob.n_el, 0.5)), chi_iron),
            f_state, [objective_builder(fes)] + cons)
        targets = [float(v) for v in lin0.values[1:]]
        result = optimize_density(prob, f_state, objective_builder(fes),
                                  cons, targets, chi_iron=chi_iron,
                                  volume_fraction=0.5, density_filter=filt,
                                  move_limit=0.1, max_iterations=30)
        verification = verify_design_iron_only(
            prob, result.density, state_builder,
            [objective_builder] + con_builders, chi_iron=chi_iron,
            density_filter=filt, gram_kwargs=dict(eps=1e-7))
    return SimpleNamespace(mesh_report=mesh_report, mesh=mesh, fes=fes,
                           result=result, verification=verification)


def test_sculpt_design_mesh_gates(hexdesign):
    r = hexdesign.mesh_report
    assert r["status"] == "ok", r
    assert r["gates"]["boundary_faces_ok"] is True
    assert r["gates"]["no_inverted_elements"] is True
    assert r["gates"]["closure_ok"] is True
    by = r["by_type"][0]
    assert by["element"].startswith("Hexahedron")
    # measured 360 hexes / min quality 0.540 at size 0.012
    assert 150 <= by["n"] <= 900
    assert by["min"] > 0.3


def test_design_loop_converges_monotone_on_hexes(hexdesign):
    hist = hexdesign.result.history
    J = [float(h["objective"]) for h in hist]
    assert len(J) >= 2
    assert all(b >= a for a, b in zip(J, J[1:])), "accepted iterates regress"
    gain = (J[-1] - J[0]) / abs(J[0])
    # measured +3.3 % in 5 accepted iterates at this coarse resolution
    assert 0.005 < gain < 0.15, gain


def test_hex_iron_only_extraction_and_bands(hexdesign):
    v = hexdesign.verification
    assert v.iron_mesh.ne > 0
    first = next(iter(v.iron_mesh.Elements(ng.VOL)))
    assert len(first.vertices) == 8, "extraction must stay hex"
    bands = np.atleast_1d(np.asarray(v.bands, dtype=float))
    # measured [+2.7 %, -4.7 %, -10.4 %]; gate 20 %
    assert np.all(np.abs(bands) < 0.20), bands
    record = dict(
        schema="radia.isochronous-topopt-hex-native-design/v1",
        design_ne=int(hexdesign.mesh.ne),
        design_ndof=int(hexdesign.fes.ndof),
        design_min_quality=float(hexdesign.mesh_report["by_type"][0]["min"]),
        iterates=len(hexdesign.result.history),
        J_start=float(hexdesign.result.history[0]["objective"]),
        J_end=float(hexdesign.result.history[-1]["objective"]),
        iron_ne=int(v.iron_mesh.ne),
        ersatz_bands=[float(b) for b in bands],
    )
    RESULTS.write_text(json.dumps(record, indent=1), encoding="utf-8")
