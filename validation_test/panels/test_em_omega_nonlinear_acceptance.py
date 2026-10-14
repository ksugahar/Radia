"""Electromagnet block, Omega route: nonlinear solve and energy against an independent route.

``calc_accel_magnet.solve_accel`` (the Simulink Electromagnet block's Omega
method) solves a soft-iron block in a uniform applied field with the sample
steel B-H table.  The independent reference minimises the co-energy of the
same reduced-potential problem with NGSolve's line-search Newton, using the
shared production B(H) law and its antiderivative at integration points.

With order-1 H1 and a uniform source, H is constant in every element, so the
panel's element-wise secant model and the pointwise law define the same
discrete problem: the two solutions must agree to round-off, at an
unsaturated drive (0.05 T) and in the vacuum-slope tail (2 T).  The energy
is checked through the same antiderivative; in the iron, 0.5*B.H is not the
energy of a saturating law.

The A-Phi method (B-conforming) is checked against the same Omega reference:
Newton convergence, energy within 1 %, and B approaching the Omega value as
the mesh is refined.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
PANELS = REPO / "src" / "radia" / "panels"
SAMPLE_BH = PANELS / "samples" / "em_sample_bh.txt"
MU0 = 4e-7 * math.pi

rad = pytest.importorskip("radia")
ng = pytest.importorskip("ngsolve")
if str(PANELS) not in sys.path:
    sys.path.insert(0, str(PANELS))

SOURCE = """
import radia as rad


class UniformSource:
    current = 1.0

    def to_wire_segments(self, n_arc=20):
        return [], None

    def to_radia(self):
        return [rad.ObjBckg(lambda point: [0.0, 0.0, {b0}])]


def build_coil():
    return UniformSource()
"""


def _block_mesh(path, maxh_iron, maxh_air):
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, Sphere

    iron = Box(Pnt(-0.01, -0.01, -0.02), Pnt(0.01, 0.01, 0.02))
    iron.mat("yoke")
    iron.maxh = maxh_iron
    ball = Sphere(Pnt(0, 0, 0), 0.2)
    ball.faces.name = "outer"
    air = ball - iron
    air.mat("air")
    with ng.TaskManager():
        OCCGeometry(Glue([air, iron])).GenerateMesh(maxh=maxh_air).Save(str(path))
    return path


@pytest.fixture(scope="module")
def mesh_file(tmp_path_factory):
    return _block_mesh(tmp_path_factory.mktemp("omega") / "block.vol", 0.008, 0.06)


@pytest.fixture(scope="module")
def fine_mesh_file(tmp_path_factory):
    return _block_mesh(tmp_path_factory.mktemp("omega_fine") / "block.vol", 0.004, 0.04)


def _reference(vol, b0, table):
    from ngsolve.solvers import NewtonMinimization
    from radia.scalar_potential_solver import _build_bh_spline_law

    mesh = ng.Mesh(str(vol))
    b_of, coenergy_of, _ = _build_bh_spline_law(table)
    fes = ng.H1(mesh, order=1, dirichlet="outer")
    u, v = fes.TnT()
    H_s = ng.CF((0.0, 0.0, b0 / MU0))
    H = H_s - ng.grad(u)
    Hn = ng.sqrt(ng.InnerProduct(H, H) + 1e-20)
    density = mesh.MaterialCF({"yoke": coenergy_of(Hn), "air": 0.5 * MU0 * ng.InnerProduct(H, H)})
    a = ng.BilinearForm(fes, symmetric=True)
    a += ng.Variation(density * ng.dx(bonus_intorder=4))
    gfu = ng.GridFunction(fes)
    with ng.TaskManager():
        NewtonMinimization(
            a,
            gfu,
            freedofs=fes.FreeDofs(),
            maxit=200,
            maxerr=1e-12,
            inverse="sparsecholesky",
            linesearch=True,
            printing=False,
        )
    Hf = H_s - ng.grad(gfu)
    Hfn = ng.sqrt(ng.InnerProduct(Hf, Hf) + 1e-20)
    B = mesh.MaterialCF({"yoke": b_of(Hfn) / Hfn * Hf, "air": MU0 * Hf})
    coenergy = mesh.MaterialCF(
        {"yoke": coenergy_of(Hfn), "air": 0.5 * MU0 * ng.InnerProduct(Hf, Hf)}
    )
    W_bh = ng.Integrate(ng.InnerProduct(B, Hf), mesh, order=10)
    W_co = ng.Integrate(coenergy, mesh, order=10)
    iron = mesh.Materials("yoke")
    iron_bh = ng.Integrate(ng.InnerProduct(B, Hf), mesh, definedon=iron, order=10)
    iron_co = ng.Integrate(coenergy_of(Hfn), mesh, definedon=iron, order=10)

    # Residual of the reference's own equation, normalised by the same
    # nonlinear residual at the zero initial potential.
    free = np.array([bool(flag) for flag in fes.FreeDofs()])
    Hs_n = ng.sqrt(ng.InnerProduct(H_s, H_s) + 1e-20)
    B_initial = mesh.MaterialCF({"yoke": b_of(Hs_n) / Hs_n * H_s, "air": MU0 * H_s})
    norms = []
    for field in (B, B_initial):
        form = ng.LinearForm(fes)
        form += ng.InnerProduct(field, ng.grad(v)) * ng.dx(bonus_intorder=4)
        form.Assemble()
        norms.append(float(np.linalg.norm(form.vec.FV().NumPy()[free])))
    assert norms[1] > 0.0
    return {
        "B_origin": np.array(B(mesh(0.0, 0.0, 0.0))),
        "W": W_bh - W_co,
        "W_coenergy": W_co,
        "iron_energy": iron_bh - iron_co,
        "iron_half_BH": 0.5 * iron_bh,
        "residual": norms[0] / norms[1],
    }


@pytest.mark.parametrize("b0", [0.05, 2.0])
def test_omega_newton_panel_equals_the_independent_route(tmp_path, mesh_file, b0):
    from calc_accel_magnet import solve_accel
    from calc_common import EMMaterial

    source = tmp_path / "source.py"
    source.write_text(SOURCE.format(b0=b0), encoding="utf-8")
    material = EMMaterial.from_name("steel", bh_file=str(SAMPLE_BH))
    table, _ = material.get_bh_curve()
    panel = solve_accel(
        coil_script=str(source),
        vol_file=str(mesh_file),
        formulation="omega",
        fes_order=1,
        mat=material,
        max_iter=40,
        tol=1e-8,
        relax=1.0,
        newton=True,
        solver="sparsecholesky",
    )
    assert "error" not in panel, panel.get("error")
    assert panel["converged"] is True
    assert panel["nonlinear_relative_residual"] < 1e-10
    assert panel["iterations"] <= 12  # Newton, not a slow fixed point

    ref = _reference(mesh_file, b0, table)
    assert ref["residual"] < 1e-10
    B_panel = np.array(panel["B_origin"])
    assert np.linalg.norm(B_panel - ref["B_origin"]) <= 1e-9 * np.linalg.norm(ref["B_origin"])
    assert panel["W_mag"] == pytest.approx(ref["W"], rel=1e-9)
    assert panel["W_coenergy"] == pytest.approx(ref["W_coenergy"], rel=1e-9)
    assert panel["W_mag"] + panel["W_coenergy"] == pytest.approx(2.0 * panel["W_half_BH"], rel=1e-9)
    assert panel["energy_convention"].startswith("reversible nonlinear")
    if b0 == 2.0:
        # Saturated iron: 0.5*B.H is not the energy of the law.
        assert abs(ref["iron_half_BH"] - ref["iron_energy"]) > 0.05 * ref["iron_energy"]


@pytest.mark.parametrize("b0", [0.05, 2.0])
def test_aphi_newton_panel_converges_toward_the_dual_omega_route(
    tmp_path, mesh_file, fine_mesh_file, b0
):
    """A-Phi is B-conforming and Omega H-conforming: they are independent
    discretizations of the same problem and approach each other under
    refinement (the outer boundary condition differs: n x A = 0 vs phi = 0)."""
    from calc_accel_magnet import solve_accel
    from calc_common import EMMaterial

    source = tmp_path / "source.py"
    source.write_text(SOURCE.format(b0=b0), encoding="utf-8")
    material = EMMaterial.from_name("steel", bh_file=str(SAMPLE_BH))
    table, _ = material.get_bh_curve()
    differences = []
    for vol in (mesh_file, fine_mesh_file):
        panel = solve_accel(
            coil_script=str(source),
            vol_file=str(vol),
            formulation="a",
            fes_order=1,
            mat=material,
            max_iter=40,
            tol=1e-8,
            relax=1.0,
            newton=True,
            solver="sparsecholesky",
        )
        assert "error" not in panel, panel.get("error")
        assert panel["converged"] is True
        assert panel["nonlinear_relative_residual"] < 1e-10
        assert panel["iterations"] <= 10
        assert panel["energy_convention"].startswith("reversible nonlinear")
        ref = _reference(vol, b0, table)
        assert ref["residual"] < 1e-10
        # Energy of the two dual discretizations agrees to 1 % (measured
        # 0.10 % at 2 T and 0.65-0.70 % at 0.05 T on these meshes).
        assert panel["W_mag"] == pytest.approx(ref["W"], rel=1e-2)
        differences.append(
            np.linalg.norm(np.array(panel["B_origin"]) - ref["B_origin"])
            / np.linalg.norm(ref["B_origin"])
        )
    # B inside the iron: the two routes approach each other under refinement.
    assert differences[1] < differences[0]


def test_omega_picard_rejects_a_relaxation_that_never_updates(tmp_path, mesh_file):
    from calc_accel_magnet import solve_accel
    from calc_common import EMMaterial

    source = tmp_path / "source.py"
    source.write_text(SOURCE.format(b0=0.05), encoding="utf-8")
    material = EMMaterial.from_name("steel", bh_file=str(SAMPLE_BH))
    with pytest.raises(ValueError, match="relax"):
        solve_accel(
            coil_script=str(source),
            vol_file=str(mesh_file),
            formulation="omega",
            fes_order=1,
            mat=material,
            relax=0.0,
            solver="sparsecholesky",
        )
