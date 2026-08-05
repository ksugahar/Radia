"""Real-solve physics golden for the IH geometry->operator assembler.

test_ih_operator_assembly locks the assembly plumbing with a FAKED
unit-current solve; this file runs the REAL electromagnetic solve and
locks the physics.  The geometry is fully authored in-test (no repo
fixtures, no PEEC cache): a 355-degree gapped rectangular-wire torus
coil -- its centerline extraction takes the ANALYTIC revolution-sweep
path, so the whole module runs in seconds -- over a small labeled
workpiece cube solved with the real BEM-SIBC workflow.

* P_wp(1 A) sits in an absolute golden band AND matches a direct
  calc_inductance run whose argv is written out independently here --
  the assembler must not mangle a single physical setting on the way
  through (frequency, conductivities, mu_r, coupling, backend, PEEC
  filament settings).
* Power closure: dot(heat_cell_weights, heat_projection) equals the
  electromagnetic P_wp -- the thermal source preserves EM power.
* Thermal operator identities on the exported CSR operators:
  sum(M) = rho*cp*V, K@1 = 0 (pure Neumann Laplacian), sum(C) = total
  surface area (h is applied at runtime), temperature_cell_weights =
  row sums of M, against the EXACT box volume/area (1 cm cube).
* The exact-response contract: n_eddy_unknown = 1, A = [1], rhs = [1],
  so the runtime's heat scales as I^2 with no surrogate in between.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

# Exact geometry of the 1 cm workpiece cube (mesh-independent).
BOX_VOLUME_M3 = 1.0e-6
BOX_AREA_M2 = 6.0e-4
RHO_CP = 7800.0 * 467.0

# Coil: gapped torus, spine radius 30 mm, 6 x 4 mm rectangular wire,
# 355-degree sweep, centered on the workpiece cube's corner region.
COIL_R_M = 0.030
COIL_W_M = 0.006
COIL_H_M = 0.004
COIL_SWEEP_DEG = 355

# Golden band recorded for a netgen maxh=6 mm cube with 21 H1 DOFs and
# the torus coil above at 7 kHz, sigma_wp=5e6, mu_r=100, weak coupling,
# and intree-dense BEM: P_wp(1 A) = 5.7718475277e-06 W.  The +-2% band
# absorbs mesher-version drift of the tiny meshes; the direct-run
# comparison below is the tight argument-routing lock.
P_WP_GOLDEN_W = 5.7718475277e-06
P_WP_BAND_W = (0.98 * P_WP_GOLDEN_W, 1.02 * P_WP_GOLDEN_W)


@pytest.fixture(scope="module")
def assembled(tmp_path_factory):
    pytest.importorskip("ngsolve")
    pytest.importorskip("cubit_mesh_export")
    from netgen.occ import (
        Axes,
        Axis,
        Box,
        Dir,
        OCCGeometry,
        Pnt,
        WorkPlane,
    )

    from radia.simulink import ih_operator_assembly as assembly

    tmp = tmp_path_factory.mktemp("ih_physics")
    shape = Box(Pnt(0, 0, 0), Pnt(0.01, 0.01, 0.01))
    shape.mat("workpiece")
    for face in shape.faces:
        face.name = "sibc"
    workpiece = tmp / "workpiece.vol"
    OCCGeometry(shape).GenerateMesh(maxh=0.006).Save(str(workpiece))

    profile = (
        WorkPlane(
            Axes(
                p=Pnt(COIL_R_M, 0, 0),
                n=Dir(0, 1, 0),
                h=Dir(0, 0, 1),
            )
        )
        .RectangleC(COIL_W_M, COIL_H_M)
        .Face()
    )
    coil_solid = profile.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), COIL_SWEEP_DEG)
    coil = tmp / "gapped_torus_coil.step"
    coil_solid.WriteStep(str(coil))

    config = assembly.assemble_ih_operators(
        workpiece,
        coil,
        output=tmp / "native.json",
        run_dir=tmp / "run",
        options=assembly.IHOperatorAssemblyOptions(),
    )
    return {"config": config, "workpiece": workpiece, "coil": coil, "tmp": tmp}


def _dense(n: int, row_ptr, col, val) -> np.ndarray:
    matrix = np.zeros((n, n))
    for row in range(n):
        for k in range(row_ptr[row], row_ptr[row + 1]):
            matrix[row, col[k]] += val[k]
    return matrix


def test_p_wp_in_golden_band_and_matches_direct_calc_run(assembled):
    from radia.panels.calc_inductance import build_argparser, run_inductance

    power = assembled["config"]["unit_current"]["electromagnetic_power_W"]
    assert (
        P_WP_BAND_W[0] <= power <= P_WP_BAND_W[1]
    ), f"P_wp(1A) = {power:.6e} W left the golden band {P_WP_BAND_W}"

    # Independently spelled-out argv for the same physics; a flag the
    # assembler forgets or mangles shows up as a power mismatch here.
    tmp = assembled["tmp"]
    argv = [
        "--coil-solver",
        "peec",
        "--vol",
        str(assembled["workpiece"]),
        "--wp-label",
        "sibc",
        "--sigma",
        "5000000",
        "--mu-r",
        "100",
        "--frequency",
        "7000",
        "--current",
        "1",
        "--coil-sigma",
        "58000000",
        "--coupling-mode",
        "weak",
        "--wp-bem-backend",
        "intree-dense",
        "--h1-order",
        "1",
        "--wp-loop-dof",
        "auto",
        "--impedance-model",
        "sibc",
        "--msh-output",
        str(tmp / "direct.msh"),
        "--coil-step",
        str(assembled["coil"]),
        "--peec-n-peri",
        "16",
        "--peec-proximity",
    ]
    payload = run_inductance(build_argparser().parse_args(argv))
    assert payload.get("status") != "error" and not payload.get("error")
    direct_power = float(payload["P_wp_W"])
    assert power == pytest.approx(direct_power, rel=1.0e-9), (
        "the assembler and a direct calc_inductance run with identical "
        f"settings disagree: {power:.10e} vs {direct_power:.10e} W"
    )


def test_thermal_source_preserves_electromagnetic_power(assembled):
    config = assembled["config"]
    unit = config["unit_current"]
    assert unit["reference_current_A"] == 1.0
    assert unit["relative_power_error"] < 1.0e-6
    closure = float(np.dot(config["heat_cell_weights"], config["heat_projection"]))
    assert closure == pytest.approx(
        unit["electromagnetic_power_W"], rel=1.0e-9
    ), "dot(w, q) no longer reproduces the electromagnetic power"
    assert unit["thermal_source_power_W"] == pytest.approx(closure, rel=1.0e-12)


def test_thermal_operator_identities(assembled):
    config = assembled["config"]
    n = config["n_temperature"]
    mass = _dense(n, config["mass_row_ptr"], config["mass_col"], config["mass_value"])
    stiffness = _dense(
        n, config["stiffness_row_ptr"], config["stiffness_col"], config["stiffness_value"]
    )
    convection = _dense(
        n, config["convection_row_ptr"], config["convection_col"], config["convection_value"]
    )

    # P1 quadrature integrates constants exactly, so these hold to
    # machine precision against the EXACT cube volume / area.
    assert mass.sum() == pytest.approx(RHO_CP * BOX_VOLUME_M3, rel=1.0e-10)
    ones = np.ones(n)
    assert (
        np.abs(stiffness @ ones).max() <= 1.0e-10 * np.abs(stiffness).max()
    ), "the conduction operator lost its constant nullspace (K@1 != 0)"
    assert convection.sum() == pytest.approx(BOX_AREA_M2, rel=1.0e-10), (
        "the raw convection operator must integrate uv over the whole "
        "surface (h is applied at runtime)"
    )
    weights = np.asarray(config["temperature_cell_weights"], dtype=float)
    assert weights == pytest.approx(
        mass.sum(axis=1), rel=1.0e-10
    ), "temperature_cell_weights must be the lumped (row-sum) mass"
    assert np.all(weights > 0.0)
    assert math.isclose(weights.sum(), RHO_CP * BOX_VOLUME_M3, rel_tol=1.0e-10)


def test_exact_single_current_response_contract(assembled):
    config = assembled["config"]
    assert config["operator_basis"] == "exact-single-current-linear-response"
    assert config["surrogate"] is False
    assert config["n_eddy_unknown"] == 1
    assert config["eddy_matrix_real"] == [1.0]
    assert config["eddy_matrix_imag"] == [0.0]
    assert config["eddy_rhs_real"] == [1.0]
    assert config["eddy_rhs_imag"] == [0.0]
    assert config["current_change_recomputes_eddy"] is False
    assert config["n_heat"] == len(config["heat_projection"])
    assert config["n_heat"] == len(config["heat_cell_weights"])
    assert all(value >= 0.0 for value in config["heat_projection"])
    assert all(value > 0.0 for value in config["heat_cell_weights"])
