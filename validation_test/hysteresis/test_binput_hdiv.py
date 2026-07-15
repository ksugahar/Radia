#!/usr/bin/env python
"""Golden validation for B-input PLAY-model hysteresis on the HDiv-VIM (RT1) charge Gram.

HDiv successor of the retired moment-path golden (test_binput_moment.py, removed together
with the collocation MMMM in 75a0086f).  Same synthetic play model, adapted physics gates:

  (A) VIRGIN-CURVE CROSS-CHECK (the trusted-reference gate): on a monotone first-
      magnetization drive the play model IS a single-valued B(H) curve, so
      vim.SolveHysteresis (one step from the virgin state) must reproduce
      vim.Solve(bh_table=<virgin curve sampled from the SAME model>) -- the production
      anhysteretic energy-Newton path -- in the volume-averaged M_z.  Runs on hex AND tet.
  (B) HYSTERETIC LOOP: a cyclic H_z drive yields a (H_internal, B) loop with nonzero
      enclosed area (dissipation, not single-valued); committed play states carry the
      branch history across steps.  Loop data saved to binput_hdiv_loop.json.
  (C) covered by parametrizing (A) over hex/tet + the loop run.
  (E) STRONG-COUPLING reversal robustness: the 4x4x4 block full loop -- the regime
      where the retired moment path's plain Picard diverged -- runs LIVE and persists
      binput_hdiv_loop_4x4x4.json.

SolveHysteresis runs the Hantila polarization iteration (constant SPD LHS W(nu0)+N built
once, lagged VECTOR polarization source nu0*M - H_mat(M)), which represents the recoil-
branch ANTI-PARALLEL H/M states -- the structural regime where a scalar secant
nu = |H|/|M| >= 0 cannot converge (the descending-branch divergence the retired moment
path hit with plain Picard).
"""

import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

import ngsolve as ng                      # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

import radia.vim as vim                   # noqa: E402  (imports radia's C++ core)

MU_0 = 4.0 * np.pi * 1e-7
HERE = Path(__file__).resolve().parent

H0 = 200.0e3          # peak applied field (A/m): drives the cube to B ~ 0.75 T (crosses eta_2)


# --------------------------------------------------------------------------
# Synthetic play model (verbatim from the retired moment golden): K=3,
# eta=[0, 0.3, 0.7] T, piecewise-linear shape functions.  f_0 (reversible,
# eta=0) sets the small-signal mu_r; f_1, f_2 (negative, the Hane-Sugahara
# Potter sign convention) add the irreversible / hysteretic part.
# --------------------------------------------------------------------------
def _synthetic_play(mu_r=200.0, irr1=0.30, irr2=0.15, eta=(0.0, 0.3, 0.7),
                    rmax=2.0, n=41):
    r = np.linspace(0.0, rmax, n)
    a0 = 1.0 / (MU_0 * mu_r)            # reversible slope dH/dB ~ a0 -> mu_r
    f0 = a0 * r
    f1 = -irr1 * a0 * r
    f2 = -irr2 * a0 * r
    eta = np.asarray(eta, dtype=float)
    tables = [(r.tolist(), f0.tolist()),
              (r.tolist(), f1.tolist()),
              (r.tolist(), f2.tolist())]
    return len(eta), eta, tables


def _hex_cube(n=3):
    return MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n,
                                mapping=lambda x, y, z: (x - 0.5, y - 0.5, z - 0.5))


def _tet_cube(maxh=0.45):
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def _virgin_bh_table(material, Bmax=1.6, n=80):
    """Sample the model's first-magnetization curve H(B) from the virgin state.

    Monotone from virgin => single-valued => a legitimate [[H, B], ...] table for the
    anhysteretic reference path.  Every sample restores the SAME virgin state, so this
    is exactly the curve the coupled solve sees on a monotone first drive.
    """
    s0 = material.state0()
    rows = [[0.0, 0.0]]
    Bs = np.linspace(0.0, Bmax, n + 1)[1:]
    Bq = np.zeros((Bs.size, 3)); Bq[:, 2] = Bs
    Hq = material.forward(Bq, np.tile(s0[None, :], (Bs.size, 1)))
    for b, h in zip(Bs, Hq[:, 2]):
        rows.append([float(h), float(b)])
    H_col = np.array([r[0] for r in rows])
    assert np.all(np.diff(H_col) > 0.0), "virgin curve must be strictly H-increasing"
    return rows


def _loop_area(h, b):
    """Shoelace area of the closed (h, b) polygon."""
    h = np.asarray(h, float); b = np.asarray(b, float)
    return 0.5 * float(np.sum(h * np.roll(b, -1) - np.roll(h, -1) * b))


# ==========================================================================
# (A) Virgin-curve cross-check vs the production anhysteretic energy-Newton.
# ==========================================================================
@pytest.mark.parametrize("mesh_builder", [_hex_cube, _tet_cube], ids=["hex", "tet"])
def test_A_virgin_ramp_matches_anhysteretic_reference(mesh_builder):
    K, eta, tables = _synthetic_play()
    with ng.TaskManager():
        mesh = mesh_builder()
        mat = vim.PlayHysteresisMaterial(K, eta, tables)
        res = vim.SolveHysteresis(mesh, [[0.0, 0.0, H0]], material=mat)
        ref = vim.Solve(mesh, bh_table=_virgin_bh_table(mat),
                        H_ext=ng.CoefficientFunction((0.0, 0.0, H0)))
    Mz = float(res["steps"][0]["M_avg"][2])
    Mz_ref = float(ref["M_avg"][2])
    rel = abs(Mz - Mz_ref) / abs(Mz_ref)
    assert Mz_ref > 1.0e5, "reference drive too weak to be meaningful (Mz_ref=%.3e)" % Mz_ref
    assert rel < 5.0e-3, (
        "virgin-ramp SolveHysteresis disagrees with the anhysteretic reference: "
        "Mz=%.6e vs ref %.6e (rel %.2e)" % (Mz, Mz_ref, rel))


# ==========================================================================
# (B) Cyclic drive -> hysteretic (H_internal, B) loop with nonzero area.
# ==========================================================================
def test_B_cyclic_loop_hysteretic():
    K, eta, tables = _synthetic_play()
    up0 = np.linspace(0.0, H0, 5)[1:]              # virgin ramp     (4 steps)
    down = np.linspace(H0, -H0, 9)[1:]             # descending      (8 steps, includes 0)
    up1 = np.linspace(-H0, H0, 9)[1:]              # ascending       (8 steps, includes 0)
    hz = np.concatenate([up0, down, up1])
    h_steps = np.zeros((hz.size, 3)); h_steps[:, 2] = hz
    with ng.TaskManager():
        mesh = _hex_cube(3)
        res = vim.SolveHysteresis(mesh, h_steps, play=(K, eta, tables))

    assert res["permanent_magnet_model"] == "simplified-play"
    assert res["permanent_magnet_level"] == 3

    Hz_int = np.array([s["H_avg"][2] for s in res["steps"]])
    Bz = np.array([s["B_avg"][2] for s in res["steps"]])

    # closed cycle = the last 16 steps (+H0 -> -H0 -> +H0)
    area = abs(_loop_area(Hz_int[4:], Bz[4:]))
    b_span = float(Bz.max() - Bz.min())
    assert b_span > 0.2, "cyclic drive barely magnetizes the cube (span %.3e T)" % b_span
    assert area > 50.0, (
        "the (H_internal, B) cycle encloses no area (area=%.3e A*T/m) -- the play states "
        "are not carrying hysteresis across steps" % area)

    # branch separation: descending vs ascending B at the same internal field near H=0
    i_desc = 4 + int(np.argmin(np.abs(Hz_int[4:12])))
    i_asc = 12 + int(np.argmin(np.abs(Hz_int[12:])))
    assert Bz[i_desc] - Bz[i_asc] > 1.0e-3, (
        "descending/ascending branches coincide near H_int=0 (dB=%.3e T)"
        % (Bz[i_desc] - Bz[i_asc]))

    payload = dict(
        description="HDiv-VIM B-input play hysteresis cube loop (gate B)",
        timestamp=datetime.now().isoformat(),
        hostname=platform.node(),
        h_applied_z=[float(v) for v in hz],
        H_internal_avg_z=[float(v) for v in Hz_int],
        B_avg_z=[float(v) for v in Bz],
        M_avg_z=[float(s["M_avg"][2]) for s in res["steps"]],
        loop_area_A_T_per_m=float(area),
        ndof=res["ndof"], n_el=res["n_el"],
        charge_gram_wall_s=res["charge_gram_wall_s"],
        t_steps_s=res["t_steps_s"],
        picard_iters=[int(s["iters"]) for s in res["steps"]],
    )
    (HERE / "binput_hdiv_loop.json").write_text(json.dumps(payload, indent=1))


# ==========================================================================
# (E) Strong-coupling regression: full loop on a 4x4x4 block.
# ==========================================================================
def test_E_strong_coupling_block_loop():
    """The retired moment path's plain Picard DIVERGED on the descending branch of a
    4x4x4 strongly-coupled block (it needed the Anderson safeguard).  The Hantila
    polarization iteration must complete the same full loop, stay hysteretic, and
    keep the per-step iteration count bounded."""
    K, eta, tables = _synthetic_play()
    up0 = np.linspace(0.0, H0, 5)[1:]
    down = np.linspace(H0, -H0, 9)[1:]
    up1 = np.linspace(-H0, H0, 9)[1:]
    hz = np.concatenate([up0, down, up1])
    h_steps = np.zeros((hz.size, 3)); h_steps[:, 2] = hz
    with ng.TaskManager():
        mesh = _hex_cube(4)
        res = vim.SolveHysteresis(mesh, h_steps, play=(K, eta, tables))

    Hz_int = np.array([s["H_avg"][2] for s in res["steps"]])
    Bz = np.array([s["B_avg"][2] for s in res["steps"]])
    area = abs(_loop_area(Hz_int[4:], Bz[4:]))
    it_max = max(s["iters"] for s in res["steps"])
    assert area > 50.0, "4x4x4 block loop lost its hysteresis (area=%.3e)" % area
    assert it_max < 150, "polarization iteration count blew up on a reversal (max %d)" % it_max

    payload = dict(
        description="HDiv-VIM B-input play hysteresis 4x4x4 strong-coupling loop (gate E)",
        timestamp=datetime.now().isoformat(),
        hostname=platform.node(),
        h_applied_z=[float(v) for v in hz],
        H_internal_avg_z=[float(v) for v in Hz_int],
        B_avg_z=[float(v) for v in Bz],
        M_avg_z=[float(s["M_avg"][2]) for s in res["steps"]],
        loop_area_A_T_per_m=float(area),
        ndof=res["ndof"], n_el=res["n_el"],
        charge_gram_wall_s=res["charge_gram_wall_s"],
        t_steps_s=res["t_steps_s"],
        picard_iters=[int(s["iters"]) for s in res["steps"]],
    )
    (HERE / "binput_hdiv_loop_4x4x4.json").write_text(json.dumps(payload, indent=1))


# ==========================================================================
# (F) INDEPENDENT anchors via the sub-threshold linear limit.
#
# Gate A's reference table is sampled from the SAME material object, so a
# material-wiring bug would cancel on both sides.  Below the first
# irreversible threshold (|B| < eta_1 = 0.3 T) the synthetic play model is
# EXACTLY the linear material mu_r = 200 (only the eta=0 reversible operator
# is active), which pins two independent checks:
#   (F1) a single sub-threshold step must match the plain LINEAR demag solve
#        vim.Solve(mesh, mu_r, H_ext) -- a chi-based path that never touches
#        the MatHys* machinery (and is itself locked against analytic demag
#        factors elsewhere);
#   (F2) a full sub-threshold cycle must enclose ~ZERO area (no hysteresis
#        below the threshold), where the full-drive loop encloses ~1173.
# ==========================================================================
H_SMALL = 40.0e3      # keeps every element's |B| below eta_1 (asserted below)


def test_F_subthreshold_linear_limit():
    K, eta, tables = _synthetic_play()
    with ng.TaskManager():
        mesh = _hex_cube(3)
        res = vim.SolveHysteresis(mesh, [[0.0, 0.0, H_SMALL]], play=(K, eta, tables))
        lin = vim.Solve(mesh, 200.0, ng.CoefficientFunction((0.0, 0.0, H_SMALL)), order=1)

    B_el = np.asarray(res["steps"][0]["B"], float)
    assert float(np.max(np.linalg.norm(B_el, axis=1))) < 0.30, (
        "precondition broken: an element exceeded the first play threshold; lower H_SMALL")
    Mz = float(res["steps"][0]["M_avg"][2])
    Mz_lin = float(lin["M_avg"][2])
    rel = abs(Mz - Mz_lin) / abs(Mz_lin)
    assert Mz_lin > 1.0e4, "linear reference too weak to be meaningful (Mz=%.3e)" % Mz_lin
    assert rel < 2.0e-3, (
        "sub-threshold SolveHysteresis disagrees with the INDEPENDENT linear demag solve: "
        "Mz=%.6e vs mu_r=200 reference %.6e (rel %.2e)" % (Mz, Mz_lin, rel))


def test_F_subthreshold_cycle_has_no_hysteresis():
    K, eta, tables = _synthetic_play()
    up0 = np.linspace(0.0, H_SMALL, 4)[1:]
    down = np.linspace(H_SMALL, -H_SMALL, 7)[1:]
    up1 = np.linspace(-H_SMALL, H_SMALL, 7)[1:]
    hz = np.concatenate([up0, down, up1])
    h_steps = np.zeros((hz.size, 3)); h_steps[:, 2] = hz
    with ng.TaskManager():
        mesh = _hex_cube(3)
        res = vim.SolveHysteresis(mesh, h_steps, play=(K, eta, tables))

    B_all = np.concatenate([np.asarray(s["B"], float) for s in res["steps"]])
    assert float(np.max(np.linalg.norm(B_all, axis=1))) < 0.30, (
        "precondition broken: the cycle crossed the first play threshold; lower H_SMALL")
    Hz_int = np.array([s["H_avg"][2] for s in res["steps"]])
    Bz = np.array([s["B_avg"][2] for s in res["steps"]])
    area = abs(_loop_area(Hz_int[3:], Bz[3:]))
    assert area < 2.0, (
        "a sub-threshold (purely reversible) cycle enclosed area %.3e J/m^3 -- the play "
        "history advanced where the material is exactly linear" % area)


# ==========================================================================
# (G) Out-of-range fail-loud guard: driving past the material's identified
# table extent (b_max) must RAISE, not silently extrapolate the shape
# functions (a demag-limited coupled solve otherwise runs M away -- a real
# incident: B_avg = 4.886 T on a low-demag rod driven too hard).
# ==========================================================================
def test_G_out_of_range_drive_fails_loud():
    K, eta, tables = _synthetic_play()
    mat = vim.PlayHysteresisMaterial(K, eta, tables)
    b_max = mat.b_max()
    assert b_max > 1.5, "synthetic play table extent unexpectedly small (%.3f T)" % b_max
    # A huge single step: the cube demag caps B, but mu_r=200 with no saturation
    # is driven well past the 2 T table extent -> the guard must fire.
    with ng.TaskManager():
        mesh = _hex_cube(3)
        with pytest.raises(RuntimeError, match="identified range"):
            vim.SolveHysteresis(mesh, [[0.0, 0.0, 8.0e5]], material=mat)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-x"]))
