#!/usr/bin/env python
"""Real-material hysteresis on the loop-free HDiv-VIM: the K=40 Potter-Schmulian play
model (validation_test/hysteresis/binput_play_fixture.npz, a literature analytical
model identified via radia.hysteresis_io) driven through the coupled demag solve.

Two results a K=3 toy cannot demonstrate and a 0D material model cannot produce:

  (1) FIDELITY -- on a low-demag rod the coupled (H_internal, B) major loop reproduces
      the demag-free intrinsic material loop (the demag only shears H_ext; the
      (H_int, B) trace IS the material law).  The intrinsic 0D loop is the anchor.
  (2) SPATIAL LOSS -- in a cube (B non-uniform from demag) the per-element loop area
      oint H.dB is a spatially-varying hysteresis-loss density; the loop-free solve
      resolves the spatial map + integrates the total, which a 0D loss curve cannot.

Correctness (numerical agreement / physical soundness), so LAB per the Benchmark
Policy; per-step timing uses hibino first or an idle-CI mdx fallback. Results saved to
real_material_*.json.
"""
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng                                    # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh         # noqa: E402

import radia.vim as vim                                 # noqa: E402

HERE = Path(__file__).resolve().parent
FIX = HERE / "binput_play_fixture.npz"
MU0 = 4.0e-7 * np.pi
RHO_STEEL = 7650.0


def _load_play():
    d = np.load(FIX, allow_pickle=True)
    K = int(d["K"]); eta = np.asarray(d["eta"], float)
    tables = [(np.asarray(d[f"r_{k}"], float), np.asarray(d[f"f_{k}"], float)) for k in range(K)]
    return K, eta, tables


def _loop_area(H, B):
    H = np.asarray(H, float); B = np.asarray(B, float)
    return 0.5 * float(abs(np.sum(H * np.roll(B, -1) - np.roll(H, -1) * B)))


def _rod(nz=16, a=0.1, aspect=10.0):
    L = aspect * a
    return MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=nz,
                                mapping=lambda x, y, z: (a * (x - 0.5), a * (y - 0.5), L * (z - 0.5)))


def _drive(H0, n_branch):
    up0 = np.linspace(0.0, H0, n_branch // 2 + 1)[1:]
    down = np.linspace(H0, -H0, n_branch + 1)[1:]
    up1 = np.linspace(-H0, H0, n_branch + 1)[1:]
    hz = np.concatenate([up0, down, up1])
    h = np.zeros((hz.size, 3)); h[:, 2] = hz
    return h, n_branch // 2


def _intrinsic_loop(mat, Bmax, nseg=300):
    st = mat.state0()[None, :].copy()
    Bseq = np.concatenate([np.linspace(0, Bmax, nseg), np.linspace(Bmax, -Bmax, 2 * nseg),
                           np.linspace(-Bmax, Bmax, 2 * nseg)])
    H = []
    for bz in Bseq:
        B = np.array([[0.0, 0.0, bz]]); H.append(mat.forward(B, st)[0, 2]); st = mat.commit(B, st)
    H = np.array(H)
    return _loop_area(H[nseg:], Bseq[nseg:])


def test_real_material_coupled_loop_matches_intrinsic():
    """(1) The coupled rod major loop reproduces the demag-free intrinsic loop."""
    K, eta, tables = _load_play()
    h_steps, n0 = _drive(H0=16.0e3, n_branch=48)     # in-range (~0.8 T, under b_max=1.9 T)
    with ng.TaskManager():
        mesh = _rod()
        res = vim.SolveHysteresis(mesh, h_steps, material=vim.PlayHysteresisMaterial(K, eta, tables))

    Hz = np.array([s["H_avg"][2] for s in res["steps"]])
    Bz = np.array([s["B_avg"][2] for s in res["steps"]])
    iters = [int(s["iters"]) for s in res["steps"]]
    peakB = float(np.abs(Bz).max())
    area = _loop_area(Hz[n0:], Bz[n0:])
    area0 = _intrinsic_loop(vim.PlayHysteresisMaterial(K, eta, tables), peakB)
    rel = abs(area - area0) / area0

    assert 0.5 < peakB < 1.9, "rod drive out of the intended in-range band (peak %.3f T)" % peakB
    assert max(iters) < 40, "polarization iterations blew up on the real material (max %d)" % max(iters)
    assert rel < 5.0e-2, (
        "coupled real-material loop disagrees with the intrinsic material loop: "
        "coupled %.1f vs intrinsic %.1f J/m^3 (rel %.2e)" % (area, area0, rel))

    (HERE / "real_material_loop.json").write_text(json.dumps(dict(
        description="Real K=40 Potter-Schmulian coupled rod loop vs demag-free intrinsic loop",
        timestamp=datetime.now().isoformat(), hostname=platform.node(),
        peakB_T=peakB, coupled_area_J_per_m3=area, intrinsic_area_J_per_m3=area0,
        loop_area_rel_diff=rel, loss_coupled_W_per_kg_50Hz=area * 50.0 / RHO_STEEL,
        picard_iters=iters, ndof=res["ndof"], n_el=res["n_el"],
        H_avg_z=[float(v) for v in Hz], B_avg_z=[float(v) for v in Bz]), indent=1))


def test_real_material_spatial_loss_nonuniform():
    """(2) The per-element loop area gives a spatially-varying loss density in the cube."""
    K, eta, tables = _load_play()
    h_steps, n0 = _drive(H0=2.0e5, n_branch=16)      # cube demag ~0.27 -> peak ~1.35 T (in range)
    N = 5
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=N, ny=N, nz=N,
                                    mapping=lambda x, y, z: (x - 0.5, y - 0.5, z - 0.5))
        res = vim.SolveHysteresis(mesh, h_steps, material=vim.PlayHysteresisMaterial(K, eta, tables))
        vol_el = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True), float)

    steps = res["steps"][n0:]
    Hseq = np.stack([np.asarray(s["H"], float) for s in steps])
    Bseq = np.stack([np.asarray(s["B"], float) for s in steps])
    dB = np.roll(Bseq, -1, axis=0) - Bseq
    Hmid = 0.5 * (Hseq + np.roll(Hseq, -1, axis=0))
    loss_el = np.abs(np.sum(np.sum(Hmid * dB, axis=2), axis=0))       # J/m^3/cycle per element
    peakB = float(np.max(np.linalg.norm(Bseq.reshape(-1, 3), axis=1)))
    spread = float(loss_el.max() / max(loss_el.min(), 1e-9))
    total = float(np.sum(loss_el * vol_el))
    mean_dens = total / float(np.sum(vol_el))

    assert peakB < 1.9, "cube drive exceeded b_max (peak %.3f T)" % peakB
    assert loss_el.min() > 0.0, "an element recorded non-positive hysteresis loss (%.3e)" % loss_el.min()
    assert spread > 2.0, (
        "the loss density is nearly uniform (spread x%.2f) -- the spatial map should show "
        "the demag-driven B non-uniformity" % spread)

    (HERE / "real_material_spatial_loss.json").write_text(json.dumps(dict(
        description="Real K=40 Potter-Schmulian per-element hysteresis loss density, cube demag",
        timestamp=datetime.now().isoformat(), hostname=platform.node(),
        peakB_T=peakB, n_el=int(loss_el.size),
        loss_density_min=float(loss_el.min()), loss_density_max=float(loss_el.max()),
        loss_density_mean=mean_dens, spatial_spread=spread, total_loss_J_per_cycle=total), indent=1))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-s"]))
