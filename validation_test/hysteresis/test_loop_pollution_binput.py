#!/usr/bin/env python
"""Controlled loop-injection experiment: loop content reaching the B-input
constitutive law corrupts the hysteresis solve (the collocation failure mode,
reproduced in a controlled way on the loop-free HDiv solver).

Mechanism under test: a discrete loop (ker of the charge map B) produces no
field, but it DOES carry magnetization |M|.  A B-input material law reads
B = mu0*(H + M), so any loop content that reaches the law shifts the material
operating point toward saturation.  The loop-free HDiv-Galerkin solve feeds
the law loop-free magnetization by construction; a non-symmetric collocation
moment operator amplifies loop content (loops are its near-null space, scaled
by chi on inversion), and the loop |M| fraction measured on such operators
was 5-40% at mu_r 1e3-1e4.

Harness: a faithful mirror of vim.SolveHysteresis's Hantila iteration with one
lever -- before every material evaluation the element-averaged magnetization
is polluted by an exact-loop field (seeded z0 in ker(B), element-averaged):

    M_fed = M_el + eps * (rms(M_el)/rms(L_el)) * L_el

so eps IS the relative loop fraction fed to the law.  The eps = 0 leg must
match the production vim.SolveHysteresis result (fidelity gate) -- this also
self-polices harness drift: if the production iteration changes, eps = 0
fails first.  Material: the committed real K=40 Potter-Schmulian play fixture.

Locked findings (3^3 cube, H0 = 200 kA/m, 20-step loop, measured 2026-07-13):
  eps = 0.05 -> peak|B| +19% spurious saturation, ~2.0x polarization iters
  eps = 0.10 -> peak|B| +52%, ~2.4x iters
  eps = 0.20 -> RUNAWAY past the material's identified range (b_max guard)
  eps = 0.30 -> polarization iteration diverges outright
This is the quantitative backing for the journal claim that complete loop
freedom is a hysteresis-correctness requirement, not solver hygiene.
Runtime ~3 min (validation_test lane).  Results: loop_pollution_binput.json.
"""
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import scipy.sparse as sp                               # noqa: E402
import ngsolve as ng                                    # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh         # noqa: E402

import radia.vim as vim                                 # noqa: E402
from radia.vim._hysteresis import (PlayHysteresisMaterial,      # noqa: E402
                                   _solve_pointwise_B, MU0)
from radia.vim._vim import build_charge_gram            # noqa: E402
from radia.vim._solve import (  # noqa: E402
    _f64,
    _h_solve_mass_riesz,
    _configure_cpp_mass,
)

HERE = Path(__file__).resolve().parent
FIX = HERE / "binput_play_fixture.npz"
NX = 3
H0 = 2.0e5
N_RAMP, N_BRANCH = 4, 8
EPS_LIST = [0.0, 0.05, 0.1, 0.2, 0.3]
SEED = 20260713


def _load_play():
    d = np.load(FIX, allow_pickle=True)
    K = int(d["K"]); eta = np.asarray(d["eta"], float)
    tables = [(np.asarray(d[f"r_{k}"], float), np.asarray(d[f"f_{k}"], float)) for k in range(K)]
    return K, eta, tables


def _drive():
    up0 = np.linspace(0.0, H0, N_RAMP + 1)[1:]
    down = np.linspace(H0, -H0, N_BRANCH + 1)[1:]
    up1 = np.linspace(-H0, H0, N_BRANCH + 1)[1:]
    hz = np.concatenate([up0, down, up1])
    h = np.zeros((hz.size, 3)); h[:, 2] = hz
    return h


def _loop_area(H, B):
    return 0.5 * float(abs(np.sum(H * np.roll(B, -1) - np.roll(H, -1) * B)))


def _run_harness(mesh, h_steps, material, eps_loop,
                 nl_tol=1e-3, nl_maxit=200, tol=1e-8, maxit=4000):
    x = MU0 * float(material.nu_bound())
    nu0 = x / (1.0 - x)
    b_max = float(material.b_max())

    fes = ng.HDiv(mesh, order=1)
    B, H, M_mass = build_charge_gram(fes, eps=1e-10, leafsize=32, eta=2.0,
                                     far_quad=None, ho_far_factor=None, nonlinear=True)
    Mm = sp.csr_matrix(M_mass)
    Bc = sp.csr_matrix(B)
    n_face = fes.ndof
    n_el = mesh.GetNE(ng.VOL)
    _configure_cpp_mass(H, Mm, int(n_face))

    uf = fes.TrialFunction()
    vf = fes.TestFunction()
    l2 = ng.L2(mesh, order=0)
    wl2 = l2.TestFunction()
    vol_el = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True), float)
    w_el = vol_el / float(np.sum(vol_el))
    blocks = []
    for c in range(3):
        blf = ng.BilinearForm(trialspace=fes, testspace=l2)
        blf += uf[c] * wl2 * ng.dx(bonus_intorder=4)
        blf.Assemble()
        r, ci, v = blf.mat.COO()
        blocks.append(sp.csr_matrix((_f64(v), (np.asarray(r), np.asarray(ci))),
                                    shape=(l2.ndof, n_face)))
    P = sp.vstack(blocks).tocsr()
    PT = P.T.tocsr()

    # exact loop field (seeded direction in ker(B)), element-averaged
    Bd = np.asarray(Bc.todense())
    s_ = np.linalg.svd(Bd, compute_uv=False)
    rank = int(np.sum(s_ > s_.max() * max(Bd.shape) * np.finfo(float).eps))
    _, _, vt = np.linalg.svd(Bd, full_matrices=True)
    Z = vt[rank:].T
    rng = np.random.default_rng(SEED)
    z0 = Z @ rng.standard_normal(Z.shape[1])
    L_el = (P @ z0).reshape(3, n_el).T / vol_el[:, None]
    L_rms = float(np.sqrt(np.mean(np.sum(L_el ** 2, axis=1))))

    def solve_W0(rhs, x0=None):
        res = _h_solve_mass_riesz(
            H, None, int(n_face), float(nu0), rhs, tol, int(maxit), x0=x0)
        return np.asarray(res["m"], float)

    states = np.tile(material.state0()[None, :], (n_el, 1))
    B_cache = None
    s_el = np.zeros((n_el, 3))
    m = np.zeros(n_face)
    m_scale = 0.0
    steps = []
    for istep in range(h_steps.shape[0]):
        hv = h_steps[istep]
        source = ng.LinearForm(fes)
        source += ng.InnerProduct(ng.CoefficientFunction(tuple(hv)), vf) * ng.dx
        source.Assemble()
        rhs_src = np.asarray(source.vec.FV().NumPy(), dtype=float).copy()
        d_prev = None
        for it in range(nl_maxit):
            m_new = solve_W0(rhs_src + PT @ s_el.T.ravel(), x0=m)
            d_now = float(np.linalg.norm(m_new - m))
            m = m_new
            m_scale = max(m_scale, float(np.linalg.norm(m_new)))
            rel = d_now / (m_scale + 1e-30)
            M_el = (P @ m).reshape(3, n_el).T / vol_el[:, None]
            if eps_loop > 0.0:
                M_rms = float(np.sqrt(np.mean(np.sum(M_el ** 2, axis=1))))
                M_fed = M_el + eps_loop * (M_rms / L_rms) * L_el
            else:
                M_fed = M_el
            B0 = MU0 * (M_fed + hv[None, :]) if B_cache is None else B_cache
            B_cache, H_el = _solve_pointwise_B(material, states, M_fed, B0)
            s_el = nu0 * M_fed - H_el
            nit = it + 1
            if it > 0 and rel < nl_tol:
                q = (d_now / d_prev) if (d_prev is not None and d_prev > 0.0) else 0.0
                if q < 1.0 and rel * max(1.0, q / (1.0 - q)) < nl_tol:
                    break
            d_prev = d_now
        else:
            return dict(status="no_converge", fail_step=istep, steps=steps)

        b_peak = float(np.max(np.linalg.norm(B_cache, axis=1)))
        if b_peak > b_max * (1.0 + 1e-6):
            return dict(status="runaway_out_of_range", fail_step=istep,
                        b_peak=b_peak, b_max=b_max, steps=steps)
        states = material.commit(B_cache, states)
        steps.append(dict(
            H_avg_z=float((w_el * H_el[:, 2]).sum()),
            B_avg_z=float((w_el * B_cache[:, 2]).sum()),
            b_peak=b_peak, iters=int(nit)))
    return dict(status="ok", steps=steps)


def test_loop_pollution_corrupts_binput_hysteresis():
    K, eta, tables = _load_play()
    h_steps = _drive()
    results = dict(
        description="Controlled loop injection into the B-input constitutive law "
                    "(collocation failure mode reproduced on the loop-free HDiv solver)",
        timestamp=datetime.now().isoformat(), hostname=platform.node(),
        nx=NX, H0_A_per_m=H0, seed=SEED, eps_cases=[])

    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=NX, ny=NX, nz=NX,
                                    mapping=lambda x, y, z: (x - .5, y - .5, z - .5))

        ref = vim.SolveHysteresis(mesh, h_steps,
                                  material=PlayHysteresisMaterial(K, eta, tables))
        Hz = np.array([s["H_avg"][2] for s in ref["steps"]])
        Bz = np.array([s["B_avg"][2] for s in ref["steps"]])
        area_ref = _loop_area(Hz[N_RAMP:], Bz[N_RAMP:])
        bpk_ref = float(max(np.max(np.linalg.norm(np.asarray(s["B"]), axis=1))
                            for s in ref["steps"]))
        it_ref = sum(int(s["iters"]) for s in ref["steps"])
        results.update(production_area=area_ref, production_b_peak=bpk_ref,
                       production_iters_sum=it_ref)

        by_eps = {}
        for eps in EPS_LIST:
            r = _run_harness(mesh, h_steps, PlayHysteresisMaterial(K, eta, tables), eps)
            case = dict(eps=eps, status=r["status"])
            if r["status"] == "ok":
                Hz = np.array([s["H_avg_z"] for s in r["steps"]])
                Bz = np.array([s["B_avg_z"] for s in r["steps"]])
                case.update(area=_loop_area(Hz[N_RAMP:], Bz[N_RAMP:]),
                            b_peak=max(s["b_peak"] for s in r["steps"]),
                            iters_sum=sum(s["iters"] for s in r["steps"]),
                            iters_max=max(s["iters"] for s in r["steps"]))
            else:
                case.update(fail_step=r.get("fail_step"), b_peak=r.get("b_peak"))
            by_eps[eps] = case
            results["eps_cases"].append(case)

    # (1) harness fidelity: the eps=0 leg IS the production iteration.  Two
    # INDEPENDENTLY converged nl_tol=1e-3 solves scatter by up to ~1e-3 in a
    # thread-schedule-dependent way (the per-region flake lesson: never gate
    # inside the outer-tolerance band), so the fidelity tolerance sits above
    # that band and far below the smallest injection effect (eps=0.05 shifts
    # the area by 1.2e-2 and peak |B| by +19%).
    c0 = by_eps[0.0]
    assert c0["status"] == "ok"
    assert abs(c0["area"] - area_ref) / area_ref < 2e-3, (
        "harness drifted from vim.SolveHysteresis (eps=0 area %.6f vs %.6f)"
        % (c0["area"], area_ref))
    assert abs(c0["b_peak"] - bpk_ref) / bpk_ref < 2e-3

    # (2) monotone harm while still converging: spurious saturation + iteration growth
    c5, c10 = by_eps[0.05], by_eps[0.1]
    assert c5["status"] == "ok" and c10["status"] == "ok"
    assert c5["b_peak"] > c0["b_peak"] * 1.08, (
        "5%% loop fraction should inflate peak |B| by >8%% (got %.3f vs %.3f)"
        % (c5["b_peak"], c0["b_peak"]))
    assert c10["b_peak"] > c0["b_peak"] * 1.25
    assert c5["iters_sum"] > 1.4 * c0["iters_sum"], (
        "5%% loop fraction should inflate polarization iterations by >1.4x "
        "(got %d vs %d)" % (c5["iters_sum"], c0["iters_sum"]))
    assert c10["iters_sum"] > 1.4 * c0["iters_sum"]

    # (3) breakdown at collocation-realistic loop fractions
    assert by_eps[0.2]["status"] != "ok", (
        "20%% loop fraction unexpectedly converged in range: %s" % by_eps[0.2])
    assert by_eps[0.3]["status"] != "ok"

    (HERE / "loop_pollution_binput.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-s"]))
