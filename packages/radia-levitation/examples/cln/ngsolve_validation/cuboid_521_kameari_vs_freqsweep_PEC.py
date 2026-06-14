"""ISOLATION TEST: Is Kameari iteration itself healthy in 3D?

Closed-PEC cuboid 5x2x1 has analytical reference:
  tau_TE_z(1,1,0) = mu_0 sigma / (pi^2 (1/a^2 + 1/b^2)) = 25.46 us

Two independent approaches on the SAME closed-PEC mesh:
  (A) NGSolve complex freq sweep at 30 freqs -> Foster fit -> tau_lead^Foster
  (B) Kameari iteration -> R_n, L_n -> tau_n^Cauer = L_n/R_n
  (C) Analytical TE modes -> tau_lead^analytical = 25.46 us

Comparison:
  - (A) vs (C): Does NG freq sweep produce the analytical answer? (sanity)
  - (B) vs (A) and (C): Does Kameari produce consistent results?

If (A)~(C) and Kameari Cauer ladder reconstructs Y(iw) matching (A),
then KAMEARI ITERATION IS HEALTHY for closed-PEC 3D.
If Kameari ladder doesn't reconstruct correct Y(iw),
then KAMEARI ITERATION IS BROKEN even for closed-PEC 3D.

This is the proper test, isolating Kameari from open-domain inner solver.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import numpy as np
import json, time
from pathlib import Path
from scipy.optimize import least_squares

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

# Analytical TE_z(1,1,0)
tau_TE_110_us = mu0 * sigma_Cu / (pi**2 * (1/ax**2 + 1/ay**2)) * 1e6
print(f"\n*** Analytical reference: tau_TE_z(1,1,0) = {tau_TE_110_us:.3f} us ***\n")

# Mesh settings
H_COND = 0.30e-3
ORDER = 2
N_STAGES = 8


def build_spanning_tree(mesh):
    nv = mesh.nv
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges


def build_geo():
    """Just the cuboid (closed PEC: A x n = 0 on conductor surface)."""
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND
    return OCCGeometry(cuboid)


def kameari_closed_PEC(mesh, n_stages):
    """Kameari iteration on closed-PEC cuboid."""
    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")
    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  Tree-cotree masked {masked} edges, active = {sum(fd)}")

    sigma_inv_cf = 1.0 / sigma_Cu
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx
    print("  Assembling closed-PEC operator...")
    with TaskManager():
        a.Assemble()
        inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

    diag = []
    J_cf = sigma_Cu * A_ext
    Apot = None

    for n in range(n_stages):
        f = LinearForm(fes)
        f += J_cf * v * dx
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes)
        with TaskManager():
            gfA.vec.data = inv * f.vec

        R_inv = float(Integrate(J_cf * J_cf * sigma_inv_cf * dx, mesh))
        if abs(R_inv) < 1e-30:
            break
        Rn = 1.0 / R_inv
        if Apot is None:
            Apot = Rn * gfA
        else:
            Apot = Apot + Rn * gfA
        Ln_int = float(Integrate(J_cf * Apot * dx, mesh))
        Ln = Rn * Ln_int
        tau_us = Ln/Rn * 1e6 if Rn > 0 and Ln > 0 else float('nan')
        signL = "+" if Ln > 0 else "-"
        diag.append({"n": n, "R_n": Rn, "L_n": Ln,
                     "tau_us": tau_us if Ln > 0 else None,
                     "sign_L": signL})
        print(f"  [{n}] R={Rn:.4e}, L={Ln:.4e} H, tau={tau_us:.3f} us")

        if Ln <= 0:
            break
        J_cf = J_cf - sigma_Cu * Apot / Ln

    return diag


def freq_sweep_closed_PEC(mesh, freqs):
    """Complex freq sweep on closed-PEC cuboid: solve at each i*omega.
    Get m(iw) and convert to Y_zz."""
    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface",
                complex=True, nograds=True)
    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False

    A_ext = CoefficientFunction((-y, x, 0)) * 0.5
    u, v = fes.TnT()
    results = []

    for k, f_hz in enumerate(freqs):
        omega = 2 * pi * f_hz
        s = 1j * omega
        a = BilinearForm(fes)
        a += (1.0/mu0) * curl(u) * curl(v) * dx
        a += s * sigma_Cu * u * v * dx
        rhs = LinearForm(fes)
        rhs += -s * sigma_Cu * A_ext * v * dx
        with TaskManager():
            a.Assemble(); rhs.Assemble()
            inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        gf = GridFunction(fes)
        with TaskManager():
            gf.vec.data = inv * rhs.vec
        # m_z from total J
        J_total = -s * sigma_Cu * (A_ext + gf)
        m_int = Integrate((x * J_total[1] - y * J_total[0]) * 0.5 * dx, mesh)
        m_re, m_im = float(m_int.real), float(m_int.imag)
        results.append({"freq": f_hz, "m_re": m_re, "m_im": m_im})
        if k % 5 == 0:
            print(f"  [{k+1}/{len(freqs)}] f={f_hz:.2e}: m={m_re:.3e}+{m_im:.3e}j")

    return results


def foster_fit(s_arr, M_arr, N):
    init = np.logspace(np.log10(0.5e-6), np.log10(50e-6), N)
    def residual(x):
        taus = np.exp(x[:N])
        a = np.exp(x[N:2*N]) * np.sign(M_arr.real[0])
        fit = np.zeros_like(s_arr, dtype=complex)
        for tk, ak in zip(taus, a):
            fit += ak / (1 + s_arr*tk)
        return np.concatenate([(fit-M_arr).real/np.abs(M_arr),
                               (fit-M_arr).imag/np.abs(M_arr)])
    M0 = abs(M_arr[0])
    x0 = np.concatenate([np.log(init), np.log(np.ones(N)*M0/N)])
    res = least_squares(residual, x0, method='trf', max_nfev=5000)
    taus = np.exp(res.x[:N])
    a = np.exp(res.x[N:2*N]) * np.sign(M_arr.real[0])
    order = np.argsort(-taus)
    M_fit = np.zeros_like(s_arr, dtype=complex)
    for tk, ak in zip(taus[order], a[order]):
        M_fit += ak / (1 + s_arr*tk)
    rms = np.sqrt(np.mean(np.abs(M_fit-M_arr)**2 / np.abs(M_arr)**2)) * 100
    return taus[order], a[order], rms


def main():
    ngsglobals.msg_level = 0
    print("=== ISOLATION TEST: Kameari vs FreqSweep on closed PEC cuboid ===")
    print(f"  Mesh: h_cond={H_COND*1000} mm, order={ORDER}\n")
    print(f"Reference: tau_TE_z(1,1,0) analytical = {tau_TE_110_us:.3f} us\n")

    geo = build_geo()
    print("Generating mesh...")
    mesh = Mesh(geo.GenerateMesh(maxh=H_COND))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}\n")

    # === Run Kameari ===
    print("=" * 60)
    print("KAMEARI iteration on closed-PEC cuboid")
    print("=" * 60)
    t0 = time.time()
    kam_diag = kameari_closed_PEC(mesh, N_STAGES)
    print(f"  Kameari done ({time.time()-t0:.1f}s)\n")

    # === Run freq sweep ===
    print("=" * 60)
    print("NGSolve complex frequency sweep on closed-PEC cuboid")
    print("=" * 60)
    # Same 30 freqs as ELF test
    freqs = np.logspace(2, 8, 30).tolist()
    t0 = time.time()
    fs_results = freq_sweep_closed_PEC(mesh, freqs)
    print(f"  Freq sweep done ({time.time()-t0:.1f}s)\n")

    # === Foster fit ===
    s_arr = np.array([1j * 2*pi*r["freq"] for r in fs_results])
    m_arr = np.array([r["m_re"] + 1j*r["m_im"] for r in fs_results])
    M_arr = -m_arr / (mu0 * s_arr)

    print("=" * 60)
    print("FOSTER FIT (closed PEC, NG freq sweep)")
    print("=" * 60)
    foster_fits = {}
    for N in [4, 6, 8]:
        taus, a, rms = foster_fit(s_arr, M_arr, N)
        foster_fits[N] = {"taus_us": (taus*1e6).tolist(),
                          "a": a.tolist(), "rms": float(rms)}
        print(f"  N={N}: tau = {[f'{t*1e6:.3f}' for t in taus]} us, rms={rms:.3f}%")

    # === Critical comparisons ===
    print("\n" + "=" * 60)
    print("CRITICAL COMPARISONS (closed PEC cuboid)")
    print("=" * 60)
    print(f"\n  Analytical:        tau_TE(1,1,0)  = {tau_TE_110_us:.3f} us")
    if 4 in foster_fits:
        tau_F = foster_fits[4]["taus_us"][0]
        print(f"  NG Foster fit:     tau_lead^F     = {tau_F:.3f} us "
              f"(N=4, rms={foster_fits[4]['rms']:.3f}%)")
        f_vs_anal = abs(tau_F - tau_TE_110_us) / tau_TE_110_us * 100
        print(f"    -> {f_vs_anal:.2f}% deviation from analytical")
    if kam_diag and kam_diag[0]['tau_us']:
        tau_K = kam_diag[0]['tau_us']
        max_K = max(d['tau_us'] for d in kam_diag if d['tau_us'])
        print(f"  Kameari Cauer:     tau_0^C        = {tau_K:.3f} us")
        print(f"  Kameari Cauer:     max(tau_n^C)   = {max_K:.3f} us")
        print(f"    -> max(tau^C) / tau_anal = {max_K/tau_TE_110_us:.4f}")
        print(f"    -> tau_0^C / tau_anal    = {tau_K/tau_TE_110_us:.4f}")

    print("\n*** CONCLUSION ***")
    if 4 in foster_fits:
        f_match = abs(foster_fits[4]["taus_us"][0] - tau_TE_110_us)/tau_TE_110_us < 0.05
        if f_match:
            print("(A) NG freq sweep + Foster fit RECOVERS the analytical answer (<5% off).")
            print("    -> Forward direct method works on closed PEC.")
        else:
            print("(A) NG freq sweep + Foster fit DOES NOT match analytical "
                  f"({foster_fits[4]['taus_us'][0]:.2f} vs {tau_TE_110_us:.2f}).")
    if kam_diag and kam_diag[0]['tau_us']:
        tau_K = kam_diag[0]['tau_us']
        if abs(tau_K - tau_TE_110_us) / tau_TE_110_us < 0.05:
            print(f"(B) Kameari Cauer tau_0 ALSO matches analytical Foster lead.")
            print("    -> Kameari iteration HEALTHY on 3D closed PEC.")
        else:
            print(f"(B) Kameari Cauer tau_0={tau_K:.2f} differs from Foster lead {tau_TE_110_us:.2f}.")
            print("    -> Cauer != Foster (different decompositions); needs Y(s) reconstruction to verify.")

    out = {
        "method": "Isolation test: Kameari vs FreqSweep on closed-PEC cuboid",
        "h_cond_mm": H_COND*1000, "order": ORDER, "ne": mesh.ne,
        "tau_TE_110_analytical_us": tau_TE_110_us,
        "kameari_stages": kam_diag,
        "freq_sweep_results": fs_results,
        "foster_fits": foster_fits,
    }
    out_path = Path(__file__).parent / "cuboid_521_kameari_vs_freqsweep_PEC.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
