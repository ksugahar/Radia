"""Closed-PEC frequency sweep and Foster fit against analytical TE modes.

Historical comparison with the Kameari iteration was retired.
The independent frequency sweep, Foster fit and analytical reference remain.
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
    print("=== Frequency sweep and Foster fit on closed PEC cuboid ===")
    print(f"  Mesh: h_cond={H_COND*1000} mm, order={ORDER}\n")
    print(f"Reference: tau_TE_z(1,1,0) analytical = {tau_TE_110_us:.3f} us\n")

    geo = build_geo()
    print("Generating mesh...")
    mesh = Mesh(geo.GenerateMesh(maxh=H_COND))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}\n")

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

    print("\n*** CONCLUSION ***")
    if 4 in foster_fits:
        f_match = abs(foster_fits[4]["taus_us"][0] - tau_TE_110_us)/tau_TE_110_us < 0.05
        if f_match:
            print("(A) NG freq sweep + Foster fit RECOVERS the analytical answer (<5% off).")
            print("    -> Forward direct method works on closed PEC.")
        else:
            print("(A) NG freq sweep + Foster fit DOES NOT match analytical "
                  f"({foster_fits[4]['taus_us'][0]:.2f} vs {tau_TE_110_us:.2f}).")

    out = {
        "method": "Frequency sweep and Foster fit on closed-PEC cuboid",
        "h_cond_mm": H_COND*1000, "order": ORDER, "ne": mesh.ne,
        "tau_TE_110_analytical_us": tau_TE_110_us,
        "freq_sweep_results": fs_results,
        "foster_fits": foster_fits,
    }
    out_path = Path(__file__).parent / "cuboid_521_kameari_vs_freqsweep_PEC.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
