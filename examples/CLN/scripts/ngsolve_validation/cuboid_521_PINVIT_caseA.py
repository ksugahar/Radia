"""5x2x1 mm Cu cuboid CLN — PINVIT eigenvalue extraction (Case A).

Avoid Kameari recursion's catastrophic cancellation by computing eigenmodes
of curl-curl directly:
  curl-curl phi_n = lambda_n^2 mu sigma phi_n,  on HCurl, with full Dirichlet.

Foster expansion (Case A, mz=0 only ansatz holds since A=(0,0,Az(x,y))):
  Y(s) = sigma * sum |beta_n|^2 / (1 + s tau_n)
  beta_n = <phi_n, J/sigma> = <phi_n, e_z>
  tau_n = mu sigma / lambda_n^2
  R_n^Foster = 1/(sigma |beta_n|^2)
  L_n^Foster = mu / (lambda_n^2 |beta_n|^2)

Then convert to Cauer-I via moments (same as Mathematica) at numerical
precision with N eigenmodes. Each new eigenmode adds one resolvable Cauer
stage (rule of thumb: need ~2N+ Cauer stages from N+ Foster modes).

Strategy:
  - Use NGSolve PINVIT to extract first N_EIGEN smallest eigenpairs
  - Compute beta_n by integration against J source
  - Build Y(s) Taylor series, extract Cauer-I via numpy polynomial division
  - Compare to Kameari v3 + Mathematica Nagamine reference
"""
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction, Preconditioner,
    CoefficientFunction, curl, dx, x, y, z, BitArray, BND,
    Integrate, TaskManager, ngsglobals,
)
from ngsolve.solvers import PINVIT
from collections import deque
from math import pi
import json
import time
from pathlib import Path

import numpy as np

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_box = ax * ay * az
N_EIGEN = 30          # first 30 eigenmodes
H_COND = 0.20e-3
ORDER = 2             # PINVIT memory cost lower at order 2
PINVIT_MAXIT = 600


def build_cuboid():
    box = Box(Pnt(0, 0, 0), Pnt(ax, ay, az))
    box.mat("conductor").bc("conductor_surface")
    box.maxh = H_COND
    return OCCGeometry(box)


def classify_vertices(mesh):
    boundary_v = set()
    for el in mesh.Elements(BND):
        for v in el.vertices:
            boundary_v.add(v.nr)
    return set(range(mesh.nv)) - boundary_v, boundary_v


def build_spanning_tree_interior(mesh):
    interior_v, _ = classify_vertices(mesh)
    nv = mesh.nv
    visited = [False] * nv
    for v in range(nv):
        if v not in interior_v:
            visited[v] = True
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        if v0 in interior_v and v1 in interior_v:
            adj[v0].append((v1, ed.nr))
            adj[v1].append((v0, ed.nr))
    for start in interior_v:
        if visited[start]:
            continue
        visited[start] = True
        queue = deque([start])
        while queue:
            v = queue.popleft()
            for vn, edn in adj[v]:
                if not visited[vn]:
                    visited[vn] = True
                    tree_edges.append(edn)
                    queue.append(vn)
    return tree_edges


def main():
    ngsglobals.msg_level = 0
    print(f"5x2x1 mm Cu cuboid - PINVIT eigenvalue extraction (Case A)")
    print(f"  h = {H_COND*1000} mm, order = {ORDER}, N_eigen = {N_EIGEN}")
    geo = build_cuboid()
    mesh = Mesh(geo.GenerateMesh(maxh=H_COND))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}")

    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")

    # Tree-cotree mask in-place into fes.FreeDofs()
    # so PINVIT (which uses fes.FreeDofs internally) sees the masked DoFs.
    tree_edges = build_spanning_tree_interior(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  Tree-cotree masked {masked} DoFs (in-place), active = {sum(fd)}")

    u, v = fes.TnT()

    # Stiffness: (1/mu) curl-curl
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx
    # Mass: sigma * (.,.) (for eigenvalue lambda^2 of -nabla^2 / mu = lambda^2 sigma A)
    # Equivalent: lambda^2 = mu sigma * generalized_eig_value
    m = BilinearForm(fes)
    m += sigma_Cu * u * v * dx

    print("  Assembling...")
    pre = Preconditioner(a, type="local")
    with TaskManager():
        a.Assemble()
        m.Assemble()

    # Initial multivector for PINVIT
    print(f"  Running PINVIT for {N_EIGEN} smallest eigenpairs...")
    t0 = time.time()
    gf_modes = GridFunction(fes, multidim=N_EIGEN, name="modes")

    with TaskManager():
        evals = PINVIT(a.mat, m.mat, pre=pre.mat, num=N_EIGEN,
                       maxit=PINVIT_MAXIT, printrates=False)

    print(f"  PINVIT done in {time.time()-t0:.1f}s")

    # Parse eigenvalues
    if isinstance(evals, tuple):
        evals_vec = evals[0]
        evec_multi = evals[1]
    else:
        evals_vec = evals
        evec_multi = gf_modes.vecs

    if hasattr(evals_vec, "NumPy"):
        evals_arr = np.array(evals_vec.NumPy())
    elif hasattr(evals_vec, "FV"):
        evals_arr = np.array(list(evals_vec.FV()))
    else:
        s = str(evals_vec).strip()
        evals_arr = np.array([float(x.strip()) for x in s.split("\n") if x.strip()])
    print(f"  evals[0..5] = {evals_arr[:6]}")

    # eigenvalue here is generalized: a v = ev * m v -> v^T a v / v^T m v = ev
    # a = (1/mu) curl-curl, m = sigma * Identity
    # So ev = (1/mu) lambda_curl_curl^2 / sigma -> lambda_curl_curl^2 = ev * mu * sigma
    # But for Foster: tau_n = mu sigma / lambda^2. So tau_n = mu sigma / (ev mu sigma) = 1/ev.
    # And lambda^2 = ev mu sigma.

    print(f"\n=== Foster coefficients (Case A: J/sigma = e_z) ===")
    J_unit = CoefficientFunction((0, 0, 1))   # J/sigma direction

    foster_R = []
    foster_L = []
    foster_tau = []
    EV_MIN = 1e-3   # filter spurious near-zero (kernel pollution) eigenvalues
    for k in range(N_EIGEN):
        ev = float(evals_arr[k])
        if ev < EV_MIN:
            print(f"  mode {k}: ev={ev:.2e} < {EV_MIN:.0e} (kernel/spurious), skip")
            continue
        gf_k = GridFunction(fes)
        gf_k.vec.data = evec_multi[k]
        # Normalize: phi^T M phi = 1 (M = sigma * I)
        nrm2 = float(Integrate(sigma_Cu * gf_k * gf_k * dx, mesh))
        if nrm2 < 1e-30:
            continue
        nrm = nrm2**0.5
        gf_k.vec.data = (1.0/nrm) * gf_k.vec   # syntax: scalar * BaseVector
        beta_k = float(Integrate(gf_k * J_unit * dx, mesh))
        beta_sq = beta_k**2
        if abs(beta_sq) < 1e-30:
            print(f"  mode {k}: ev={ev:.4e}, beta^2={beta_sq:.2e} (decoupled), skip")
            continue
        # Foster R, L
        # tau_n = 1/ev (from generalized eig: a v = ev * m v with M=sigma*I)
        # R_n^F = 1/(sigma * beta^2)
        # L_n^F = R_n^F * tau_n
        Rn = 1.0 / (sigma_Cu * beta_sq)
        tau_n = 1.0 / ev
        Ln = Rn * tau_n
        foster_R.append(Rn)
        foster_L.append(Ln)
        foster_tau.append(tau_n)
        print(f"  mode {k}: ev={ev:.4e}, beta^2={beta_sq:.4e}, "
              f"R_F={Rn:.4e}, L_F={Ln*1e9:.4e} nH, tau={tau_n*1e6:.4f} us")

    print(f"\nTotal {len(foster_R)} modes coupled to J_z")
    print(f"Y(0) = sum 1/R = {sum(1.0/r for r in foster_R):.4e}")
    print(f"Reference (analytic): sigma*V = {sigma_Cu * V_box:.4e}")

    # === Build Y(s) Taylor and extract Cauer-I via numpy polynomial division ===
    # M_k = sum_n |beta_n|^2 sigma * tau_n^k = sum (1/R_n) tau_n^k
    Kmax = 25
    inv_R = np.array([1.0/r for r in foster_R])
    tau_arr = np.array(foster_tau)
    M = np.zeros(Kmax + 1)
    for k in range(Kmax + 1):
        M[k] = np.sum(inv_R * tau_arr**k)
    print(f"\nMoments M_k:")
    for k in range(min(6, Kmax+1)):
        print(f"  M_{k} = {M[k]:.10e}")

    # Build Z(s) = 1/Y(s) where Y(s) = sum (-s)^k M_k (truncated)
    # Z(s) = R_0 + s L_1 / (1 + s ... ) ladder
    # Use Mathematica-equivalent algorithm: Cauer extraction from rational
    # via numpy polynomials (poly1d) — high precision via mpmath if needed
    try:
        import mpmath
        mpmath.mp.dps = 60
        # Convert M_k to mpmath
        Mmp = [mpmath.mpf(m) for m in M]
        # Y(s) polynomial coefficients (low-order first):
        # Y(s) = M_0 - s M_1 + s^2 M_2 - ...
        # For numpy.polynomial.polynomial, coefficients in increasing degree
        ycoef = [Mmp[k] * (-1)**k for k in range(Kmax+1)]

        # Z(s) = 1/Y(s) - work as rational num/den
        # We compute Z as polynomial-division Cauer-I extraction
        # using mpmath rational handling (or just polynomial of degree Kmax)
        # Cauer-I algorithm:
        #   Input: Z(s) = p(s)/q(s) (initialize Z = 1/Y so num=1, den=ycoef)
        #   Stage k: extract R_k = p(0)/q(0), L_k = ...
        #   ...
        from mpmath import mpf, mpc, mp
        def cauer_extract(num_coefs, den_coefs, n_stages):
            # num, den are lists of coefs, low-order first
            # Iterate
            stages = []
            num = list(num_coefs)
            den = list(den_coefs)
            for kk in range(n_stages):
                # Trim trailing zeros
                while len(num) > 1 and num[-1] == 0:
                    num.pop()
                while len(den) > 1 and den[-1] == 0:
                    den.pop()
                if not den or den[0] == 0:
                    print(f"  Stage {kk}: den(0)~0, stop")
                    break
                if not num:
                    print(f"  Stage {kk}: num=[], stop")
                    break
                R_k = num[0] / den[0]
                stages.append(("R", float(R_k)))
                # presid = num - R_k * den
                presid = [num[i] - R_k * den[i] for i in range(max(len(num), len(den)))]
                # pad to same length
                while len(presid) < max(len(num), len(den)):
                    presid.append(mpf(0))
                # pas = presid / s = shift down (drop first if zero)
                if abs(presid[0]) > 1e-50:
                    print(f"  Stage {kk}: presid(0) = {float(presid[0]):.2e}, not zero")
                pas = presid[1:]
                if not pas or pas[0] == 0:
                    print(f"  Stage {kk}: pas(0)~0, stop")
                    break
                L_k = pas[0] / den[0]
                stages.append(("L", float(L_k)))
                # nresid = den * L_k - pas
                nresid = [den[i]*L_k - (pas[i] if i < len(pas) else mpf(0))
                          for i in range(max(len(den), len(pas)))]
                while len(nresid) < max(len(den), len(pas)):
                    nresid.append(mpf(0))
                if abs(nresid[0]) > 1e-50:
                    print(f"  Stage {kk}: nresid(0) = {float(nresid[0]):.2e}, not zero")
                # nresid_shifted = nresid / s
                nresid_shifted = nresid[1:]
                # New: num = L_k * pas, den = nresid_shifted
                num = [L_k * pas[i] for i in range(len(pas))]
                den = nresid_shifted
            return stages

        # Initial: Z = 1/Y -> num = [1], den = ycoef
        stages = cauer_extract([mpf(1)], ycoef, 24)
        cauer_R = [s[1] for s in stages if s[0] == "R"]
        cauer_L = [s[1] for s in stages if s[0] == "L"]
    except ImportError:
        print("mpmath not available, skipping Cauer extraction")
        cauer_R, cauer_L = [], []

    print(f"\n=== Cauer-I from PINVIT Foster ({len(cauer_R)} stages) ===")
    for n in range(min(len(cauer_R), len(cauer_L))):
        print(f"  n={n}: R={cauer_R[n]:.6e},  L={cauer_L[n]*1e9:.6e} nH")

    out = {
        "method": "PINVIT eigenvalue extraction -> Foster -> Cauer-I via polynomial division",
        "geometry_mm": [ax*1000, ay*1000, az*1000],
        "h_mm": H_COND*1000, "order": ORDER, "ne": mesh.ne,
        "n_eigen_used": N_EIGEN,
        "n_modes_coupled": len(foster_R),
        "Y0_numerical": float(sum(inv_R)),
        "Y0_analytic": sigma_Cu * V_box,
        "foster_R": foster_R,
        "foster_L": foster_L,
        "foster_tau": foster_tau,
        "cauer_R": cauer_R,
        "cauer_L": cauer_L,
    }
    out_path = Path(__file__).parent / "cuboid_521_PINVIT_caseA_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
