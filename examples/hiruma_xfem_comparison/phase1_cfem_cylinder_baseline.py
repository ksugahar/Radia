"""Phase 1 — CFEM baselines for the Hiruma 2023 cylinder benchmark.

Reference:  Hiruma 2023, IEEE TMag 59(5), DOI 10.1109/TMAG.2023.3246629,
            §III.A "Skin Effect in Cylindrical Conductor", Figs. 2-4.

Geometry / material (Hiruma's §III.A):
    2D disk Omega, radius r = 0.01 m
    sigma = 1e6 S/m, mu_r = 1
    Uniform current density J0 = sigma (so E_z = 1 V/m on the boundary)
    Governing: -div(nu*grad A_z) + j*omega*sigma*A_z = J0
    BC: A_z = 0 on outer circle Gamma_D
    Sweep: r/delta in [0.1, 15] (i.e. omega = 2*delta^-2/(sigma*mu))

Goal of this phase:
    Reproduce Hiruma Fig. 4 (Joule loss density and magnetic energy
    density vs r/delta) using NGSolve CFEM on TWO meshes:
        - coarse (target ~37 DOF: matches Hiruma's 37 in Fig. 2(b))
        - fine reference (target ~9489 DOF: matches Hiruma)
    Validate against the analytical Bessel solution for the disk.

Bessel reference (closed form for cylinder, J0 = sigma, A_z = 0 at r=R):
    A_z(r) = (1/(j*omega)) * (1 - I_0(gamma*r) / I_0(gamma*R))
    where gamma = sqrt(j*omega*sigma*mu_0) = (1+j)/delta
    (Derived: pde becomes (1/r)d/dr(r dA/dr) - gamma^2 A = -j*omega*mu*sigma * (sigma * 1)/(-j*omega) = ...
     simpler: H_phi from -dA_z/dr, J_z = -j*omega*sigma*A_z + J0,
              integrate J_z^2/sigma over disk for Joule loss.)

Output:
    results_phase1.npz : r_over_delta, P_coarse, P_fine, P_bessel,
                         E_coarse, E_fine, E_bessel
    plus log file with iteration counts and DOF stats
"""

import math
import numpy as np
import json
from pathlib import Path

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction, Preconditioner,
    grad, dx, ds, CoefficientFunction, Integrate, x, y, sqrt,
    Conj,
)
from netgen.geom2d import SplineGeometry
from scipy.special import iv  # modified Bessel I_n (cp932-safe)


# ---------- physics constants ----------
R     = 0.01          # m, cylinder radius
SIGMA = 1.0e6         # S/m
MU0   = 4.0e-7 * math.pi
NU    = 1.0 / MU0     # magnetic reluctivity
J0    = SIGMA         # uniform external current density (A/m^2)


# =====================================================================
# Bessel analytical reference (no FE)
# =====================================================================

def bessel_reference(r_over_delta):
    """Closed-form Joule loss and magnetic energy in the disk.

    The PDE is:
        -nu * Laplacian(A) + j*omega*sigma*A = J0
    Rewriting with k^2 = -j*omega*sigma*mu_0 = -gamma^2 (gamma = (1+j)/delta):
        Laplacian(A) - gamma^2 A = -mu_0 * J0
        => A_p = mu_0 * J0 / gamma^2 = J0 / (j*omega*sigma)    (particular)
        => A_h = c * I_0(gamma * r)                            (homogeneous)
    BC A(R) = 0 forces c = -A_p / I_0(gamma R).

    Joule loss density (per unit length in z):
        P = (1/2) int_disk sigma |E_z|^2 dA, E_z = -j*omega*A_z (no V drop, 2D)
        but Hiruma uses J0 fixed so we just integrate sigma*|j*omega*A - J0/sigma|^2 / sigma.
    For our purposes report only the trend; Phase 2-3 compare absolute scales.

    Returns: dict with P_density_avg [W/m^3], E_density_avg [J/m^3].
    """
    delta = R / r_over_delta
    omega = 2.0 / (delta * delta * SIGMA * MU0)
    gamma = (1 + 1j) / delta

    A_p = J0 / (1j * omega * SIGMA)          # particular
    c   = -A_p / iv(0, gamma * R)             # so that A(R) = 0

    # Evaluate radial profile on a Gauss-Legendre quadrature
    n_quad = 200
    xi, wi = np.polynomial.legendre.leggauss(n_quad)
    r_pts  = R * 0.5 * (xi + 1)
    r_wts  = R * 0.5 * wi

    A      = A_p + c * iv(0, gamma * r_pts)   # complex A_z(r)
    E      = -1j * omega * A                   # E_z = -j w A (Coulomb gauge, V=0)
    J_tot  = SIGMA * E                          # J_z = sigma E
    # Joule loss density (per unit volume) at this r: |J|^2 / (2 sigma)
    P_r    = np.abs(J_tot) ** 2 / (2 * SIGMA)
    # Magnetic field H_phi = -(1/mu_0) dA_z/dr; B_phi = -dA_z/dr ; energy density |B|^2/(2 mu_0)
    # dA/dr = c * gamma * I_1(gamma r)
    dA_dr  = c * gamma * iv(1, gamma * r_pts)
    B      = -dA_dr
    W_r    = np.abs(B) ** 2 / (2 * MU0)

    # Disk averages: integral 2*pi*r * f(r) dr / (pi R^2)
    P_avg  = np.sum(2 * math.pi * r_pts * P_r * r_wts) / (math.pi * R * R)
    W_avg  = np.sum(2 * math.pi * r_pts * W_r * r_wts) / (math.pi * R * R)
    return {"P_density_avg": float(P_avg),
            "E_density_avg": float(W_avg),
            "omega": float(omega),
            "delta": float(delta)}


# =====================================================================
# NGSolve CFEM disk problem
# =====================================================================

def make_disk_mesh(maxh, name="disk"):
    """Generate a unit-disk mesh of radius R with prescribed maxh."""
    geo = SplineGeometry()
    geo.AddCircle((0, 0), R, leftdomain=1, rightdomain=0, bc="outer")
    return Mesh(geo.GenerateMesh(maxh=maxh))


def cfem_solve(mesh, omega):
    """Standard CFEM solve of  -nu*Laplace A + j*omega*sigma A = J0  with A|outer = 0.

    Returns (P_density_avg, E_density_avg, ndof, dict_diag).
    """
    fes = H1(mesh, order=1, complex=True, dirichlet="outer")
    u, v = fes.TnT()

    a = BilinearForm(fes, symmetric=True)
    a += NU * grad(u) * grad(v) * dx
    a += 1j * omega * SIGMA * u * v * dx
    a.Assemble()

    f = LinearForm(fes)
    f += J0 * v * dx
    f.Assemble()

    A = GridFunction(fes)
    A.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # Post-process: Joule loss + magnetic energy as area-averages
    E_z = -1j * omega * A
    J_z = SIGMA * E_z              # the eddy current (J0 is the source, but the actual J = -j w sigma A + 0 if E_z accounts for it; double-check)
    # Wait: -nu*Laplace(A) + j*omega*sigma*A = J0
    # => the field equation is sigma*E_z = J_z, where E_z = -j*omega*A (since V=0 in 2D).
    # However J0 is the SOURCE, NOT the eddy current.  The total current is
    #   J_total = J0 + sigma * E_z = J0 - j w sigma A.
    # Joule loss per volume = (1/2) sigma |E_z|^2 = (1/2) |J_total - J0|^2 / sigma
    # In Hiruma's setup with E_z = 1 prescribed at boundary, J0 = sigma * 1 = sigma,
    # and we want sigma |E_z|^2 / 2.
    P_cf = 0.5 * SIGMA * (E_z * Conj(E_z))   # CF for Joule loss density [W/m^3]
    B_x  = grad(A)[1]                         # B = curl(A_z e_z) = (dA/dy, -dA/dx)
    B_y  = -grad(A)[0]
    W_cf = (B_x * Conj(B_x) + B_y * Conj(B_y)) / (2 * MU0)

    # Disk area = pi*R^2; .real because the CF includes |.|^2 already (complex but real-valued)
    area = math.pi * R * R
    P_avg = Integrate(P_cf, mesh).real / area
    W_avg = Integrate(W_cf, mesh).real / area
    return P_avg, W_avg, fes.ndof, {}


# =====================================================================
# Main: sweep r/delta and tabulate
# =====================================================================

def main():
    out_dir = Path(__file__).parent
    print(f"output dir: {out_dir}", flush=True)
    print(f"R = {R*1000:.2f} mm, sigma = {SIGMA:.2e} S/m, mu_r = 1", flush=True)
    print()

    # Build meshes once, reuse for all frequencies
    print("=== mesh generation ===", flush=True)
    mesh_coarse = make_disk_mesh(maxh=R / 3.5,  name="coarse")  # target ~37 DOF (Hiruma 2023 §III.A)
    mesh_fine   = make_disk_mesh(maxh=R / 30.0, name="fine")    # target ~10k DOF
    print(f"  coarse: ne={mesh_coarse.ne}, nv={mesh_coarse.nv}", flush=True)
    print(f"  fine:   ne={mesh_fine.ne},   nv={mesh_fine.nv}",   flush=True)

    # Frequency sweep (15 points to start, log-spaced through skin transition)
    r_over_delta_list = np.array([0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0,
                                   7.0, 10.0, 12.0, 15.0])
    results = []
    print()
    print("=== frequency sweep (CFEM coarse + fine + Bessel) ===", flush=True)
    print(f"{'r/d':>6} {'omega':>12} {'P_bessel':>12} {'P_coarse':>12} {'P_fine':>12} "
          f"{'err_coarse':>10} {'err_fine':>10}", flush=True)
    for rd in r_over_delta_list:
        # Bessel reference
        ref = bessel_reference(rd)
        omega = ref["omega"]
        # CFEM solves
        P_c, E_c, n_c, _ = cfem_solve(mesh_coarse, omega)
        P_f, E_f, n_f, _ = cfem_solve(mesh_fine,   omega)
        # Joule loss errors vs Bessel reference
        err_c = (P_c - ref["P_density_avg"]) / ref["P_density_avg"] * 100
        err_f = (P_f - ref["P_density_avg"]) / ref["P_density_avg"] * 100
        results.append({
            "r_over_delta": float(rd),
            "omega":       float(omega),
            "P_bessel":    ref["P_density_avg"],
            "E_bessel":    ref["E_density_avg"],
            "P_coarse":    P_c, "E_coarse": E_c, "ndof_coarse": n_c,
            "P_fine":      P_f, "E_fine":   E_f, "ndof_fine":   n_f,
            "err_coarse_pct": float(err_c),
            "err_fine_pct":   float(err_f),
        })
        print(f"{rd:6.2f} {omega:12.3e} {ref['P_density_avg']:12.4e} "
              f"{P_c:12.4e} {P_f:12.4e} {err_c:+9.2f}% {err_f:+9.2f}%",
              flush=True)

    # Save
    out_npz = out_dir / "results_phase1.npz"
    np.savez(out_npz,
             r_over_delta = np.array([r["r_over_delta"] for r in results]),
             omega        = np.array([r["omega"]        for r in results]),
             P_bessel     = np.array([r["P_bessel"]     for r in results]),
             E_bessel     = np.array([r["E_bessel"]     for r in results]),
             P_coarse     = np.array([r["P_coarse"]     for r in results]),
             P_fine       = np.array([r["P_fine"]       for r in results]),
             E_coarse     = np.array([r["E_coarse"]     for r in results]),
             E_fine       = np.array([r["E_fine"]       for r in results]),
             ndof_coarse  = np.array([r["ndof_coarse"]  for r in results]),
             ndof_fine    = np.array([r["ndof_fine"]    for r in results]))
    print()
    print(f"saved: {out_npz}", flush=True)

    out_json = out_dir / "results_phase1.json"
    with open(out_json, "w") as f:
        json.dump({"meta": {"R_m": R, "sigma": SIGMA, "mu_r": 1,
                            "ndof_coarse": int(results[0]["ndof_coarse"]),
                            "ndof_fine":   int(results[0]["ndof_fine"])},
                   "results": results}, f, indent=2)
    print(f"saved: {out_json}", flush=True)

    print()
    print("=== Phase 1 summary ===", flush=True)
    print(f"  CFEM coarse DOF: {results[0]['ndof_coarse']} "
          f"(Hiruma 37)", flush=True)
    print(f"  CFEM fine   DOF: {results[0]['ndof_fine']} "
          f"(Hiruma 9489)", flush=True)
    err_fine_at_15 = next(r['err_fine_pct'] for r in results
                          if r['r_over_delta'] == 15.0)
    err_coarse_at_15 = next(r['err_coarse_pct'] for r in results
                            if r['r_over_delta'] == 15.0)
    print(f"  at r/d=15: CFEM-coarse err = {err_coarse_at_15:+.2f}%, "
          f"CFEM-fine err = {err_fine_at_15:+.2f}% (vs Bessel)", flush=True)


if __name__ == "__main__":
    main()
