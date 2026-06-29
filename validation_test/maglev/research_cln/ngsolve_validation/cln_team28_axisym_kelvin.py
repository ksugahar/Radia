"""cln_team28_axisym_kelvin.py — COMSOL TEAM 28 CLN with Kelvin two-half-disc
in axisymmetric NGSolve (standard A-formulation u = r·A_θ).

Geometry (Z-offset Kelvin two-half-disc, following the maintained
axisymmetric Kelvin convention in docs/kelvin/CONVENTION.md and the
validated material-factor API in radia.kelvin_source):
  Interior half-disc at (0, 0): radius R_K, contains conductor + air
  Exterior Kelvin half-disc at (0, z_offset): radius R_K, represents ∞
  z_offset = 5 × R_K (well-separated)
  Periodic BC identifies inner and outer kelvin arcs
  GND vertex at (0, z_offset, 0) → image of physical infinity

Conductor (cylinder R_DISK=10mm, T_DISK=2mm Cu) inside interior half-disc.
Source: A_imposed_θ = B_0 r / 2  →  u_imposed = r·A_θ = B_0 r²/2

Equation (standard axisym A-form, u = r·A_θ):
  ∫(nu/r) grad(u) · grad(v) dx = ∫J_θ · v dx

COMSOL TEAM 28 iteration:
  G_n = ∫_cond (J_n²/σ) × 2πr dr dz
  R_n = 1/G_n
  Solve K · u_n = ∫J_n · v dx → A_θ_n = u_n / r
  A_acc_θ += R_n × A_θ_n
  τ_n = ∫_cond J_n × A_acc_θ × 2πr dr dz
  L_n = R_n × τ_n
  J_{n+1} = J_n - σ × A_acc_θ / L_n
"""
from __future__ import annotations

import json
import sys
import time
from math import pi, sqrt
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "src" / "radia"))

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from netgen.occ import (
    WorkPlane, MoveTo, Vertex, Glue, OCCGeometry, Pnt, Axes, X, Z,
    IdentificationType,
)
from ngsolve import (
    Mesh, H1, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, x, y, dx, IfPos, sqrt as ngs_sqrt, grad,
    Integrate, TaskManager, ngsglobals,
)

from kelvin_source import (
    kelvin_nu_factor_axisym_cf, build_material_cf,
)


# === Physical parameters ====================================================
mu0 = 4 * pi * 1e-7
nu0 = 1 / mu0
SIGMA_CU = 5.8e7
B0 = 1.0
R_DISK = 10e-3   # cylinder radius
T_DISK = 2e-3    # cylinder thickness

# Kelvin parameters
R_K = 50e-3            # Kelvin sphere radius (= inner half-disc radius)
Z_OFFSET = 5 * R_K     # exterior half-disc center z-position

# FE parameters
ORDER = 2
MAXH_INNER = 2e-3
MAXH_KELVIN = 5e-3

N_STAGES = 6


def create_geometry():
    """Build axisym Z-offset Kelvin two-half-disc with cylinder conductor.

    Note: NGSolve axisym convention r=x, z=y. We use 2D (x, y) coordinates.
    """
    # Interior FULL disc at origin
    wp_inner = WorkPlane()
    inner_full = wp_inner.Circle(R_K).Face()
    inner_full.name = "air_inner"

    # Conductor: rectangle in (r, z) = cylinder cross-section
    cond_full = MoveTo(0, -T_DISK / 2).Rectangle(R_DISK, T_DISK).Face()
    cond_full.name = "conductor"
    cond_full.maxh = MAXH_INNER / 4

    # Exterior Kelvin FULL disc at (0, z_offset)
    # Name MUST contain "outer" so build_material_cf applies Kelvin factor
    wp_outer = WorkPlane(Axes((0, Z_OFFSET, 0), n=Z, h=X))
    outer_full = wp_outer.Circle(R_K).Face()
    outer_full.name = "air_outer"

    # GND vertex at exterior center (image of infinity)
    gnd_point = Vertex(Pnt(0, Z_OFFSET, 0))
    gnd_point.name = "GND"

    # Cut to half-disc: keep r >= 0 (x >= 0)
    margin = 0.1
    cutter_inner = MoveTo(-R_K - margin, -R_K - margin).Rectangle(
        R_K + margin, 2 * R_K + 2 * margin).Face()
    cutter_outer = MoveTo(-R_K - margin, Z_OFFSET - R_K - margin).Rectangle(
        R_K + margin, 2 * R_K + 2 * margin).Face()

    inner_half = inner_full - cutter_inner
    cond_half = cond_full - cutter_inner
    outer_half = outer_full - cutter_outer

    inner_air_half = inner_half - cond_half
    inner_air_half.name = "air_inner"

    # Identify inner Kelvin arc (radius R_K from origin)
    kelvin_inner_edges = []
    for edge in inner_air_half.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x**2 + v0.p.y**2)
            d1 = sqrt(v1.p.x**2 + v1.p.y**2)
            is_kelvin = (abs(d0 - R_K) < 1e-3 and abs(d1 - R_K) < 1e-3
                         and cx > 1e-3)
        except Exception:
            is_kelvin = False
        if cx < 1e-3:
            edge.name = "axis"
        elif is_kelvin:
            edge.name = "kelvin_int"
            kelvin_inner_edges.append(edge)
        else:
            edge.name = "default"

    for edge in cond_half.edges:
        cx = edge.center.x
        if cx < 1e-3:
            edge.name = "axis"
        else:
            edge.name = "cond_bnd"

    # Identify outer Kelvin arc (radius R_K from offset center)
    kelvin_outer_edges = []
    for edge in outer_half.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x**2 + (v0.p.y - Z_OFFSET)**2)
            d1 = sqrt(v1.p.x**2 + (v1.p.y - Z_OFFSET)**2)
            is_kelvin = (abs(d0 - R_K) < 1e-3 and abs(d1 - R_K) < 1e-3
                         and cx > 1e-3)
        except Exception:
            is_kelvin = False
        if cx < 1e-3:
            edge.name = "axis_ext"
        elif is_kelvin:
            edge.name = "kelvin_ext"
            kelvin_outer_edges.append(edge)
        else:
            edge.name = "default"

    # Periodic identification (inner ↔ outer kelvin arcs by sign of z)
    for int_edge in kelvin_inner_edges:
        int_y = int_edge.center.y
        for ext_edge in kelvin_outer_edges:
            ext_y = ext_edge.center.y - Z_OFFSET
            if (int_y > 0 and ext_y > 0) or (int_y < 0 and ext_y < 0):
                int_edge.Identify(ext_edge, "kelvin",
                                   IdentificationType.PERIODIC)
                break

    shape = Glue([inner_air_half, cond_half, outer_half, gnd_point])
    geo = OCCGeometry(shape, dim=2)
    return geo


def run_one(maxh_inner, maxh_kelvin, label=""):
    global MAXH_INNER, MAXH_KELVIN
    MAXH_INNER = maxh_inner
    MAXH_KELVIN = maxh_kelvin
    return run_main(label=label)


def run_main(label=""):
    print("=" * 78, flush=True)
    print(f" COMSOL TEAM 28 CLN (axisym + Kelvin) — {label}", flush=True)
    print("=" * 78, flush=True)
    print(f"  Cylinder R_DISK={R_DISK*1000} mm, T_DISK={T_DISK*1000} mm",
          flush=True)
    print(f"  Kelvin R_K={R_K*1000} mm, Z_OFFSET={Z_OFFSET*1000} mm",
          flush=True)
    print(f"  ORDER={ORDER}", flush=True)
    print()

    print(" Building geometry + mesh ...", flush=True)
    t0 = time.time()
    geo = create_geometry()
    ngmesh = geo.GenerateMesh(maxh=MAXH_KELVIN)
    mesh = Mesh(ngmesh)
    print(f"  ne={mesh.ne}, materials={mesh.GetMaterials()}",
          flush=True)
    print(f"  boundaries={mesh.GetBoundaries()}", flush=True)
    print(f"  build mesh: {time.time()-t0:.1f} s", flush=True)

    # FE space (standard H1 + Periodic for Kelvin BC + Dirichlet on axis & GND)
    fes_pre = H1(mesh, order=ORDER, dirichlet="axis|axis_ext",
                  dirichlet_bbnd="GND")
    fes = Periodic(fes_pre)
    print(f"  fes ndof={fes.ndof}", flush=True)
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  free DOFs = {n_free}", flush=True)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    # r-weight for axisym
    r_weight = IfPos(x - 1e-10, x, 1e-10)

    # Kelvin nu factor (only in air_kelvin region)
    nu_kelvin_factor = kelvin_nu_factor_axisym_cf(z_offset=Z_OFFSET, R=R_K)
    nu_cf = build_material_cf(mesh, nu0, nu_kelvin_factor)

    # K matrix: ∫(nu/r) grad(u) · grad(v) dx
    a_form = BilinearForm(fes, symmetric=True)
    a_form += nu_cf / r_weight * grad(u) * grad(v) * dx
    print(" Assembling K ...", flush=True)
    t0 = time.time()
    with TaskManager(): a_form.Assemble()
    print(f"  K assemble: {time.time()-t0:.1f} s", flush=True)

    # M matrix: ∫σ × u × v / r dx (for J²/σ computation in u=r·A_θ form,
    #   J = -σ ∂A_θ/∂t and A_θ = u/r,
    #   J²/σ = σ × A_θ² (= σ × u²/r²),
    #   ∫J²/σ × 2πr dr dz = 2π σ × ∫u²/r dr dz   -- weight is 1/r in axisym!)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)
    m_form = BilinearForm(fes, symmetric=True)
    m_form += sigma_cf / r_weight * u * v * dx
    print(" Assembling M ...", flush=True)
    t0 = time.time()
    with TaskManager(): m_form.Assemble()
    print(f"  M assemble: {time.time()-t0:.1f} s", flush=True)

    # Convert to scipy
    def to_csr(bf_mat, n):
        rs, cs, vs = bf_mat.COO()
        K = sp.csr_matrix(
            (np.asarray(vs), (np.asarray(rs), np.asarray(cs))),
            shape=(n, n))
        return (K + K.T) * 0.5

    K_csr = to_csr(a_form.mat, fes.ndof)
    M_csr = to_csr(m_form.mat, fes.ndof)

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                     dtype=int)
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]

    # Identify cond DOFs (M_red diag != 0)
    M_diag = M_red.diagonal()
    cond_local = np.where(np.abs(M_diag) > 1e-30)[0]
    print(f"  cond DOFs = {len(cond_local)}", flush=True)
    M_cond = M_red[cond_local[:, None], cond_local[None, :]]

    # Add small Tikhonov regularization to K to handle any null mode
    K_reg = K_red + 1e-10 * sp.eye(K_red.shape[0], format="csr")
    K_factor = spla.factorized(K_reg.tocsc())
    M_cond_factor = spla.factorized(M_cond.tocsc())

    # Initial: u_imposed = r × A_θ_imposed = B_0 r²/2 (in conductor)
    u_imposed_cf = B0 * x**2 / 2  # = B_0 x²/2 in NGSolve coords (x=r)
    # Project u_imposed to fes via L2:
    # b vector with σ × A_θ_imposed × v / r weighting
    # b_i = ∫σ × (u_imposed/r) × v_i dx = ∫σ × u_imposed × v_i / r dx
    b_form = LinearForm(fes)
    b_form += sigma_cf * u_imposed_cf / r_weight * v * dx
    with TaskManager(): b_form.Assemble()
    b_full = np.array(b_form.vec)
    b_red = b_full[free]
    b_cond = b_red[cond_local]

    # Initial u_dof_cond from σ-weighted projection
    u_dof_cond = M_cond_factor(b_cond)

    A_acc_full = np.zeros(len(free))

    Rn, Ln, tau_us = [], [], []
    diag = []

    for n in range(N_STAGES):
        t_stage = time.time()

        # All in NGSolve 2D coordinates (no 2π factor — they cancel in τ = L/R)
        G_n = float(u_dof_cond @ (M_cond @ u_dof_cond))
        if G_n < 1e-50 or not np.isfinite(G_n):
            print(f"  Stage {n}: G_n = {G_n}, breaking", flush=True)
            break
        R_n = 1.0 / G_n
        Rn.append(R_n)

        # Loading b_n: M_cond × u_dof_cond on cond, 0 elsewhere
        b_n_full = np.zeros(len(free))
        b_n_full[cond_local] = M_cond @ u_dof_cond

        # Solve K × w_n = b_n
        w_n = K_factor(b_n_full)

        # Accumulate A_acc_full
        A_acc_full += R_n * w_n

        # τ_n = b_n_cond · A_acc_cond (Galerkin pairing)
        A_acc_cond = A_acc_full[cond_local]
        tau_n = float((M_cond @ u_dof_cond) @ A_acc_cond)
        tau_n_us_val = tau_n * 1e6
        L_n = R_n * tau_n
        Ln.append(L_n)
        tau_us.append(tau_n_us_val)

        # Update u_dof_cond -= A_acc_cond / L_n
        u_dof_cond = u_dof_cond - A_acc_cond / L_n

        diag.append({
            "stage": n,
            "G_n": G_n,
            "R_n": R_n,
            "tau_n_us": tau_n_us_val,
            "L_n": L_n,
            "stage_time_s": time.time() - t_stage,
        })
        print(f"  Stage {n}: G={G_n:.4e}, R={R_n:.4e},"
              f" τ={tau_n_us_val:.4f} μs, L={L_n:.4e}"
              f"  ({time.time()-t_stage:.2f} s)",
              flush=True)

    print("\n" + "=" * 78, flush=True)
    print(" Summary (axisym Kelvin, COMSOL TEAM 28)", flush=True)
    print("=" * 78, flush=True)
    print("\n  Reference (v23 axisym Hiruma 3-term):", flush=True)
    v23_ref = [218.6222, 78.2068, 39.5928, 23.1196, 16.0682, 13.3002]
    for k in range(min(N_STAGES, len(v23_ref))):
        match = ""
        if k < len(tau_us):
            gap = (tau_us[k] - v23_ref[k]) / v23_ref[k] * 100
            match = (f"  τ_{k} (Kelvin) = {tau_us[k]:.4f} μs"
                     f"  gap = {gap:+.1f}%")
        print(f"    k={k}: v23 τ = {v23_ref[k]:.4f} us{match}",
              flush=True)

    out = {"Rn": Rn, "Ln": Ln, "tau_n_us": tau_us, "stages": diag}
    out_path = (Path(__file__).parent
                 / "2026-05-10-cln_team28_axisym_kelvin_results.json")
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    main()
