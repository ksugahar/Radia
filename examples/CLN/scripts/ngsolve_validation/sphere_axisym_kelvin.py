"""sphere_axisym_kelvin.py — Sphere in uniform B_z field, axisymmetric NGSolve.

Assembles the sphere operators and reports the Stoll analytical spectrum
for a Cu sphere of radius a in uniform B_0 = 1 T applied (z direction).
Historical Hiruma and Kameari method comparisons were retired; these
attributions are retained without the field-recurrence implementations.

Stoll analytical Foster pole spectrum:
  τ_n = μ_0 σ a² / (n² π²),  n = 1, 2, 3, ...

For Cu (σ = 5.8e7 S/m), a = 10 mm:
  τ_1 = 738.7 μs (leading)
  τ_2 = 184.7 μs
  τ_3 =  82.1 μs
  τ_4 =  46.2 μs
  τ_5 =  29.6 μs
  τ_6 =  20.5 μs

Geometry (axisym, NGSolve coords r=x, z=y):
  Conductor: half-disc r² + z² ≤ a², r ≥ 0
  Inner air: rest of inner half-disc r² + z² ≤ R_K²
  Outer Kelvin half-disc at (0, z_offset = 5*R_K)
  Periodic BC between kelvin_int and kelvin_ext arcs
"""
from __future__ import annotations

import json
import sys
import time
from math import pi, sqrt
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from netgen.occ import (
    WorkPlane, MoveTo, Vertex, Glue, OCCGeometry, Pnt, Axes, X, Z,
    IdentificationType,
)
from ngsolve import (
    Mesh, H1, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, x, y, dx, IfPos, grad,
    Integrate, TaskManager, ngsglobals,
)

from kelvin_source import (
    kelvin_nu_factor_axisym_cf, build_material_cf,
)


mu0 = 4 * pi * 1e-7
nu0 = 1 / mu0
SIGMA_CU = 5.8e7
B0 = 1.0
SPHERE_R = 10e-3   # sphere radius

R_K = 50e-3
Z_OFFSET = 5 * R_K

ORDER = 3   # P3 for curved boundary
MAXH_INNER = 0.5e-3   # finer mesh near sphere
MAXH_KELVIN = 3e-3
N_STAGES = 30


def stoll_tau_us(n, a, sigma):
    """Stoll analytical: τ_n = μ_0 σ a² / (n² π²) in seconds."""
    return mu0 * sigma * a**2 / (n**2 * pi**2) * 1e6


def create_geometry():
    """Half-disc sphere conductor + inner air + Kelvin two-half-disc."""
    # Inner full disc (will be cut to half)
    wp_inner = WorkPlane()
    inner_full = wp_inner.Circle(R_K).Face()
    inner_full.name = "air_inner"

    # Conductor: half-disc sphere of radius SPHERE_R at origin
    wp_sphere = WorkPlane()
    sphere_full = wp_sphere.Circle(SPHERE_R).Face()
    sphere_full.name = "conductor"
    sphere_full.maxh = MAXH_INNER

    # Outer Kelvin disc
    wp_outer = WorkPlane(Axes((0, Z_OFFSET, 0), n=Z, h=X))
    outer_full = wp_outer.Circle(R_K).Face()
    outer_full.name = "air_outer"

    gnd_point = Vertex(Pnt(0, Z_OFFSET, 0))
    gnd_point.name = "GND"

    margin = 0.1
    cutter_inner = MoveTo(-R_K - margin, -R_K - margin).Rectangle(
        R_K + margin, 2 * R_K + 2 * margin).Face()
    cutter_outer = MoveTo(-R_K - margin, Z_OFFSET - R_K - margin).Rectangle(
        R_K + margin, 2 * R_K + 2 * margin).Face()

    inner_half = inner_full - cutter_inner
    sphere_half = sphere_full - cutter_inner
    outer_half = outer_full - cutter_outer

    inner_air_half = inner_half - sphere_half
    inner_air_half.name = "air_inner"

    # Kelvin arc identification
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

    for edge in sphere_half.edges:
        cx = edge.center.x
        if cx < 1e-3:
            edge.name = "axis"
        else:
            edge.name = "sphere_bnd"

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

    for int_edge in kelvin_inner_edges:
        int_y = int_edge.center.y
        for ext_edge in kelvin_outer_edges:
            ext_y = ext_edge.center.y - Z_OFFSET
            if (int_y > 0 and ext_y > 0) or (int_y < 0 and ext_y < 0):
                int_edge.Identify(ext_edge, "kelvin",
                                   IdentificationType.PERIODIC)
                break

    shape = Glue([inner_air_half, sphere_half, outer_half, gnd_point])
    return OCCGeometry(shape, dim=2)


def to_csr(bf_mat, n):
    rs, cs, vs = bf_mat.COO()
    K = sp.csr_matrix(
        (np.asarray(vs), (np.asarray(rs), np.asarray(cs))),
        shape=(n, n))
    return (K + K.T) * 0.5






def main():
    print("=" * 78, flush=True)
    print(" Sphere in uniform B_z (axisym + Kelvin)", flush=True)
    print("=" * 78, flush=True)
    print(f"  Sphere radius a = {SPHERE_R*1000} mm",
          flush=True)
    print(f"  σ = {SIGMA_CU:.2e} S/m, B_0 = {B0} T",
          flush=True)
    print(f"  Kelvin R_K = {R_K*1000} mm, Z_OFFSET = {Z_OFFSET*1000} mm",
          flush=True)
    print()

    # Stoll analytical reference
    print(" Stoll Bessel analytical Foster pole spectrum:", flush=True)
    stoll_taus = []
    for n in range(1, N_STAGES + 1):
        t = stoll_tau_us(n, SPHERE_R, SIGMA_CU)
        stoll_taus.append(t)
        print(f"   τ_{n} = μ_0 σ a² / ({n}² π²) = {t:.4f} μs",
              flush=True)
    print()

    # Build geometry & mesh
    print(" Building geometry + mesh ...", flush=True)
    t0 = time.time()
    geo = create_geometry()
    ngmesh = geo.GenerateMesh(maxh=MAXH_KELVIN)
    mesh = Mesh(ngmesh)
    # Curve mesh to match sphere boundary at high order
    mesh.Curve(ORDER)
    print(f"  ne={mesh.ne}, materials={mesh.GetMaterials()}",
          flush=True)
    print(f"  Curved with order {ORDER}", flush=True)
    print(f"  build mesh: {time.time()-t0:.1f} s", flush=True)

    # FE space (H1 P2 + Periodic for Kelvin)
    fes_pre = H1(mesh, order=ORDER, dirichlet="axis|axis_ext",
                  dirichlet_bbnd="GND")
    fes = Periodic(fes_pre)
    n_free = sum(1 for f in fes.FreeDofs() if f)
    print(f"  fes ndof={fes.ndof}, free={n_free}", flush=True)

    u, v = fes.TnT()
    r_weight = IfPos(x - 1e-10, x, 1e-10)
    nu_kelvin_factor = kelvin_nu_factor_axisym_cf(z_offset=Z_OFFSET, R=R_K)
    nu_cf = build_material_cf(mesh, nu0, nu_kelvin_factor)
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)

    a_form = BilinearForm(fes, symmetric=True)
    a_form += nu_cf / r_weight * grad(u) * grad(v) * dx
    print(" Assembling K ...", flush=True)
    with TaskManager(): a_form.Assemble()

    m_form = BilinearForm(fes, symmetric=True)
    m_form += sigma_cf / r_weight * u * v * dx
    print(" Assembling M ...", flush=True)
    with TaskManager(): m_form.Assemble()

    A_imposed = B0 * x / 2
    b_form = LinearForm(fes)
    b_form += sigma_cf * A_imposed * v * dx
    with TaskManager(): b_form.Assemble()
    b_full = np.array(b_form.vec)

    K_csr = to_csr(a_form.mat, fes.ndof)
    M_csr = to_csr(m_form.mat, fes.ndof)

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                     dtype=int)
    K_red = K_csr[free[:, None], free[None, :]]
    M_red = M_csr[free[:, None], free[None, :]]
    b_red = b_full[free]

    M_diag = M_red.diagonal()
    cond_local = np.where(np.abs(M_diag) > 1e-30)[0]
    M_cond = M_red[cond_local[:, None], cond_local[None, :]]

    print(f"  free DOFs={len(free)}, cond DOFs={len(cond_local)}",
          flush=True)

    # Independent assembly and analytical spectrum only.
    out = {"sphere_R_mm": SPHERE_R*1000, "sigma": SIGMA_CU, "B0": B0,
           "stoll_tau_us": stoll_taus,
           "mesh": {"ne": mesh.ne, "ndof": fes.ndof,
                    "n_cond_dof": len(cond_local)},
           "operators": {"K_nnz": K_red.nnz, "M_nnz": M_red.nnz},
           "status": "assembled; no reduced model extracted"}
    out_path = Path(__file__).parent / "sphere_axisym_kelvin_assembly.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Results: {out_path}", flush=True)


if __name__ == "__main__":
    main()
