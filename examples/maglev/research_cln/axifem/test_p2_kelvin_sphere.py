"""Phase B3 axifem + Kelvin transformation sphere benchmark.

Compares Cu sphere leading time-constant tau_1 vs Stoll analytical
(mu_0 sigma R^2 / pi^2 = 738.48 us at R = 10 mm) using:

  - canonical z-offset Kelvin geometry (per radia-mcp kelvin knowledge
    base; same pattern as examples/CLN/scripts/ngsolve_validation/
    sphere_axisym_kelvin.py which uses NGSolve H1 + SymbolicBFI)
  - axifem P2 triangle FE (AxiHenrotteStiffnessBFI + Periodic)
  - mu_eff = mu_0 * (R_K / rho')^2 in the "kelvin" material region
    via radia.kelvin_source.kelvin_mu_factor_axisym_cf

KNOWN axifem caveat (logged for Phase B3.2): AxiHenrotteStiffnessBFI
samples mu at the element CENTROID (axi_henrotte_integrators.cpp:710),
so the (R_K/rho')^2 modulation -- which diverges at rho'=0 (=GND vertex)
-- has discretisation error proportional to the cell variance of
(R_K/rho')^2. To upper-bound that error we (a) put the Kelvin center
at z_offset = 5*R_K (i.e. far from the conductor), (b) refine the
"kelvin" face with maxh_kelvin so cells near the GND vertex are small.
If even refined results stay > a few percent off Stoll, the resolution
will need the integrator API extended to per-quadpt mu evaluation
(Phase B3.2 work item).

Mesh strategy: conductor (sphere half-disc, R_SPHERE) sits in
"air_inner" half-disc of radius R_K centered at origin. add_kelvin_
2d_axisym appends the "kelvin" half-disc at z_offset = 5*R_K, names
"axis"/"axis_ext"/"GND" boundaries, and Periodic-identifies the two
arcs.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")

from math import pi, sqrt as msqrt
import time
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

from netgen.occ import (
    WorkPlane, MoveTo, Vertex, Glue, OCCGeometry, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, Periodic, TaskManager, ngsglobals,
)
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)
from radia.panels.add_kelvin import add_kelvin_2d_axisym
from radia.kelvin_source import kelvin_mu_factor_axisym_cf, build_material_cf

R_SPHERE = 10.0e-3
R_K = 50.0e-3
SIGMA_CU = 5.8e7
MU0 = 4 * pi * 1e-7
STOLL_TAU_1 = MU0 * SIGMA_CU * R_SPHERE ** 2 / (pi ** 2)


def build_interior_half_disk(maxh_cond=1.0e-3):
    """Build interior half-disc (conductor + air_inner) on x >= 0.

    Outer arc named "kelvin_int", axis "axis", sphere boundary "sphere_bnd".
    Conductor face material is "conductor"; surrounding annulus is
    "air_inner".
    """
    # Full inner disc (kelvin radius) and full sphere disc (both centered
    # at origin). Cut with rectangle at x < 0 to get half-discs.
    wp_inner = WorkPlane()
    inner_full = wp_inner.Circle(R_K).Face()

    wp_sphere = WorkPlane()
    sphere_full = wp_sphere.Circle(R_SPHERE).Face()
    sphere_full.name = "conductor"
    sphere_full.maxh = maxh_cond

    # Half-plane cutter (x < 0)
    margin = 0.1
    cutter = MoveTo(-R_K - margin, -R_K - margin).Rectangle(
        R_K + margin, 2 * R_K + 2 * margin).Face()

    inner_half = inner_full - cutter
    sphere_half = sphere_full - cutter

    # air_inner = inner_half minus sphere_half
    air_inner = inner_half - sphere_half
    air_inner.name = "air_inner"

    # Name conductor edges (axis or sphere boundary)
    for edge in sphere_half.edges:
        cx = edge.center.x
        edge.name = "axis" if cx < 1e-4 else "sphere_bnd"

    # Name air_inner edges: x=0 -> axis, outer arc -> kelvin_int
    for edge in air_inner.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = msqrt(v0.p.x ** 2 + v0.p.y ** 2)
            d1 = msqrt(v1.p.x ** 2 + v1.p.y ** 2)
            is_outer = (abs(d0 - R_K) < 0.01 * R_K and
                        abs(d1 - R_K) < 0.01 * R_K and cx > 0.01 * R_K)
        except Exception:
            is_outer = False
        if cx < 1e-4:
            edge.name = "axis"
        elif is_outer:
            edge.name = "kelvin_int"
        else:
            edge.name = "default"

    # Glue conductor + annulus into single interior shape
    interior = Glue([air_inner, sphere_half])
    return interior


def to_csr(mat, n):
    rs, cs, vs = mat.COO()
    return sp.coo_matrix(
        (np.array(vs), (np.array(rs), np.array(cs))),
        shape=(n, n),
    ).toarray()


def run_one(maxh, maxh_kelvin, curve_order=0):
    ngsglobals.msg_level = 0

    interior = build_interior_half_disk(maxh_cond=maxh)
    shape, info = add_kelvin_2d_axisym(
        interior, R=R_K, z_offset=5 * R_K, maxh_kelvin=maxh_kelvin)
    geo = OCCGeometry(shape, dim=2)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh_kelvin))
    if curve_order >= 2:
        mesh.Curve(curve_order)

    materials = mesh.GetMaterials()
    print(f"  Materials: {materials}, boundaries: {mesh.GetBoundaries()}")

    # Wrap H1Henrotte with Periodic to couple kelvin_int <-> kelvin_ext DOFs.
    # AxiHenrotteFESpace subclasses FESpace, so ngsolve.Periodic works on it
    # the same way as on a stock H1 space.
    fes_pre = H1Henrotte(
        mesh, order=2,
        dirichlet="axis|axis_ext",
        dirichlet_bbnd="GND",
    )
    fes = Periodic(fes_pre)

    mu_factor = kelvin_mu_factor_axisym_cf(z_offset=5 * R_K, R=R_K)
    mu_cf = build_material_cf(
        mesh, MU0, mu_factor, outer_keyword="kelvin")
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)

    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager():
        a.Assemble()
    m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
    m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager():
        m_bf.Assemble()

    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]],
                    dtype=int)
    K = to_csr(a.mat, fes.ndof)[np.ix_(free, free)]
    M = to_csr(m_bf.mat, fes.ndof)[np.ix_(free, free)]
    K = 0.5 * (K + K.T)
    M = 0.5 * (M + M.T)

    diag_M = np.diag(M)
    diag_K = np.diag(K)
    cond_mask = diag_M > 1e-30 * (np.max(diag_M) + 1e-30)
    nonzero_K = diag_K > 1e-30 * (np.max(diag_K) + 1e-30)
    air_mask = (~cond_mask) & nonzero_K
    cond_idx = np.where(cond_mask)[0]
    air_idx = np.where(air_mask)[0]

    if len(air_idx) > 0:
        K_cc = K[np.ix_(cond_idx, cond_idx)]
        K_ca = K[np.ix_(cond_idx, air_idx)]
        K_aa = K[np.ix_(air_idx, air_idx)]
        X = sla.solve(K_aa, K_ca.T, assume_a="sym")
        S = K_cc - K_ca @ X
    else:
        S = K[np.ix_(cond_idx, cond_idx)]
    M_cc = M[np.ix_(cond_idx, cond_idx)]
    S = 0.5 * (S + S.T)
    M_cc = 0.5 * (M_cc + M_cc.T)

    eigvals = sla.eigvalsh(
        S, M_cc, subset_by_index=[0, min(3, S.shape[0] - 1)])
    pos = eigvals[eigvals > 0]
    tau = 1.0 / pos[0] if len(pos) > 0 else np.nan
    gap = (tau - STOLL_TAU_1) / STOLL_TAU_1 * 100
    return mesh.ne, fes.ndof, tau, gap


def main():
    print("=" * 70)
    print(f"Phase B3 axifem + Kelvin sphere benchmark vs Stoll "
          f"{STOLL_TAU_1*1e6:.4f} us")
    print(f"  Cu sphere R={R_SPHERE*1e3} mm, Kelvin R_K={R_K*1e3} mm, "
          f"z_offset={5*R_K*1e3} mm")
    print(f"  z-offset Kelvin: nu_kelvin = nu_0 (rho'/R_K)^2 (Nagamine "
          f"CEFC 2026 canonical)")
    print(f"  axifem AxiHenrotteStiffnessBFI: mu sampled at element "
          f"CENTROID")
    print("=" * 70)

    levels = [
        # (maxh_cond, maxh_kelvin)
        (2.0e-3, 8.0e-3),
        (1.0e-3, 5.0e-3),
        (0.7e-3, 3.0e-3),
        (0.5e-3, 2.0e-3),
    ]
    for maxh_c, maxh_k in levels:
        for cv in [0, 2]:
            tag = f"Curve({cv})"
            print(f"\nmaxh_cond={maxh_c:.1e}  maxh_kelvin={maxh_k:.1e}  {tag}")
            try:
                ne, ndof, tau, gap = run_one(maxh_c, maxh_k, curve_order=cv)
                print(f"  ne={ne}  ndof={ndof}  tau_1={tau*1e6:.4f} us  "
                      f"gap={gap:+.3f}%")
            except Exception as e:
                print(f"  FAIL: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
