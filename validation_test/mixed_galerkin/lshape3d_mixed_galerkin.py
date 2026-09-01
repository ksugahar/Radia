"""
3D L-shape mixed Galerkin: rank-1 FEM Krylov bulk + rank-1 distance-function
surface envelope, coupled via Schur complement.

L-shape regions:
  horiz = [0, 2L] x [0, L] x [0, L]
  vert  = [0,  L] x [0, L] x [L, 2L]

Outer boundary = 8 face pieces (L-shape has one concave dihedral edge at x=L,
z=L, y in [0, L]).

Bulk basis: q1 = K0^{-1} b  (FEM Krylov at s=0, single mode)
            K0 = -Laplace,  b = integrate v dV
Surface basis: psi(r, s) = cosh(gamma (d_max - d(r))) / cosh(gamma d_max) - 1
            d(r) = signed distance to nearest outer face (L-shape, 8 pieces)
            d_max = L/2 (half of smallest body thickness)

Mixed system per s:
  [K_bb(s)  K_bs(s)] [xi_b]   [-sMS b_b]
  [K_bs    K_ss(s)] [xi_s] = [-sMS b_s]

Y(s) = sigma V_L + sigma (xi_b b_b + xi_s b_s)

Compared against full FEM Y(s) ground truth on the same mesh.
"""
import cmath
import math
import time

import numpy as np

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, Integrate, TaskManager,
    SetNumThreads, IfPos, cosh, sinh, cos, sin, InnerProduct,
    x as X, y as Y, z as Z,
)
from netgen.occ import OCCGeometry, Box, Pnt, Glue

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA

V_L = 3 * L**3
S_L = 14 * L**2
Y_DC = SIGMA * V_L
K_SIBC = S_L * math.sqrt(SIGMA / MU)


def build_lshape(maxh):
    horiz = Box(Pnt(0, 0, 0), Pnt(2*L, L, L))
    vert  = Box(Pnt(0, 0, L), Pnt(L, L, 2*L))
    shape = Glue([horiz, vert])
    shape.mat("Lshape")
    for f in shape.faces:
        f.name = "outer"
    return Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))


def _min_cf(*args):
    """Min over many real CFs via repeated IfPos."""
    r = args[0]
    for a in args[1:]:
        r = IfPos(r - a, a, r)
    return r


def distance_to_boundary_L():
    """Signed distance from interior point to nearest outer face piece.

    L-shape = horiz ∪ vert (Glue). Interior face at z=L, x ∈ [0, L], y ∈ [0, L]
    is NOT outer.

    Region split (z < L -> horiz, z > L -> vert):
      In horiz: 6 face pieces x=0, x=2L, y=0, y=L, z=0, and z=L for x > L.
      In vert:  5 face pieces x=0, x=L (concave outer), y=0, y=L, z=2L.
    """
    BIG = 1e10  # mark face as "not present in this region"

    d_x0 = X            # always relevant (full L-shape face)
    d_y0 = Y            # always relevant
    d_yL = L - Y        # always relevant

    # Horiz region (z < L)
    d_x2L_h = 2*L - X
    d_z0    = Z
    d_zL_overhang = IfPos(X - L, L - Z, BIG)  # only x > L has z=L outer face
    d_horiz = _min_cf(d_x0, d_x2L_h, d_y0, d_yL, d_z0, d_zL_overhang)

    # Vert region (z > L)
    d_xL_v = L - X      # concave outer face of vert
    d_z2L  = 2*L - Z
    d_vert = _min_cf(d_x0, d_xL_v, d_y0, d_yL, d_z2L)

    return IfPos(L - Z, d_horiz, d_vert)


def _planar_1d_factor(xi, half_thick, gamma_val):
    """Complex CF: cosh((xi - half_thick) gamma)/cosh(half_thick gamma) - 1.

    Vanishes at xi=0 and xi=2*half_thick (boundary of 1D slab).
    For xi inside [0, 2*half_thick], gives planar SIBC profile.

    Handles complex gamma by splitting cosh into Re/Im parts.
    """
    gre = gamma_val.real
    gim = gamma_val.imag

    a = gre * (xi - half_thick)
    b = gim * (xi - half_thick)
    num_re = cosh(a) * cos(b)
    num_im = sinh(a) * sin(b)

    den = cmath.cosh(gamma_val * half_thick)
    den_sq = den.real**2 + den.imag**2
    re = (num_re * den.real + num_im * den.imag) / den_sq - 1.0
    im = (num_im * den.real - num_re * den.imag) / den_sq

    return CoefficientFunction(re + 1j * im)


def psi_cf_boundingbox(gamma_val):
    """Bounding-box tensor SIBC envelope for L-shape (legacy, BC mismatch).

    Bounding box = [0, 2L] x [0, L] x [0, 2L]
    Satisfies the BC on 6 out of 8 face pieces. Does NOT satisfy BC on
    the 2 concave-step faces (x=L for z>L, z=L for x>L) - that BC
    violation is the L-shape-specific error vs the cube case.
    """
    fx = _planar_1d_factor(X, L,     gamma_val)
    fy = _planar_1d_factor(Y, L / 2, gamma_val)
    fz = _planar_1d_factor(Z, L,     gamma_val)
    return fx * fy * fz


def psi_horiz_only(gamma_val):
    """Per-region tensor SIBC envelope restricted to horiz region (z < L).

    horiz = [0, 2L] x [0, L] x [0, L]; half-thicknesses (L, L/2, L/2).
    Vanishes at x=0, x=2L (horiz outer), y=0, y=L, z=0, AND z=L (= top
    of horiz; this is interior for x<L but boundary for x>L overhang).

    Multiplied by IfPos(L-Z, 1, 0) to be zero in vert region. Continuous
    across z=L interface (both sides are 0 there).
    """
    fx = _planar_1d_factor(X, L,         gamma_val)
    fy = _planar_1d_factor(Y, L / 2,     gamma_val)
    fz = _planar_1d_factor(Z, L / 2,     gamma_val)  # axis [0, L], vanishes at 0, L
    return IfPos(L - Z, fx * fy * fz, CoefficientFunction(0.0 + 0.0j))


def psi_vert_only(gamma_val):
    """Per-region tensor SIBC envelope restricted to vert region (z > L).

    vert = [0, L] x [0, L] x [L, 2L]; half-thicknesses (L/2, L/2, L/2).
    Vanishes at x=0, x=L (concave outer), y=0, y=L, z=2L, AND z=L (= bottom
    of vert; interior for x<L; matches horiz at interface).

    Shifted z-axis: use zeta = z - L for the z-factor.
    """
    fx = _planar_1d_factor(X,         L / 2, gamma_val)
    fy = _planar_1d_factor(Y,         L / 2, gamma_val)
    fz = _planar_1d_factor(Z - L,     L / 2, gamma_val)  # axis zeta in [0, L]
    return IfPos(Z - L, fx * fy * fz, CoefficientFunction(0.0 + 0.0j))


def psi_per_region_sum(gamma_val):
    """Sum of horiz and vert per-region envelopes — single 1-DOF basis.

    Continuous on L-shape, vanishes on ALL 8 outer face pieces
    (including the 2 concave-step faces). Gradient has a jump at z=L
    interface for x<L (interior interface) but |grad psi|^2 is L^2
    integrable so the variational form is well-posed.
    """
    return psi_horiz_only(gamma_val) + psi_vert_only(gamma_val)


def solve_fem_ground_truth(mesh, f, order):
    s = 1j * 2 * math.pi * f
    sMS = s * MS
    fes = H1(mesh, order=order, complex=True, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    a += sMS * u * v * dx
    F = LinearForm(fes)
    F += -sMS * v * dx
    a.Assemble()
    F.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * F.vec
    u_avg = Integrate(gfu, mesh) / V_L
    return Y_DC * (1 + u_avg)


def precompute_bulk_krylov_modes(mesh, order, N_bulk):
    """FEM Krylov-at-s=0 modes:
      q_1 = K0^{-1} b   (b = integrate v dV)
      q_n = K0^{-1} (MS * q_{n-1})  for n >= 2

    Returns: list of N_bulk GridFunctions in real FE space.
    """
    fes_r = H1(mesh, order=order, dirichlet="outer")
    u, v = fes_r.TnT()
    a = BilinearForm(fes_r, symmetric=True)
    a += grad(u) * grad(v) * dx
    a.Assemble()
    inv = a.mat.Inverse(freedofs=fes_r.FreeDofs(), inverse="sparsecholesky")
    F = LinearForm(fes_r)
    F += v * dx
    F.Assemble()

    # Mass matrix for MS * q step
    m = BilinearForm(fes_r, symmetric=True)
    m += MS * u * v * dx
    m.Assemble()

    qs = []
    q1 = GridFunction(fes_r)
    q1.vec.data = inv * F.vec
    qs.append(q1)
    for n in range(2, N_bulk + 1):
        qn = GridFunction(fes_r)
        # qn = K0^{-1} (M_MS * q_{n-1})
        tmp = qs[-1].vec.CreateVector()
        tmp.data = m.mat * qs[-1].vec
        qn.vec.data = inv * tmp
        qs.append(qn)
    return qs

# Keep old single-mode name for backward compat
def precompute_bulk_krylov_mode(mesh, order):
    return precompute_bulk_krylov_modes(mesh, order, 1)[0]


def solve_mixed_rankN_plus_surface(mesh, f, qs, surface_mode="none"):
    """Rank-N bulk (qs from precompute) + optional surface envelope.

    surface_mode:
      "none"        - bulk only
      "boundingbox" - 1-DOF, bounding-box tensor (BC mismatch at concave step)
      "perregion"   - 1-DOF, sum of horiz + vert per-region tensors (per-face PoU)
      "split"       - 2-DOF, horiz + vert as separate basis functions
    """
    s = 1j * 2 * math.pi * f
    sMS = s * MS
    gamma_val = cmath.sqrt(s * MS)

    N = len(qs)

    # Build list of surface CF basis functions
    if surface_mode == "none":
        psi_list = []
    elif surface_mode == "boundingbox":
        psi_list = [psi_cf_boundingbox(gamma_val)]
    elif surface_mode == "perregion":
        psi_list = [psi_per_region_sum(gamma_val)]
    elif surface_mode == "split":
        psi_list = [psi_horiz_only(gamma_val), psi_vert_only(gamma_val)]
    else:
        raise ValueError(f"Unknown surface_mode: {surface_mode}")

    n_surf = len(psi_list)
    M = N + n_surf

    def Igrad_grad(a, b):
        return Integrate(InnerProduct(a, b), mesh)

    def Imass(a, b):
        return Integrate(a * b, mesh)

    K = np.zeros((M, M), dtype=complex)
    b_vec = np.zeros(M, dtype=complex)

    grads_q = [grad(q) for q in qs]

    # Bulk-bulk
    for i in range(N):
        for j in range(i, N):
            val = Igrad_grad(grads_q[i], grads_q[j]) + sMS * Imass(qs[i], qs[j])
            K[i, j] = val
            K[j, i] = val
        b_vec[i] = Imass(qs[i], CoefficientFunction(1.0))

    # Surface gradients via Diff(X,Y,Z) on the complex CF
    grads_psi = []
    for psi in psi_list:
        grads_psi.append(CoefficientFunction((psi.Diff(X), psi.Diff(Y), psi.Diff(Z))))

    # Bulk-surface and surface-surface
    for k in range(n_surf):
        psi = psi_list[k]
        gpsi = grads_psi[k]
        for i in range(N):
            val = Igrad_grad(grads_q[i], gpsi) + sMS * Imass(qs[i], psi)
            K[i, N + k] = val
            K[N + k, i] = val
        for l in range(k, n_surf):
            val = Igrad_grad(gpsi, grads_psi[l]) + sMS * Imass(psi, psi_list[l])
            K[N + k, N + l] = val
            K[N + l, N + k] = val
        b_vec[N + k] = Imass(psi, CoefficientFunction(1.0))

    rhs_vec = -sMS * b_vec
    xi = np.linalg.solve(K, rhs_vec)
    u_avg = (xi @ b_vec) / V_L
    return Y_DC * (1 + u_avg)


# Backwards-compat alias used by earlier main()
def solve_mixed_rankN_plus_1(mesh, f, qs, with_surface=True):
    return solve_mixed_rankN_plus_surface(
        mesh, f, qs, "boundingbox" if with_surface else "none")


# Keep simpler entry point used by main()
def solve_mixed_rank1_plus_1(mesh, f, order, q1_real, d_cf, d_max):
    Y_mix = solve_mixed_rankN_plus_1(mesh, f, [q1_real], with_surface=True)
    Y_b1  = solve_mixed_rankN_plus_1(mesh, f, [q1_real], with_surface=False)
    return Y_mix, Y_b1


def main():
    SetNumThreads(8)
    print("=== 3D L-shape rank-1 + 1 mixed Galerkin ===")
    print(f"L = {L*1e3} mm, V_L = {V_L:.4e} m^3, S_L = {S_L:.4e} m^2")
    print(f"Y_DC = {Y_DC:.4e} S, K_SIBC^L = {K_SIBC:.4e}")
    print()

    maxh = L / 8
    order = 3
    d_max = L / 2  # half thickness of bulkiest direction

    print(f"Mesh: maxh = L/8, P{order}")
    mesh = build_lshape(maxh)
    print(f"  ne = {mesh.ne}")

    print(f"d_max = L/2 = {d_max*1e3:.3f} mm (envelope decay length scale)")
    print()

    # Precompute FEM Krylov bulk modes at s=0 (up to rank N_max)
    N_max = 10
    print(f"Precomputing q1..q{N_max} = FEM Krylov modes at s=0...")
    with TaskManager():
        t0 = time.time()
        qs = precompute_bulk_krylov_modes(mesh, order, N_max)
        print(f"  {N_max} modes done in {time.time()-t0:.2f}s")
    print()

    N_RANKS = [3, 5, 7, 10]
    F_LIST = [1e3, 1e4, 1e5, 1e6]
    SURFACE_MODES = ["none", "boundingbox"]
    MODE_LABEL = {
        "none":        "bulk-only",
        "boundingbox": "+bbox 1DOF",
    }

    # 1. Ground truth FEM Y(s)
    print("=== Step A: FEM ground truth ===")
    Y_fem_list = {}
    for f in F_LIST:
        with TaskManager():
            Y_fem = solve_fem_ground_truth(mesh, f, order)
        Y_fem_list[f] = Y_fem
        s = 1j * 2 * math.pi * f
        gL = abs(cmath.sqrt(s * MS) * L)
        print(f"  f={f:.0e}  |gL|={gL:6.2f}  |Y_FEM|={abs(Y_fem):.6e}")
    print()

    # 2. Sweep: for each rank N, compare 4 surface modes at all frequencies
    for N in N_RANKS:
        print(f"=== Step B: rank-N={N} bulk + various surface envelopes (% err vs FEM) ===")
        header = f"  {'f (Hz)':>10}  {'|gL|':>6}  " + "  ".join(
            f"{MODE_LABEL[m]:>13}" for m in SURFACE_MODES)
        print(header)
        for f in F_LIST:
            s = 1j * 2 * math.pi * f
            gL = abs(cmath.sqrt(s * MS) * L)
            line = f"  {f:10.2e}  {gL:6.2f}"
            for mode in SURFACE_MODES:
                with TaskManager():
                    Y = solve_mixed_rankN_plus_surface(mesh, f, qs[:N], mode)
                err = (abs(Y) - abs(Y_fem_list[f])) / abs(Y_fem_list[f]) * 100
                line += f"  {err:+11.3f}%"
            print(line)
        print()

    print("Notes:")
    print("- bulk-only is straight CLN (rank-N Krylov-at-0), Q(s) rational, integer tail.")
    print("- +bbox 1DOF: bounding-box tensor SIBC (BC mismatch at 2 concave step faces).")
    print("- +region 1DOF: sum of horiz + vert per-region tensors, vanishes on ALL 8 faces.")
    print("- +split 2DOF: horiz and vert per-region tensors as 2 separate basis functions.")
    print("- Wall band ~ f=1e4 (|gL| ~ 10); deep skin ~ f=1e6 (|gL| ~ 100).")


if __name__ == "__main__":
    main()
