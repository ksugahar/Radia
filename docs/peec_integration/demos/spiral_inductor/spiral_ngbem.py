"""
spiral_ngbem.py

BEM (Boundary Element Method) model for a 14-turn rectangular spiral inductor.

Each turn is modeled as a closed rectangular ring surface mesh.
Self-inductance per turn: ngsbem LaplaceSL + Hodge decomposition.
Mutual inductance between turns: Grover analytical formula (parallel filaments).

Comparison target:
  - PEEC (spiral_peec.py): L_DC = 13.91 uH
  - Measurement: L_DC = 17.39 uH, R_DC = 1.05 Ohm
"""

import numpy as np
import time
from scipy.linalg import null_space

MU_0 = 4.0 * np.pi * 1e-7


# ============================================================
# Rectangular ring surface mesh (OCC)
# ============================================================

def create_rectangular_ring_mesh(a_half, b_half, trace_width, maxh,
                                  label="conductor"):
    """Create rectangular ring (frame) surface mesh using Netgen OCC.

    The ring is centered at origin in z=0 plane.
    Outer rectangle: (2*a_half + w) x (2*b_half + w)
    Inner rectangle: (2*a_half - w) x (2*b_half - w)

    Args:
        a_half: Half-width of ring centerline (x direction) [m]
        b_half: Half-height of ring centerline (y direction) [m]
        trace_width: Width of the ring trace [m]
        maxh: Maximum mesh element size [m]
        label: Boundary label for BEM

    Returns:
        NGSolve Mesh (surface triangulation in 3D)
    """
    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import Mesh

    w2 = trace_width / 2.0

    # Outer rectangle: corner at (-(a+w/2), -(b+w/2))
    a_out = a_half + w2
    b_out = b_half + w2
    wp_o = WorkPlane(Axes(Pnt(-a_out, -b_out, 0), Dir(0, 0, 1), Dir(1, 0, 0)))
    outer = wp_o.Rectangle(2 * a_out, 2 * b_out).Face()

    # Inner rectangle
    a_in = a_half - w2
    b_in = b_half - w2
    wp_i = WorkPlane(Axes(Pnt(-a_in, -b_in, 0), Dir(0, 0, 1), Dir(1, 0, 0)))
    inner = wp_i.Rectangle(2 * a_in, 2 * b_in).Face()

    ring = outer - inner
    ring.faces.name = label

    geo = OCCGeometry(ring)
    return Mesh(geo.GenerateMesh(maxh=maxh))


# ============================================================
# BEM self-inductance (LaplaceSL + Hodge decomposition)
# ============================================================

def _extract_dense(mat, n):
    """Extract dense matrix from NGSolve sparse matrix."""
    ei = mat.CreateColVector()
    col = mat.CreateColVector()
    M = np.zeros((n, n), dtype=complex)
    for i in range(n):
        ei[:] = 0; ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(n):
            M[j, i] = col[j]
    return M


def _extract_rect(mat, nrow, ncol):
    """Extract rectangular dense matrix from NGSolve matrix."""
    ei = mat.CreateRowVector()
    col = mat.CreateColVector()
    M = np.zeros((nrow, ncol))
    for i in range(ncol):
        ei[:] = 0; ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(nrow):
            M[j, i] = col[j]
    return M


def compute_ring_inductance(mesh, perimeter, trace_width,
                             sigma, thickness, label="conductor"):
    """Compute self-inductance of a rectangular ring via BEM.

    Uses LaplaceSL on HDivSurface with Hodge decomposition to extract
    the harmonic (loop) mode, then L/R ratio for normalization-independent
    inductance.

    Args:
        mesh: NGSolve surface mesh (rectangular ring)
        perimeter: Ring perimeter [m] = 4*(a + b)
        trace_width: Ring trace width [m]
        sigma: Conductivity [S/m]
        thickness: Conductor thickness [m]
        label: Boundary label

    Returns:
        L: Self-inductance [H]
        info: dict with diagnostics
    """
    from ngsolve import (HDivSurface, SurfaceL2, BilinearForm,
                         ds, TaskManager, InnerProduct)
    from ngsolve import div as ng_div
    from ngsolve.bem import LaplaceSL

    # --- Function spaces ---
    fes_J = HDivSurface(mesh, order=0)
    fes_L2 = SurfaceL2(mesh, order=0)
    n_J = fes_J.ndof
    n_v = mesh.nv
    n_f = fes_L2.ndof

    # --- D: divergence matrix (n_face x n_edge) ---
    u_J = fes_J.TrialFunction()
    q_L2 = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += ng_div(u_J.Trace()) * q_L2 * ds
    bf_D.Assemble()
    D = _extract_rect(bf_D.mat, n_f, n_J)

    # --- C: signed vertex-edge incidence ---
    C = np.zeros((n_J, n_v))
    for e_idx, edge in enumerate(mesh.edges):
        verts = list(edge.vertices)
        C[e_idx, verts[0].nr] = -1
        C[e_idx, verts[1].nr] = +1

    # --- M_J: edge mass matrix ---
    u2, v2 = fes_J.TnT()
    bf_M = BilinearForm(fes_J)
    bf_M += InnerProduct(u2.Trace(), v2.Trace()) * ds
    bf_M.Assemble()
    M_J = np.real(_extract_dense(bf_M.mat, n_J))

    # --- Hodge decomposition: harmonic mode ---
    constraint = np.vstack([D, C.T @ M_J])
    null_h = null_space(constraint, rcond=1e-10)

    if null_h.shape[1] == 0:
        raise RuntimeError("No harmonic mode found (genus must be >= 1).")

    c_h = null_h[:, 0]
    energy = c_h @ M_J @ c_h

    # --- V_A: LaplaceSL on HDivSurface ---
    j_trial = fes_J.TrialFunction()
    j_test = fes_J.TestFunction()
    with TaskManager():
        V_op = LaplaceSL(
            j_trial.Trace() * ds(label)
        ) * j_test.Trace() * ds(label)

    V_A = np.real(_extract_dense(V_op.mat, n_J))
    V_A_proj = c_h @ V_A @ c_h

    # --- Inductance via L/R ratio ---
    R_sheet = 1.0 / (sigma * thickness)
    R_loop = R_sheet * perimeter / trace_width

    LR_ratio = MU_0 * V_A_proj / (R_sheet * energy)
    L = LR_ratio * R_loop

    info = {
        'n_edges': n_J,
        'n_faces': n_f,
        'n_harmonic': null_h.shape[1],
        'div_check': np.linalg.norm(D @ c_h),
        'energy': energy,
        'V_A_proj': V_A_proj,
        'R_loop': R_loop,
    }
    return L, info


# ============================================================
# Grover mutual inductance (parallel filament pairs)
# ============================================================

def _psi(u, d):
    """Grover kernel for parallel filament mutual inductance.

    psi(u, d) = u * ln(u + sqrt(u^2 + d^2)) - sqrt(u^2 + d^2)

    For d=0: psi(u, 0) = u * ln(|u|) - |u|  (limiting case)
    """
    if d < 1e-15:
        # Degenerate case: overlapping filaments
        au = np.abs(u)
        if au < 1e-15:
            return 0.0
        return u * np.log(au) - au
    r = np.sqrt(u * u + d * d)
    return u * np.log(u + r) - r


def neumann_M_parallel_filaments(x1a, x1b, y1, x2a, x2b, y2):
    """Mutual inductance between two parallel filaments along x-axis.

    Segment 1: from (x1a, y1) to (x1b, y1) along x
    Segment 2: from (x2a, y2) to (x2b, y2) along x

    Both segments are along x-direction (same direction: x1b > x1a, x2b > x2a).
    Perpendicular distance: d = |y1 - y2|.

    Returns:
        M [H] (positive if same direction, can be negative if anti-parallel)
    """
    d = abs(y1 - y2)
    # Four-term Grover formula (Neumann integral of parallel filaments)
    # I = F(x1b-x2a) - F(x1a-x2a) - F(x1b-x2b) + F(x1a-x2b)
    M = (MU_0 / (4.0 * np.pi)) * (
        _psi(x1b - x2a, d)
        - _psi(x1a - x2a, d)
        - _psi(x1b - x2b, d)
        + _psi(x1a - x2b, d)
    )
    return M


def _mutual_side_pair(xa1, xa2, ya, xb1, xb2, yb, w, n_sub, sign):
    """Multi-filament mutual between two parallel sides.

    Distributes n_sub filaments across trace width w perpendicular to the
    side direction, computes all pairwise M, and averages.

    Args:
        xa1, xa2: start/end along side axis for side a
        ya: perpendicular coordinate (centerline) of side a
        xb1, xb2: start/end along side axis for side b
        yb: perpendicular coordinate (centerline) of side b
        w: trace width [m]
        n_sub: number of sub-filaments across width
        sign: +1 (same direction) or -1 (anti-parallel)

    Returns:
        M contribution [H] (including direction sign)
    """
    offsets = np.linspace(-w / 2 + w / (2 * n_sub),
                           w / 2 - w / (2 * n_sub), n_sub)
    M_sum = 0.0
    for dy_a in offsets:
        for dy_b in offsets:
            M_sum += neumann_M_parallel_filaments(
                xa1, xa2, ya + dy_a,
                xb1, xb2, yb + dy_b)
    return sign * M_sum / (n_sub * n_sub)


def compute_mutual_inductance_rect_loops(a_n, b_n, a_m, b_m, w=0.0,
                                          n_sub=1):
    """Mutual inductance between two coplanar rectangular loops.

    Loop n: sides at x=+-a_n, y=+-b_n (centerline)
    Loop m: sides at x=+-a_m, y=+-b_m (centerline)

    Each loop has 4 sides with circulation direction:
      Bottom: +x at y=-b   (from (-a,-b) to (+a,-b))
      Right:  +y at x=+a   (from (+a,-b) to (+a,+b))
      Top:    -x at y=+b   (from (+a,+b) to (-a,+b))
      Left:   -y at x=-a   (from (-a,+b) to (-a,-b))

    Args:
        w: trace width [m] (0 = single filament)
        n_sub: number of sub-filaments across trace width

    Returns:
        M [H]: mutual inductance (includes direction signs)
    """
    if n_sub <= 1 or w <= 0:
        # Single filament per side (original)
        M_total = 0.0
        M_total += neumann_M_parallel_filaments(-a_n, +a_n, -b_n,
                                                 -a_m, +a_m, -b_m)
        M_total -= neumann_M_parallel_filaments(-a_n, +a_n, -b_n,
                                                 -a_m, +a_m, +b_m)
        M_total -= neumann_M_parallel_filaments(-a_n, +a_n, +b_n,
                                                 -a_m, +a_m, -b_m)
        M_total += neumann_M_parallel_filaments(-a_n, +a_n, +b_n,
                                                 -a_m, +a_m, +b_m)
        M_total += neumann_M_parallel_filaments(-b_n, +b_n, +a_n,
                                                 -b_m, +b_m, +a_m)
        M_total -= neumann_M_parallel_filaments(-b_n, +b_n, +a_n,
                                                 -b_m, +b_m, -a_m)
        M_total -= neumann_M_parallel_filaments(-b_n, +b_n, -a_n,
                                                 -b_m, +b_m, +a_m)
        M_total += neumann_M_parallel_filaments(-b_n, +b_n, -a_n,
                                                 -b_m, +b_m, -a_m)
        return M_total

    M_total = 0.0

    # Horizontal pairs (along x, offsets in y)
    # Bottom-Bottom (same dir)
    M_total += _mutual_side_pair(-a_n, +a_n, -b_n, -a_m, +a_m, -b_m,
                                  w, n_sub, +1)
    # Bottom-Top (anti-par)
    M_total += _mutual_side_pair(-a_n, +a_n, -b_n, -a_m, +a_m, +b_m,
                                  w, n_sub, -1)
    # Top-Bottom (anti-par)
    M_total += _mutual_side_pair(-a_n, +a_n, +b_n, -a_m, +a_m, -b_m,
                                  w, n_sub, -1)
    # Top-Top (same dir)
    M_total += _mutual_side_pair(-a_n, +a_n, +b_n, -a_m, +a_m, +b_m,
                                  w, n_sub, +1)

    # Vertical pairs (along y, offsets in x)
    # Right-Right (same dir)
    M_total += _mutual_side_pair(-b_n, +b_n, +a_n, -b_m, +b_m, +a_m,
                                  w, n_sub, +1)
    # Right-Left (anti-par)
    M_total += _mutual_side_pair(-b_n, +b_n, +a_n, -b_m, +b_m, -a_m,
                                  w, n_sub, -1)
    # Left-Right (anti-par)
    M_total += _mutual_side_pair(-b_n, +b_n, -a_n, -b_m, +b_m, +a_m,
                                  w, n_sub, -1)
    # Left-Left (same dir)
    M_total += _mutual_side_pair(-b_n, +b_n, -a_n, -b_m, +b_m, -a_m,
                                  w, n_sub, +1)

    return M_total


# ============================================================
# Main
# ============================================================

def main():
    import sys

    print("=" * 60)
    print("14-Turn Rectangular Spiral - BEM Inductance Model")
    print("=" * 60)

    # --- Parameters (same as spiral_peec.py) ---
    N_turns = 14
    pitch = 1.253e-3
    w = 0.80e-3
    sigma = 5.8e7
    a0 = (6.475 + w * 1e3 / 2) * 1e-3   # 6.875 mm
    b0 = (27.225 + w * 1e3 / 2) * 1e-3   # 27.625 mm

    # Conductor thickness from R_DC
    from spiral_peec import compute_wire_length, generate_rectangular_spiral
    waypoints, segments = generate_rectangular_spiral(a0, b0, N_turns, pitch)
    L_wire = compute_wire_length(waypoints, segments)
    R_dc_meas = 1.05
    h = (1.0 / sigma) * L_wire / (w * R_dc_meas)

    maxh = 0.4e-3
    n_sub = 5   # sub-filaments per side for inter-turn mutual

    print(f"\nGeometry: {N_turns} turns, pitch={pitch*1e3:.3f} mm")
    print(f"  a0={a0*1e3:.3f} mm, b0={b0*1e3:.3f} mm")
    print(f"  w={w*1e3:.2f} mm, h={h*1e6:.1f} um")
    print(f"  maxh={maxh*1e3:.1f} mm, n_sub={n_sub}")

    # --- Step 1: BEM self-inductance for each turn ---
    # Use cached values from previous BEM run to avoid 2+ hour recomputation
    L_self_cached_nH = [113.58, 125.86, 138.04, 150.18, 162.32, 174.43,
                        186.63, 198.85, 211.11, 223.43, 235.77, 248.21,
                        260.68, 273.22]

    use_cache = '--no-cache' not in sys.argv
    if use_cache and len(L_self_cached_nH) == N_turns:
        print(f"\n--- BEM self-inductance per turn (cached) ---")
        L_self = np.array(L_self_cached_nH) * 1e-9
        for n in range(N_turns):
            a_n = a0 + n * pitch
            b_n = b0 + n * pitch
            print(f"  Turn {n:2d}: a={a_n*1e3:7.3f} mm, b={b_n*1e3:7.3f} mm, "
                  f"L_self={L_self[n]*1e9:8.2f} nH")
    else:
        print(f"\n--- BEM self-inductance per turn ---")
        L_self = np.zeros(N_turns)

        t0 = time.perf_counter()
        for n in range(N_turns):
            a_n = a0 + n * pitch
            b_n = b0 + n * pitch
            perimeter = 4.0 * (a_n + b_n)

            mesh = create_rectangular_ring_mesh(a_n, b_n, w, maxh)
            L_val, info = compute_ring_inductance(
                mesh, perimeter, w, sigma, h)
            L_self[n] = L_val

            print(f"  Turn {n:2d}: a={a_n*1e3:7.3f} mm, b={b_n*1e3:7.3f} mm, "
                  f"L_self={L_val*1e9:8.2f} nH  "
                  f"({info['n_edges']} edges, {info['n_harmonic']} harm)")

        dt_bem = time.perf_counter() - t0
        print(f"  BEM total: {dt_bem:.1f} s")

    # --- Step 2: Mutual inductance (Grover, multi-filament) ---
    print(f"\n--- Mutual inductance (Grover, n_sub={n_sub}) ---")
    M_matrix = np.zeros((N_turns, N_turns))
    t0 = time.perf_counter()

    for n in range(N_turns):
        a_n = a0 + n * pitch
        b_n = b0 + n * pitch
        for m in range(n + 1, N_turns):
            a_m = a0 + m * pitch
            b_m = b0 + m * pitch
            M_nm = compute_mutual_inductance_rect_loops(
                a_n, b_n, a_m, b_m, w=w, n_sub=n_sub)
            M_matrix[n, m] = M_nm
            M_matrix[m, n] = M_nm

    dt_mutual = time.perf_counter() - t0

    # Print nearest-neighbor and far couplings
    for dn in [1, 2, 5, 13]:
        if dn < N_turns:
            M_avg = np.mean([M_matrix[n, n + dn]
                             for n in range(N_turns - dn)])
            k_avg = M_avg / np.mean(L_self)
            print(f"  M(dn={dn:2d}) avg = {M_avg*1e9:8.3f} nH  "
                  f"(k = {k_avg:.3f})")
    print(f"  Grover calc: {dt_mutual:.2f} s")

    # --- Step 3: Total inductance ---
    L_self_total = np.sum(L_self)
    M_total = 2.0 * np.sum(np.triu(M_matrix, k=1))
    L_dc = L_self_total + M_total

    print(f"\n--- Total inductance ---")
    print(f"  Sum L_self   = {L_self_total*1e6:.3f} uH")
    print(f"  Sum 2*M      = {M_total*1e6:+.3f} uH")
    print(f"  L_DC (BEM)   = {L_dc*1e6:.3f} uH")
    print(f"  L_DC (PEEC)  = 13.907 uH")
    print(f"  L_DC (meas)  = 17.388 uH")
    print(f"  BEM vs meas  = {(L_dc*1e6-17.388)/17.388*100:+.1f}%")
    print(f"  PEEC vs meas = {(13.907-17.388)/17.388*100:+.1f}%")

    # --- R_DC ---
    R_sheet = 1.0 / (sigma * h)
    R_dc_total = 0.0
    for n in range(N_turns):
        a_n = a0 + n * pitch
        b_n = b0 + n * pitch
        perimeter = 4.0 * (a_n + b_n)
        R_dc_total += R_sheet * perimeter / w
    print(f"\n  R_DC = {R_dc_total:.4f} Ohm  (meas: 1.05 Ohm)")

    print("\nDone.")


if __name__ == '__main__':
    main()
