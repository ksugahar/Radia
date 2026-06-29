"""
verify_loop_peec_vs_ngsbem.py

Cross-verification of loop inductance: PEEC vs NGSBEM vs Analytical.

Geometry: 10mm x 10mm square loop, trace width 1mm, flat (2D).
Same geometry as compute_L_final.py.

Reference: FastHenry L ~ 21.6 nH, Grover L ~ 24 nH

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad  # DLL loader

MU_0 = 4e-7 * np.pi

# Geometry
SIDE = 0.01       # 10 mm loop side length (centerline)
TRACE_W = 1e-3    # 1 mm trace width
THICKNESS = 35e-6  # 35 um (for PEEC resistance, not inductance)
SIGMA = 5.8e7     # copper


def analytical_square_loop(side, trace_w):
    """Grover formula for square loop self-inductance.

    L = (2*mu_0*a/pi) * [ln(2*a/GMD) - 0.774]
    where a = side length, GMD = geometric mean distance of cross-section.
    For thin strip: GMD ~= 0.2235 * trace_w
    """
    a = side
    gmd = 0.2235 * trace_w  # thin strip approximation
    L = (2 * MU_0 * a / np.pi) * (np.log(2 * a / gmd) - 0.774)
    return L


def test_peec_loop():
    """PEEC loop inductance via PEECBuilder."""
    from peec_matrices import PEECBuilder
    from peec_topology import PEECCircuitSolver

    builder = PEECBuilder()

    # Square loop: 4 corners at centerline of traces
    # Use 5 nodes: split one side into two segments with a port gap
    half = SIDE / 2
    n1 = builder.add_node_at(-half, -half, 0)     # port node A
    n1b = builder.add_node_at(-half, -half, 0)    # port node B (same location, different node)
    n2 = builder.add_node_at(+half, -half, 0)
    n3 = builder.add_node_at(+half, +half, 0)
    n4 = builder.add_node_at(-half, +half, 0)

    # 4 sides: n1b -> n2 -> n3 -> n4 -> n1
    builder.add_connected_segment(n1b, n2, TRACE_W, THICKNESS, sigma=SIGMA)
    builder.add_connected_segment(n2, n3, TRACE_W, THICKNESS, sigma=SIGMA)
    builder.add_connected_segment(n3, n4, TRACE_W, THICKNESS, sigma=SIGMA)
    builder.add_connected_segment(n4, n1, TRACE_W, THICKNESS, sigma=SIGMA)

    builder.add_port(n1, n1b)  # Port across the gap

    topo = builder.build_topology()
    solver = PEECCircuitSolver(topo)
    Z = solver.compute_port_impedance(1e3)
    L = np.imag(Z) / (2 * np.pi * 1e3)
    R = np.real(Z)
    return L, R


def test_ngsbem_loop():
    """NGSBEM loop inductance via Hodge decomposition (same as compute_L_final.py)."""
    try:
        from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
        from ngsolve import (Mesh, HDivSurface, SurfaceL2, BilinearForm,
                             ds, BND, TaskManager, InnerProduct)
        from ngsolve import div as ng_div
        from ngsolve.bem import HelmholtzSL
        from scipy.linalg import null_space
    except ImportError as e:
        return None, str(e)

    hw = TRACE_W / 2.0

    # Frame geometry (same as compute_L_final.py)
    wp_o = WorkPlane(Axes(Pnt(-hw, -hw, 0), Dir(0, 0, 1), Dir(1, 0, 0)))
    outer = wp_o.Rectangle(SIDE + TRACE_W, SIDE + TRACE_W).Face()
    wp_i = WorkPlane(Axes(Pnt(hw, hw, 0), Dir(0, 0, 1), Dir(1, 0, 0)))
    inner = wp_i.Rectangle(SIDE - TRACE_W, SIDE - TRACE_W).Face()
    frame = outer - inner
    frame.faces.name = "conductor"
    geo = OCCGeometry(frame)

    mesh = Mesh(geo.GenerateMesh(maxh=0.001))
    n_el = sum(1 for _ in mesh.Elements(BND))

    fes_J = HDivSurface(mesh, order=0)
    fes_L2 = SurfaceL2(mesh, order=0)
    n_J = fes_J.ndof
    n_f = fes_L2.ndof
    n_v = mesh.nv
    print(f"     NGSBEM: {n_el} el, {n_J} edges, {n_f} faces")

    # Helper: extract dense matrix
    def extract_dense_mat(mat, n):
        ei = mat.CreateColVector()
        col = mat.CreateColVector()
        M = np.zeros((n, n), dtype=complex)
        for i in range(n):
            ei[:] = 0; ei[i] = 1.0
            mat.Mult(ei, col)
            for j in range(n):
                M[j, i] = col[j]
        return M

    def extract_rect_mat(mat, nrow, ncol):
        ei = mat.CreateRowVector()
        col = mat.CreateColVector()
        M = np.zeros((nrow, ncol))
        for i in range(ncol):
            ei[:] = 0; ei[i] = 1.0
            mat.Mult(ei, col)
            for j in range(nrow):
                M[j, i] = col[j]
        return M

    # D (divergence) and C (incidence) matrices
    u_J = fes_J.TrialFunction()
    q = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += ng_div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D = extract_rect_mat(bf_D.mat, n_f, n_J)

    C = np.zeros((n_J, n_v))
    for e_idx, edge in enumerate(mesh.edges):
        verts = list(edge.vertices)
        C[e_idx, verts[0].nr] = -1
        C[e_idx, verts[1].nr] = +1

    # M_J (mass matrix)
    u2, v2 = fes_J.TnT()
    bf_M = BilinearForm(fes_J)
    bf_M += InnerProduct(u2.Trace(), v2.Trace()) * ds
    bf_M.Assemble()
    M_J = np.real(extract_dense_mat(bf_M.mat, n_J))

    # Harmonic mode (Hodge decomposition)
    constraint = np.vstack([D, C.T @ M_J])
    null_h = null_space(constraint, rcond=1e-10)
    if null_h.shape[1] == 0:
        return None, "No harmonic mode found"
    c_h = null_h[:, 0]
    energy = c_h @ M_J @ c_h

    # V_A via HelmholtzSL (Laplace limit)
    kappa = 0.01
    u3, v3 = fes_J.TnT()
    with TaskManager():
        V1_op = HelmholtzSL(
            u3.Trace() * ds("conductor", bonus_intorder=6), kappa
        ) * v3.Trace() * ds("conductor", bonus_intorder=6)
    V1_mat = extract_dense_mat(V1_op.mat, n_J)
    V_A = np.real(c_h @ V1_mat @ c_h)

    # L via L/R ratio (normalization-independent, from compute_L_final.py)
    R_sheet = 1.0 / (SIGMA * THICKNESS)
    perimeter = 4 * SIDE
    R_loop = R_sheet * perimeter / TRACE_W
    LR_ratio = MU_0 * V_A / (R_sheet * energy)
    L = LR_ratio * R_loop

    return float(L), None


def main():
    print("=" * 70)
    print("Loop Inductance: PEEC vs NGSBEM vs Analytical")
    print(f"Geometry: {SIDE*1e3:.0f}mm x {SIDE*1e3:.0f}mm square loop, "
          f"trace={TRACE_W*1e3:.1f}mm, t={THICKNESS*1e6:.0f}um")
    print("=" * 70)
    print()

    # Analytical
    L_ana = analytical_square_loop(SIDE, TRACE_W)
    print(f"1. Analytical (Grover):  L = {L_ana*1e9:.2f} nH")

    # PEEC
    try:
        L_peec, R_peec = test_peec_loop()
        print(f"2. PEEC (Neumann):       L = {L_peec*1e9:.2f} nH, R = {R_peec*1e3:.4f} mOhm")
    except Exception as e:
        L_peec = None
        print(f"2. PEEC: FAILED ({e})")

    # NGSBEM
    try:
        L_ngsbem, err = test_ngsbem_loop()
        if L_ngsbem is not None:
            print(f"3. NGSBEM (Laplace SL):  L = {L_ngsbem*1e9:.2f} nH")
        else:
            print(f"3. NGSBEM: SKIPPED ({err})")
    except Exception as e:
        L_ngsbem = None
        print(f"3. NGSBEM: FAILED ({e})")

    # Reference
    print(f"\n   FastHenry reference:  L ~ 21.6 nH")

    # Comparison
    print()
    if L_peec and L_ana:
        print(f"   PEEC/Analytical = {L_peec/L_ana:.4f}")
    if L_ngsbem and L_ana:
        print(f"   NGSBEM/Analytical = {L_ngsbem/L_ana:.4f}")
    if L_peec and L_ngsbem:
        print(f"   PEEC/NGSBEM = {L_peec/L_ngsbem:.4f}")


if __name__ == "__main__":
    main()
