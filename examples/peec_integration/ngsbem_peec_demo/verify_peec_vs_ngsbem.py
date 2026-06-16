"""
verify_peec_vs_ngsbem.py

Cross-verification of inductance computation:
  1. PEEC (geometric Neumann formula via peec_matrices)
  2. NGSBEM (Laplace SL BEM via bem_inductance)
  3. Radia MSC (surface charge method via MMMBuilder)
  4. Analytical (for straight wire: mu_0*l/(2pi)*(ln(2l/a)-1))

Geometry: 100mm straight wire with 1mm x 1mm cross-section

This validates that the MSC kernel fix (FieldFromChargedTriangle)
produces correct interaction matrices by comparing multiple methods.

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad  # DLL loader for MKL etc.

MU_0 = 4e-7 * np.pi


def analytical_wire_inductance(length, width, height):
    """Analytical self-inductance of a straight rectangular conductor.

    Grover formula: L = (mu_0 * l) / (2*pi) * (ln(2*l / GMD) - 1 + GMD/(3*l))
    where GMD is the geometric mean distance of the cross-section.
    For square: GMD ~= 0.2235 * (w + h) (Rosa's approximation)
    """
    # GMD for rectangular cross-section (Hoer & Love, 1965)
    w, h = width, height
    if abs(w - h) < 1e-15:
        # Square cross-section
        gmd = 0.447 * w  # exp(-1/4) * sqrt(2) * w / (sqrt(pi)) ~= 0.447*w
    else:
        gmd = np.exp(-0.25) * np.sqrt(w * h)  # approximate

    L = (MU_0 * length) / (2 * np.pi) * (np.log(2 * length / gmd) - 1)
    return L


def test_peec_inductance():
    """PEEC geometric inductance (Neumann formula)."""
    from peec_matrices import PEECBuilder
    from peec_topology import PEECCircuitSolver

    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.1, 0, 0)  # 100mm
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder.add_port(n1, n2)
    topo = builder.build_topology()

    solver = PEECCircuitSolver(topo)
    Z = solver.compute_port_impedance(1e3)  # low freq for L extraction
    L = np.imag(Z) / (2 * np.pi * 1e3)
    return L


def test_radia_delta_L():
    """Radia MSC: Delta_L from magnetic core coupling.

    Uses the same geometry as demo_fasthenry_magnetic_core.py scenario 2.
    """
    import radia as rad
    from peec_msc_schur import biot_savart_finite_filament

    core_verts = np.array([
        [0.02, 0.005, -0.005], [0.08, 0.005, -0.005],
        [0.08, 0.015, -0.005], [0.02, 0.015, -0.005],
        [0.02, 0.005,  0.005], [0.08, 0.005,  0.005],
        [0.08, 0.015,  0.005], [0.02, 0.015,  0.005],
    ])

    seg_p1 = np.array([0.0, 0.0, 0.0])
    seg_p2 = np.array([0.1, 0.0, 0.0])
    seg_center = (seg_p1 + seg_p2) / 2.0
    seg_dir = np.array([1.0, 0.0, 0.0])
    seg_len = 0.1

    chi = 999.0

    rad.UtiDelAll()
    core = rad.ObjHexahedron(core_verts.tolist(), [0, 0, 0])
    rad.MatApl(core, rad.MatLin(chi))

    def h_bkg(pt):
        H = biot_savart_finite_filament(seg_p1, seg_p2, pt)
        return [MU_0 * H[0], MU_0 * H[1], MU_0 * H[2]]

    bkg = rad.ObjBckg(h_bkg)
    cnt = rad.ObjCnt([core, bkg])
    rad.Solve(cnt, 1e-6, 5000, 0)

    A = np.array(rad.Fld(core, 'a', seg_center.tolist()))
    Delta_L = np.dot(A, seg_dir) * seg_len
    M = np.array(rad.Fld(core, 'm', [0.05, 0.01, 0]))

    return Delta_L, M


def test_mmm_delta_L():
    """MMMBuilder MSC: Delta_L via standalone N matrix solve."""
    import radia as rad
    from mmm_core import MMMBuilder, MMMSolver
    from peec_msc_schur import biot_savart_finite_filament

    core_verts = np.array([
        [0.02, 0.005, -0.005], [0.08, 0.005, -0.005],
        [0.08, 0.015, -0.005], [0.02, 0.015, -0.005],
        [0.02, 0.005,  0.005], [0.08, 0.005,  0.005],
        [0.08, 0.015,  0.005], [0.02, 0.015,  0.005],
    ])

    seg_p1 = np.array([0.0, 0.0, 0.0])
    seg_p2 = np.array([0.1, 0.0, 0.0])
    seg_center = (seg_p1 + seg_p2) / 2.0
    seg_dir = np.array([1.0, 0.0, 0.0])
    seg_len = 0.1

    chi = 999.0

    builder = MMMBuilder()
    builder.add_hexahedron(core_verts)
    N_mmm, dof_off = builder.build()
    N_mmm = np.array(N_mmm)
    dof_off = np.array(dof_off)
    inv_chi = np.full(6, 1.0 / chi)

    # System matrix: (1/chi + N) sigma = H_ext
    K_msc = np.diag(inv_chi) + N_mmm

    # Face normals (MMMBuilder order)
    normals = np.array([
        [0, 0, -1], [0, 0, 1], [0, -1, 0], [0, 1, 0], [-1, 0, 0], [1, 0, 0]
    ], dtype=float)

    # Face eval points (Yano midpoint)
    v = core_verts
    face_verts = [
        [v[0], v[3], v[2], v[1]], [v[4], v[5], v[6], v[7]],
        [v[0], v[1], v[5], v[4]], [v[2], v[3], v[7], v[6]],
        [v[0], v[4], v[7], v[3]], [v[1], v[2], v[6], v[5]],
    ]
    elem_center = core_verts.mean(axis=0)
    eval_pts = np.array([(np.mean(fv, axis=0) + elem_center) / 2.0 for fv in face_verts])

    # B_pm: H normal at eval points from unit current
    B_pm = np.zeros(6)
    for i in range(6):
        H = biot_savart_finite_filament(seg_p1, seg_p2, eval_pts[i])
        B_pm[i] = np.dot(H, normals[i])

    # Solve MSC
    sigma = np.linalg.solve(K_msc, B_pm)

    # Reconstruct M and compute Delta_L via Radia A-field
    NtN = normals.T @ normals
    M_recon = np.linalg.solve(NtN, normals.T @ sigma)

    rad.UtiDelAll()
    core_m = rad.ObjHexahedron(core_verts.tolist(), M_recon.tolist())
    A = np.array(rad.Fld(core_m, 'a', seg_center.tolist()))
    Delta_L = np.dot(A, seg_dir) * seg_len

    return Delta_L, M_recon, sigma


def test_ngsbem_inductance():
    """NGSBEM: Direct BEM inductance (air, no core).

    Uses bem_inductance.py with a simple wire mesh.
    """
    try:
        import ngsolve
        from ngsolve import Mesh, BND
        from ngsolve.bem import LaplaceSL
        from ngsolve import HDivSurface, SurfaceL2, BilinearForm, LinearForm, ds, div, GridFunction, TaskManager
        from netgen.occ import Box, Pnt, Glue, OCCGeometry
    except ImportError as e:
        return None, str(e)

    # Create wire geometry: 100mm x 1mm x 1mm
    # IMPORTANT: Use Glue(faces) for surface-only mesh (no volume).
    # Volume mesh boundary extraction causes degenerate HDivSurface DOFs.
    wire = Box(Pnt(0, -0.5e-3, -0.5e-3), Pnt(0.1, 0.5e-3, 0.5e-3))

    for f in wire.faces:
        ctr = f.center
        if abs(ctr.x) < 1e-6:
            f.name = "source"
        elif abs(ctr.x - 0.1) < 1e-6:
            f.name = "sink"
        else:
            f.name = "wire"

    wire_surf = Glue(wire.faces)  # Surface-only, no volume mesh
    geo = OCCGeometry(wire_surf)
    # maxh must be <= min cross-section / 2 for well-shaped elements
    ngmesh = geo.GenerateMesh(maxh=0.5e-3)
    mesh = Mesh(ngmesh)

    nse = sum(1 for _ in mesh.Elements(BND))
    print(f"     NGSBEM mesh: {nse} surface elements, {mesh.nv} vertices")
    bnds = set(mesh.GetBoundaries())
    print(f"     Boundaries: {bnds}")

    # Direct BEM computation (same as bem_inductance but with error handling)
    fes_J = HDivSurface(mesh, order=0)
    fes_L2 = SurfaceL2(mesh, order=0)
    n_J = fes_J.ndof
    n_f = fes_L2.ndof
    print(f"     n_J={n_J}, n_f={n_f}")

    # Divergence matrix
    u_J = fes_J.TrialFunction()
    q = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()

    from scipy.sparse import coo_matrix
    rows, cols, vals = bf_D.mat.COO()
    D = coo_matrix((vals, (rows, cols)), shape=(bf_D.mat.height, bf_D.mat.width)).toarray()

    # LaplaceSL
    jt, jv = fes_J.TnT()
    with TaskManager():
        V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
    rows_sl, cols_sl, vals_sl = V_op.mat.COO()
    SL = coo_matrix((vals_sl, (rows_sl, cols_sl)),
                     shape=(V_op.mat.height, V_op.mat.width)).toarray()

    # Source/sink RHS
    f_src = LinearForm(fes_L2)
    f_src += q * ds("source")
    f_src.Assemble()
    g_src = f_src.vec.FV().NumPy().copy()
    A_src = np.sum(g_src)

    f_snk = LinearForm(fes_L2)
    f_snk += q * ds("sink")
    f_snk.Assemble()
    g_snk = f_snk.vec.FV().NumPy().copy()
    A_snk = np.sum(g_snk)

    print(f"     A_source={A_src:.6e}, A_sink={A_snk:.6e}")

    if A_src < 1e-30 or A_snk < 1e-30:
        return None, f"Source/sink faces not found (A_src={A_src}, A_snk={A_snk})"

    g = g_src / A_src - g_snk / A_snk

    # Saddle point: remove one row from D for regularity
    D_red = D[:-1, :]
    g_red = g[:-1]
    n_constraint = n_f - 1

    K = np.block([
        [SL,              D_red.T],
        [D_red, np.zeros((n_constraint, n_constraint))]
    ])

    rhs = np.zeros(n_J + n_constraint)
    rhs[n_J:] = g_red

    print(f"     K size: {K.shape}, cond: {np.linalg.cond(K):.2e}")

    from scipy.linalg import solve as scipy_solve
    x = scipy_solve(K, rhs)
    J = x[:n_J]

    L = MU_0 * J @ SL @ J
    residual = np.max(np.abs(D @ J - g))
    print(f"     residual = {residual:.2e}")

    return float(L), None


def main():
    print("=" * 70)
    print("Cross-Verification: PEEC vs NGSBEM vs Radia MSC vs MMMBuilder")
    print("=" * 70)
    print()

    # Wire parameters
    length = 0.1     # 100 mm
    width = 1e-3     # 1 mm
    height = 1e-3    # 1 mm

    # 1. Analytical
    L_analytical = analytical_wire_inductance(length, width, height)
    print(f"1. Analytical (Grover):  L = {L_analytical*1e9:.2f} nH")

    # 2. PEEC
    try:
        L_peec = test_peec_inductance()
        print(f"2. PEEC (Neumann):       L = {L_peec*1e9:.2f} nH")
    except Exception as e:
        L_peec = None
        print(f"2. PEEC: FAILED ({e})")

    # 3. NGSBEM
    try:
        L_ngsbem, err = test_ngsbem_inductance()
        if L_ngsbem is not None:
            print(f"3. NGSBEM (Laplace SL):  L = {L_ngsbem*1e9:.2f} nH")
        else:
            print(f"3. NGSBEM: SKIPPED ({err})")
    except Exception as e:
        L_ngsbem = None
        print(f"3. NGSBEM: FAILED ({e})")

    print()
    print("-" * 70)
    print("Delta_L from magnetic core (1 hex, mu_r=1000, 60x10x10mm, 10mm above)")
    print("-" * 70)

    # 4. Radia direct
    try:
        import radia as rad
        DL_radia, M_radia = test_radia_delta_L()
        print(f"4. Radia (direct Solve): Delta_L = {DL_radia*1e9:.4f} nH, "
              f"Mz = {M_radia[2]:.2f} A/m")
    except Exception as e:
        DL_radia = None
        print(f"4. Radia: FAILED ({e})")

    # 5. MMMBuilder MSC
    try:
        DL_mmm, M_mmm, sigma_mmm = test_mmm_delta_L()
        print(f"5. MMMBuilder (MSC):     Delta_L = {DL_mmm*1e9:.4f} nH, "
              f"Mz = {M_mmm[2]:.2f} A/m")
    except Exception as e:
        DL_mmm = None
        print(f"5. MMMBuilder: FAILED ({e})")

    # Summary
    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)

    if DL_radia is not None and DL_mmm is not None:
        ratio = DL_mmm / DL_radia
        print(f"Delta_L ratio (MMMBuilder/Radia) = {ratio:.6f}")
        if abs(ratio - 1.0) < 0.01:
            print("  --> PASS: MMMBuilder MSC matches Radia within 1%")
        else:
            print(f"  --> FAIL: {abs(ratio-1.0)*100:.2f}% difference")

    if L_peec is not None and L_ngsbem is not None:
        ratio_air = L_peec / L_ngsbem
        print(f"L_air ratio (PEEC/NGSBEM) = {ratio_air:.4f}")

    print()


if __name__ == "__main__":
    main()
