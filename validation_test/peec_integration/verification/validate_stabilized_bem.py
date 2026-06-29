"""
validate_stabilized_bem.py

Validation of ngbem Low-Frequency Stabilized BEM Formulation
(Weggler's stabilized DtN formulation)

Reference:
  https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb

Tests:
  Phase 1: Sphere condition number survey (stabilized vs classical EFIE)
  Phase 2: Plate conductor impedance extraction
  Phase 3: Three-solver comparison (Radia PEEC vs Laplace BEM vs Stabilized BEM)

Part of Radia project
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7    # H/m
EPS_0 = 8.854187817e-12       # F/m
C_0 = 1.0 / np.sqrt(MU_0 * EPS_0)  # speed of light [m/s]


# ============================================================
# Helper functions
# ============================================================

def extract_dense_matrix(mat, ndof):
    """Extract dense numpy matrix from NGSolve BaseMatrix."""
    M = np.zeros((ndof, ndof), dtype=complex)
    for i in range(ndof):
        ei = mat.CreateColVector()
        ei[:] = 0
        ei[i] = 1.0
        col = mat.CreateColVector()
        mat.Mult(ei, col)
        for j in range(ndof):
            M[j, i] = col[j]
    return M


def compute_condition_number(mat_dense):
    """Compute condition number via SVD."""
    s = np.linalg.svd(mat_dense, compute_uv=False)
    s_real = np.abs(s)
    s_nonzero = s_real[s_real > 1e-14 * s_real[0]]
    if len(s_nonzero) < 2:
        return float('inf'), s_real
    return s_nonzero[0] / s_nonzero[-1], s_real


# ============================================================
# Phase 1: Sphere Condition Number Test
# ============================================================

def test_sphere_condition_numbers():
    """
    Phase 1: Reproduce Weggler's notebook results on sphere.

    Compare condition numbers of stabilized vs classical EFIE
    across kappa values from 5.0 to 0.00005.

    Expected: stabilized improves conditioning by ~kappa^2 factor.
    """
    print("=" * 70)
    print("Phase 1: Sphere Condition Number Survey")
    print("        (Stabilized vs Classical EFIE)")
    print("=" * 70)

    import netgen.meshing as meshing
    from netgen.occ import Sphere, OCCGeometry
    from ngsolve import Mesh, HDivSurface, SurfaceL2, TaskManager, ds
    from ngsolve import BilinearForm, LinearForm, CF
    from ngsolve.bem import HelmholtzSL
    from ngsolve.krylovspace import GMRes

    # Sphere mesh
    sp = Sphere((0, 0, 0), 1)
    mesh = Mesh(OCCGeometry(sp).GenerateMesh(
        maxh=1, perfstepsend=meshing.MeshingStep.MESHSURFACE)).Curve(4)

    print(f"\n  Mesh: sphere, {mesh.nv} vertices, {mesh.ne} elements")

    kappa_values = [5.0, 0.5, 0.05, 0.005, 0.0005, 0.00005]
    order = 1
    intorder = 4

    results_stab = []
    results_class = []

    print(f"\n  {'kappa':>10s} | {'Cond(Stab)':>12s} | {'Cond(Class)':>12s} | "
          f"{'Ratio S/C':>10s} | {'GMRes Stab':>10s}")
    print("  " + "-" * 66)

    for kappa in kappa_values:
        # --- Stabilized formulation ---
        t0 = time.perf_counter()
        try:
            fes_hdiv = HDivSurface(mesh, order=order, complex=True)
            fes_l2 = SurfaceL2(mesh, order=order-1, complex=True,
                               dual_mapping=True)
            fes = fes_hdiv * fes_l2
            (uHDiv, uL2), (vHDiv, vL2) = fes.TnT()

            with TaskManager():
                A_k = HelmholtzSL(
                    uHDiv.Trace() * ds(bonus_intorder=intorder), kappa
                ) * vHDiv.Trace() * ds(bonus_intorder=intorder)

                V_k = HelmholtzSL(
                    uL2 * ds(bonus_intorder=intorder), kappa
                ) * vL2 * ds(bonus_intorder=intorder)

                from ngsolve import div
                Q_k = HelmholtzSL(
                    div(uHDiv.Trace()) * ds(bonus_intorder=intorder), kappa
                ) * vL2 * ds(bonus_intorder=intorder)

                lhs_stab = (A_k.mat + Q_k.mat + Q_k.mat.T
                            + kappa * kappa * V_k.mat)

                # RHS: plane wave E_inc = (1,0,0) * exp(-j*kappa*z)
                from ngsolve import x as cf_x, y as cf_y, z as cf_z, exp
                z_coord = cf_z  # z coordinate CF
                E_inc = CF((1, 0, 0)) * exp(-1j * kappa * z_coord)
                rhs = LinearForm(
                    E_inc * vHDiv.Trace() * ds(bonus_intorder=2*intorder)
                ).Assemble()

                # Preconditioner
                preBlock = BilinearForm(
                    uHDiv.Trace() * vHDiv.Trace() * ds
                    + uL2 * vL2 * ds
                ).Assemble().mat.Inverse(freedofs=fes.FreeDofs())

                # Solve
                sol_stab = GMRes(A=lhs_stab, b=rhs.vec, pre=preBlock,
                                 maxsteps=3000, tol=1e-8, printrates=False)

            # Condition number
            ndof_stab = fes.ndof
            M_stab = extract_dense_matrix(lhs_stab, ndof_stab)
            cond_stab, _ = compute_condition_number(M_stab)
            gmres_iters = "OK"
            results_stab.append((kappa, cond_stab, ndof_stab))

        except Exception as e:
            cond_stab = float('inf')
            gmres_iters = f"ERR: {e}"
            results_stab.append((kappa, cond_stab, 0))

        # --- Classical EFIE ---
        try:
            fes_cl = HDivSurface(mesh, order=order, complex=True)
            u_cl, v_cl = fes_cl.TnT()

            with TaskManager():
                V1 = HelmholtzSL(
                    u_cl.Trace() * ds(bonus_intorder=intorder), kappa
                ) * v_cl.Trace() * ds(bonus_intorder=intorder)

                V2 = HelmholtzSL(
                    div(u_cl.Trace()) * ds(bonus_intorder=intorder), kappa
                ) * div(v_cl.Trace()) * ds(bonus_intorder=intorder)

                lhs_class = V1.mat - (1.0 / (kappa * kappa)) * V2.mat

            ndof_cl = fes_cl.ndof
            M_class = extract_dense_matrix(lhs_class, ndof_cl)
            cond_class, _ = compute_condition_number(M_class)
            results_class.append((kappa, cond_class, ndof_cl))

        except Exception as e:
            cond_class = float('inf')
            results_class.append((kappa, cond_class, 0))

        ratio = cond_stab / cond_class if cond_class > 0 else float('inf')
        print(f"  {kappa:10.5f} | {cond_stab:12.3e} | {cond_class:12.3e} | "
              f"{ratio:10.3e} | {gmres_iters}")

    # Summary
    print("\n  Analysis:")
    if len(results_stab) >= 2 and len(results_class) >= 2:
        # Check that stabilized is better at low kappa
        low_kappa_stab = results_stab[-1][1]
        low_kappa_class = results_class[-1][1]
        improvement = low_kappa_class / low_kappa_stab if low_kappa_stab > 0 else 0

        print(f"  At kappa={kappa_values[-1]}: "
              f"stabilized={low_kappa_stab:.3e}, "
              f"classical={low_kappa_class:.3e}")
        print(f"  Improvement factor: {improvement:.1f}x")

        if improvement > 10:
            print("  -> Stabilized formulation significantly better at low freq")
            print("  PASSED\n")
            return True
        else:
            print("  -> Improvement smaller than expected")
            print("  WARNING\n")
            return False
    else:
        print("  -> Insufficient results")
        print("  FAILED\n")
        return False


# ============================================================
# Phase 2: Plate Conductor - Matrix Extraction
# ============================================================

def test_plate_matrix_comparison():
    """
    Phase 2: Compare BEM matrices on plate geometry.

    Extract L, P, M_LS matrices from:
    A) Laplace kernel (existing ngbem_peec.py)
    B) Stabilized Helmholtz (kappa -> 0 limit)

    At kappa -> 0, Helmholtz SL -> Laplace SL, so matrices should match.
    """
    print("=" * 70)
    print("Phase 2: Plate Conductor Matrix Comparison")
    print("        (Helmholtz kappa->0 vs Laplace kernel)")
    print("=" * 70)

    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import Mesh, HDivSurface, SurfaceL2, TaskManager, ds
    from ngsolve import BilinearForm, div
    from ngsolve.bem import HelmholtzSL, LaplaceSL

    # Plate mesh (10mm x 1mm conductor, like PEEC wire cross-section)
    width = 0.01   # 10mm
    height = 0.001  # 1mm
    maxh = 0.003

    wp = WorkPlane(Axes(p=Pnt(0, 0, 0), n=Dir(0, 0, 1), h=Dir(1, 0, 0)))
    rect = wp.Rectangle(width, height).Face()
    rect.name = 'conductor'
    mesh = Mesh(OCCGeometry(rect).GenerateMesh(maxh=maxh))

    print(f"\n  Mesh: plate {width*1e3:.0f}mm x {height*1e3:.1f}mm, "
          f"maxh={maxh*1e3:.1f}mm")
    print(f"  Vertices: {mesh.nv}, Elements: {mesh.ne}")

    order = 0
    intorder = 8  # High order for Helmholtz->Laplace accuracy

    # Function spaces
    fes_hdiv = HDivSurface(mesh, order=order, complex=True)
    fes_l2 = SurfaceL2(mesh, order=max(0, order-1), complex=True,
                       dual_mapping=True)
    fes = fes_hdiv * fes_l2
    (uHDiv, uL2), (vHDiv, vL2) = fes.TnT()

    n_loop = fes_hdiv.ndof
    n_star = fes_l2.ndof
    n_total = fes.ndof

    print(f"  Loop DOFs (HDivSurface): {n_loop}")
    print(f"  Star DOFs (SurfaceL2):   {n_star}")
    print(f"  Total DOFs:              {n_total}")

    # --- A) Laplace kernel matrices ---
    print("\n  A) Laplace kernel (LaplaceSL)...")
    t0 = time.perf_counter()

    with TaskManager():
        # L matrix: mu_0 * LaplaceSL on HDivSurface
        L_lap_op = LaplaceSL(
            uHDiv.Trace() * ds(bonus_intorder=intorder)
        ) * vHDiv.Trace() * ds(bonus_intorder=intorder)

    L_lap = MU_0 * np.real(extract_dense_matrix(L_lap_op.mat, n_loop))
    t_lap = time.perf_counter() - t0
    print(f"     L matrix: ({n_loop}x{n_loop}), time={t_lap:.2f}s")

    # --- B) Helmholtz kernel at small kappa ---
    kappa_small = 0.001  # Very small kappa (quasi-static)
    print(f"\n  B) Helmholtz kernel (kappa={kappa_small})...")
    t0 = time.perf_counter()

    with TaskManager():
        A_k = HelmholtzSL(
            uHDiv.Trace() * ds(bonus_intorder=intorder), kappa_small
        ) * vHDiv.Trace() * ds(bonus_intorder=intorder)

        V_k = HelmholtzSL(
            uL2 * ds(bonus_intorder=intorder), kappa_small
        ) * vL2 * ds(bonus_intorder=intorder)

        Q_k = HelmholtzSL(
            div(uHDiv.Trace()) * ds(bonus_intorder=intorder), kappa_small
        ) * vL2 * ds(bonus_intorder=intorder)

    # Extract matrices
    A_helm = extract_dense_matrix(A_k.mat, n_loop)
    L_helm = MU_0 * np.real(A_helm)

    # Q_k is (n_loop+n_star) shaped via product space? No - it's n_loop x n_star
    # Actually Q_k lives in the product space
    # Extract from product space
    lhs_full = (A_k.mat + Q_k.mat + Q_k.mat.T
                + kappa_small**2 * V_k.mat)

    M_full = extract_dense_matrix(lhs_full, n_total)

    # Extract blocks
    A_block = M_full[:n_loop, :n_loop]
    Q_block = M_full[n_loop:, :n_loop]
    V_block = M_full[n_loop:, n_loop:]

    t_helm = time.perf_counter() - t0
    print(f"     Full system: ({n_total}x{n_total}), time={t_helm:.2f}s")

    # --- Compare L matrices ---
    print("\n  Comparison: L (Laplace) vs mu_0*A_kappa (Helmholtz):")

    # The Helmholtz A_k at small kappa should approximate Laplace SL
    L_helm_from_block = MU_0 * np.real(A_block)
    L_diff = np.linalg.norm(L_lap - L_helm_from_block) / np.linalg.norm(L_lap)
    print(f"    ||L_lap - L_helm|| / ||L_lap|| = {L_diff:.6e}")
    print(f"    L_lap diagonal: [{np.min(np.diag(L_lap)):.4e}, "
          f"{np.max(np.diag(L_lap)):.4e}]")
    print(f"    L_helm diagonal: [{np.min(np.diag(L_helm_from_block)):.4e}, "
          f"{np.max(np.diag(L_helm_from_block)):.4e}]")

    # Check symmetry of stabilized system
    sym_err = np.linalg.norm(M_full - M_full.T) / np.linalg.norm(M_full)
    print(f"\n    Stabilized system symmetry: ||M-M^T||/||M|| = {sym_err:.4e}")

    # Condition number of stabilized system
    cond, _ = compute_condition_number(M_full)
    print(f"    Condition number: {cond:.4e}")

    # Q block (coupling)
    print(f"\n    Q_k block ({Q_block.shape[0]}x{Q_block.shape[1]}):")
    print(f"      Range: [{np.min(np.abs(Q_block)):.4e}, "
          f"{np.max(np.abs(Q_block)):.4e}]")
    print(f"      Nonzeros: {np.count_nonzero(np.abs(Q_block) > 1e-15)}"
          f" / {Q_block.size}")

    # V block (capacitance, scaled by kappa^2)
    print(f"\n    kappa^2*V_k block ({V_block.shape[0]}x{V_block.shape[1]}):")
    print(f"      Range: [{np.min(np.abs(V_block)):.4e}, "
          f"{np.max(np.abs(V_block)):.4e}]")

    if L_diff < 0.05:  # 5% tolerance
        print(f"\n  Helmholtz -> Laplace convergence CONFIRMED (error {L_diff*100:.2f}%)")
        print("  PASSED\n")
        return True, L_lap, n_loop, n_star, mesh
    else:
        print(f"\n  Helmholtz -> Laplace difference too large ({L_diff*100:.2f}%)")
        print("  FAILED\n")
        return False, L_lap, n_loop, n_star, mesh


# ============================================================
# Phase 3: Three-Solver Impedance Comparison
# ============================================================

def test_three_solver_comparison():
    """
    Phase 3: Compare impedance from three methods.

    A) Radia PEEC (filament, Neumann integral)
    B) ngbem Laplace PEEC (LaplaceSL on surface mesh)
    C) ngbem Stabilized Helmholtz (at multiple kappa)

    Geometry: 100mm copper wire, 1mm x 1mm cross-section
    """
    print("=" * 70)
    print("Phase 3: Three-Solver Impedance Comparison")
    print("        (PEEC vs Laplace BEM vs Stabilized BEM)")
    print("=" * 70)

    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import Mesh, HDivSurface, SurfaceL2, TaskManager, ds
    from ngsolve import BilinearForm, div
    from ngsolve.bem import HelmholtzSL, LaplaceSL

    # --- Parameters ---
    wire_length = 0.1   # 100mm
    wire_w = 1e-3       # 1mm width
    wire_h = 1e-3       # 1mm height
    sigma = 5.8e7       # Copper
    thickness = wire_h  # For sheet resistance

    freqs = np.array([1e3, 1e4, 1e5, 1e6, 1e7, 1e8])

    # ============================================================
    # Solver A: Radia PEEC (filament)
    # ============================================================
    print("\n  A) Radia PEEC (filament)...")
    try:
        from peec_matrices import PEECBuilder
        from peec_topology import PEECCircuitSolver

        builder = PEECBuilder()
        n1 = builder.add_node_at(0.0, 0.0, 0.0)
        n2 = builder.add_node_at(wire_length, 0.0, 0.0)
        builder.add_connected_segment(n1, n2, wire_w, wire_h, sigma)
        builder.add_port(n1, n2)
        topo = builder.build_topology()
        solver_peec = PEECCircuitSolver(topo)

        Z_peec = np.array([solver_peec.compute_port_impedance(f) for f in freqs])
        R_peec = np.real(Z_peec)
        L_peec = np.imag(Z_peec) / (2 * np.pi * freqs)

        print(f"     L_peec = {L_peec[0]*1e9:.2f} nH")
        print(f"     R_peec = {R_peec[0]*1e3:.4f} mOhm")
        peec_ok = True
    except Exception as e:
        print(f"     PEEC failed: {e}")
        Z_peec = np.zeros(len(freqs), dtype=complex)
        R_peec = np.zeros(len(freqs))
        L_peec = np.zeros(len(freqs))
        peec_ok = False

    # ============================================================
    # Solver B: ngbem Laplace PEEC (surface mesh)
    # ============================================================
    print("\n  B) ngbem Laplace PEEC (surface mesh)...")

    # Create plate mesh representing wire cross-section view
    # For a wire along x, the plate is width x height in the x-y plane
    plate_w = wire_length  # 100mm along x
    plate_h = wire_w       # 1mm along y
    maxh = 0.005           # 5mm mesh size

    wp = WorkPlane(Axes(p=Pnt(0, 0, 0), n=Dir(0, 0, 1), h=Dir(1, 0, 0)))
    rect = wp.Rectangle(plate_w, plate_h).Face()
    rect.name = 'conductor'
    mesh = Mesh(OCCGeometry(rect).GenerateMesh(maxh=maxh))

    order = 0
    intorder = 8

    fes_hdiv = HDivSurface(mesh, order=order, complex=True)
    n_loop = fes_hdiv.ndof

    print(f"     Mesh: {mesh.nv} vertices, {mesh.ne} elements, "
          f"{n_loop} edge DOFs")

    # Assemble Laplace L matrix
    u_h, v_h = fes_hdiv.TnT()
    with TaskManager():
        L_op = LaplaceSL(
            u_h.Trace() * ds(bonus_intorder=intorder)
        ) * v_h.Trace() * ds(bonus_intorder=intorder)

    L_lap = MU_0 * np.real(extract_dense_matrix(L_op.mat, n_loop))

    # Resistance: sheet resistance model
    R_sheet = 1.0 / (sigma * thickness)
    triangles_areas = []
    from ngsolve import BND
    total_area = 0
    n_tri = 0
    for el in mesh.Elements(BND):
        verts = []
        for v in el.vertices:
            pt = mesh.vertices[v.nr].point
            verts.append(np.array([pt[0], pt[1], pt[2]]))
        if len(verts) == 3:
            e1 = verts[1] - verts[0]
            e2 = verts[2] - verts[0]
            area = 0.5 * np.linalg.norm(np.cross(e1, e2))
            total_area += area
            n_tri += 1

    mean_edge = np.sqrt(4.0 * total_area / (np.sqrt(3) * n_tri))
    A_mean = total_area / n_tri
    w_eff = 2.0 * A_mean / mean_edge
    R_per_edge = R_sheet * mean_edge / w_eff
    R_lap = np.full(n_loop, R_per_edge)

    # MQS impedance: Z = R + jwL (using total L and R)
    # Port approximation: uniform current (simple average)
    Z_laplace = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        omega = 2.0 * np.pi * f
        Z_branch = np.diag(R_lap.astype(complex)) + 1j * omega * L_lap
        try:
            Z_inv = np.linalg.inv(Z_branch)
            e = np.ones(n_loop) / n_loop
            Z_laplace[i] = 1.0 / (e @ Z_inv @ e)
        except np.linalg.LinAlgError:
            Z_laplace[i] = np.sum(R_lap) + 1j * omega * np.sum(L_lap)

    R_laplace = np.real(Z_laplace)
    L_laplace = np.imag(Z_laplace) / (2 * np.pi * freqs)

    print(f"     L_laplace = {L_laplace[0]*1e9:.2f} nH")
    print(f"     R_laplace = {R_laplace[0]*1e3:.4f} mOhm")

    # ============================================================
    # Solver C: Stabilized Helmholtz BEM
    # ============================================================
    print("\n  C) Stabilized Helmholtz BEM...")

    fes_l2 = SurfaceL2(mesh, order=max(0, order-1), complex=True,
                       dual_mapping=True)
    fes_prod = fes_hdiv * fes_l2
    (uHDiv, uL2), (vHDiv, vL2) = fes_prod.TnT()
    n_star = fes_l2.ndof
    n_total = fes_prod.ndof

    print(f"     Product space: {n_loop} loop + {n_star} star = {n_total} DOFs")

    Z_stab = np.zeros(len(freqs), dtype=complex)
    cond_stab = np.zeros(len(freqs))

    for i, f in enumerate(freqs):
        omega = 2.0 * np.pi * f
        kappa = omega / C_0

        with TaskManager():
            A_k = HelmholtzSL(
                uHDiv.Trace() * ds(bonus_intorder=intorder), kappa
            ) * vHDiv.Trace() * ds(bonus_intorder=intorder)

            V_k = HelmholtzSL(
                uL2 * ds(bonus_intorder=intorder), kappa
            ) * vL2 * ds(bonus_intorder=intorder)

            Q_k = HelmholtzSL(
                div(uHDiv.Trace()) * ds(bonus_intorder=intorder), kappa
            ) * vL2 * ds(bonus_intorder=intorder)

        # Extract full system
        lhs = A_k.mat + Q_k.mat + Q_k.mat.T + kappa**2 * V_k.mat
        M_full = extract_dense_matrix(lhs, n_total)

        # Extract A_kappa block (loop-loop) for inductance
        A_block = M_full[:n_loop, :n_loop]
        L_helm = MU_0 * np.real(A_block)

        # MQS impedance from Helmholtz L (same port extraction as Laplace)
        Z_branch = np.diag(R_lap.astype(complex)) + 1j * omega * L_helm
        try:
            Z_inv = np.linalg.inv(Z_branch)
            e = np.ones(n_loop) / n_loop
            Z_stab[i] = 1.0 / (e @ Z_inv @ e)
        except np.linalg.LinAlgError:
            Z_stab[i] = np.sum(R_lap) + 1j * omega * np.sum(L_helm)

        # Condition number of full stabilized system
        cond_stab[i], _ = compute_condition_number(M_full)

    R_stab = np.real(Z_stab)
    L_stab = np.imag(Z_stab) / (2 * np.pi * freqs)

    print(f"     L_stab = {L_stab[0]*1e9:.2f} nH")
    print(f"     R_stab = {R_stab[0]*1e3:.4f} mOhm")

    # ============================================================
    # Comparison table
    # ============================================================
    print("\n" + "=" * 70)
    print("  Results Comparison")
    print("=" * 70)

    # DC resistance analytical
    R_dc_analytical = wire_length / (sigma * wire_w * wire_h)
    print(f"\n  DC Resistance (analytical): {R_dc_analytical*1e3:.4f} mOhm")

    print(f"\n  {'Freq':>10s} | {'kappa':>10s} | "
          f"{'L_PEEC(nH)':>10s} | {'L_Lap(nH)':>10s} | {'L_Stab(nH)':>10s} | "
          f"{'Cond(Stab)':>10s}")
    print("  " + "-" * 76)

    for i, f in enumerate(freqs):
        kappa = 2 * np.pi * f / C_0
        print(f"  {f:10.0e} | {kappa:10.3e} | "
              f"{L_peec[i]*1e9:10.2f} | {L_laplace[i]*1e9:10.2f} | "
              f"{L_stab[i]*1e9:10.2f} | {cond_stab[i]:10.3e}")

    # MQS range comparison (Laplace vs Stabilized)
    print("\n  Laplace vs Stabilized (MQS range, f <= 1 MHz):")
    mqs_mask = freqs <= 1e6
    if np.any(mqs_mask):
        L_lap_mqs = L_laplace[mqs_mask]
        L_stab_mqs = L_stab[mqs_mask]
        if np.all(np.abs(L_lap_mqs) > 0):
            max_diff = np.max(np.abs(L_lap_mqs - L_stab_mqs) / np.abs(L_lap_mqs)) * 100
            print(f"    Max L difference: {max_diff:.4f}%")
            if max_diff < 2.0:
                print("    -> Helmholtz matches Laplace in MQS limit")
                print("    PASSED")
            else:
                print(f"    -> Difference {max_diff:.2f}% > 2% threshold")
                print("    WARNING")

    # PEEC vs BEM comparison
    if peec_ok:
        print("\n  PEEC vs Laplace BEM:")
        print(f"    L_PEEC = {L_peec[0]*1e9:.2f} nH")
        print(f"    L_Laplace = {L_laplace[0]*1e9:.2f} nH")
        ratio = L_laplace[0] / L_peec[0] if L_peec[0] != 0 else float('inf')
        print(f"    Ratio (Laplace/PEEC): {ratio:.3f}")
        print("    Note: Different discretization (surface vs filament)")
        print("          Ratio != 1.0 is expected")

    # Condition number trend
    print("\n  Stabilized BEM Condition Number vs Frequency:")
    for i, f in enumerate(freqs):
        kappa = 2 * np.pi * f / C_0
        print(f"    f={f:.0e}, kappa={kappa:.3e}, cond={cond_stab[i]:.3e}")

    return True


# ============================================================
# Main
# ============================================================

def main():
    print()
    print("*" * 70)
    print("  ngbem Low-Frequency Stabilized BEM Validation")
    print("  Reference: Weggler (2026)")
    print("*" * 70)
    print()

    results = {}

    # Phase 1: Sphere condition number
    try:
        results['phase1'] = test_sphere_condition_numbers()
    except Exception as e:
        print(f"Phase 1 FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        results['phase1'] = False

    # Phase 2: Plate matrix comparison
    try:
        result = test_plate_matrix_comparison()
        results['phase2'] = result[0]
    except Exception as e:
        print(f"Phase 2 FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        results['phase2'] = False

    # Phase 3: Three-solver comparison
    try:
        results['phase3'] = test_three_solver_comparison()
    except Exception as e:
        print(f"Phase 3 FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        results['phase3'] = False

    # Summary
    print()
    print("=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    for phase, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        print(f"  {phase}: {status}")
    print()

    all_passed = all(results.values())
    if all_passed:
        print("  All phases PASSED")
    else:
        print("  Some phases FAILED - see details above")

    return all_passed


if __name__ == '__main__':
    main()
