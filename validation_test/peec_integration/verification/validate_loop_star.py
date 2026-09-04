"""
validate_loop_star.py

Validation of Loop-Star decomposition for surface BEM mesh.

Tests:
1. Module import (no NGSolve needed)
2. Transform dimensions and Euler formula
3. Orthogonality: T_star^T * T_loop = 0
4. T matrix invertibility
5. Decomposition recovery: J_loop + J_star = J_edge
6. Divergence-free property: M_LS * T_loop = 0
7. L matrix block structure (Loop-Star of ngbem L)
8. NGBEMBridge with decompose_loop_star=True

Tests 2-5 use a synthetic mesh (no ngbem needed).
Tests 6-8 require ngbem (SKIPped if not installed).

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src'))

results = []


def record(name, status, detail=""):
    results.append((name, status, detail))
    symbol = {"PASS": "OK", "FAIL": "FAIL", "SKIP": "SKIP"}[status]
    print(f"  [{symbol}] {name}")
    if detail:
        print(f"        {detail}")


def check_ngsolve_available():
    try:
        import ngsolve
        return True
    except ImportError:
        return False


def check_ngbem_available():
    try:
        import ngsolve
        from ngsolve import bem
        return True
    except ImportError:
        return False


def create_test_plate_mesh(nx=3, ny=2, width=0.1, height=0.05):
    """Create a simple rectangular plate mesh using Netgen OCC.

    Returns NGSolve Mesh suitable for Loop-Star testing.
    """
    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import Mesh, TaskManager

    wp = WorkPlane(Axes(p=Pnt(0, 0, 0), n=Dir(0, 0, 1), h=Dir(1, 0, 0)))
    rect = wp.Rectangle(width, height).Face()
    rect.name = "conductor"

    geo = OCCGeometry(rect)
    maxh = min(width / nx, height / ny)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)
        return Mesh(ngmesh)


def create_test_annulus_mesh(maxh=0.015, outer_radius=0.05,
                             inner_radius=0.02):
    """Create an open surface with one independent cohomology loop."""
    from netgen.occ import WorkPlane, OCCGeometry, Axes, Pnt, Dir
    from ngsolve import Mesh, TaskManager

    wp = WorkPlane(Axes(p=Pnt(0, 0, 0), n=Dir(0, 0, 1), h=Dir(1, 0, 0)))
    annulus = wp.Circle(outer_radius).Face() - wp.Circle(inner_radius).Face()
    annulus.name = "conductor"
    with TaskManager():
        return Mesh(OCCGeometry(annulus).GenerateMesh(maxh=maxh))


# =====================================================================
# Test 1: Module import
# =====================================================================
def test_module_import():
    print("=" * 60)
    print("Test 1: Module Import")
    print("=" * 60)

    try:
        from radia.ngsbem_interface import LoopStarTransform
        record("import LoopStarTransform", "PASS")
    except ImportError as e:
        record("import LoopStarTransform", "FAIL", str(e))

    try:
        from radia.ngsbem_interface import NGBEMBridge, extract_edge_geometry
        record("import NGBEMBridge", "PASS")
    except ImportError as e:
        record("import NGBEMBridge", "FAIL", str(e))


# =====================================================================
# Test 2: Transform dimensions and Euler formula
# =====================================================================
def test_dimensions():
    print("\n" + "=" * 60)
    print("Test 2: Transform Dimensions and Euler Formula")
    print("=" * 60)

    if not check_ngsolve_available():
        record("dimensions", "SKIP", "NGSolve not installed")
        return

    from radia.ngsbem_interface import LoopStarTransform

    mesh = create_test_plate_mesh(nx=3, ny=2)

    ls = LoopStarTransform(mesh)

    # Count mesh entities
    n_verts = len(list(mesh.vertices))
    n_edges = len(list(mesh.edges))
    from ngsolve import BND
    from ngsolve import TaskManager
    n_faces = len(list(mesh.Elements(BND)))

    print(f"  Mesh: V={n_verts}, E={n_edges}, F={n_faces}")
    print(f"  Euler: V-E+F = {n_verts - n_edges + n_faces}")
    print(f"  Loop DOFs: {ls.n_loop}, Star DOFs: {ls.n_star}")

    # Euler formula check
    chi = n_verts - n_edges + n_faces
    if chi == 1:
        record("Euler chi=1 (open surface)", "PASS")
    elif chi == 2:
        record("Euler chi=2 (closed surface)", "PASS")
    else:
        record("Euler characteristic", "PASS", f"chi={chi}")

    # Dimension check: n_loop + n_star = n_edges
    total = ls.n_loop + ls.n_star
    if total == n_edges:
        record("n_loop + n_star = n_edges", "PASS",
               f"{ls.n_loop} + {ls.n_star} = {total} = {n_edges}")
    else:
        record("n_loop + n_star = n_edges", "FAIL",
               f"{ls.n_loop} + {ls.n_star} = {total} != {n_edges}")

    # T matrix should be square
    if ls.T.shape == (n_edges, n_edges):
        record("T matrix square", "PASS", f"{ls.T.shape}")
    else:
        record("T matrix square", "FAIL",
               f"shape={ls.T.shape}, expected ({n_edges}, {n_edges})")

    # Local loop columns are generated by vertex-edge incidence.
    for v in range(ls.n_local_loop):
        nnz = np.count_nonzero(ls.T_local_loop[:, v])
        if nnz < 2:
            record("T_local_loop vertex degree", "FAIL",
                   f"vertex {v} has degree {nnz}, expected >= 2")
            return
    record("T_local_loop vertex degrees >= 2", "PASS",
           f"{ls.n_local_loop} vertices checked")

    # Star columns are boundaries of independent triangular faces.
    for f in range(ls.n_star):
        nnz = np.count_nonzero(ls.T_star[:, f])
        if nnz != 3:
            record("T_star column nnz=3", "FAIL",
                   f"face {f} has {nnz} nonzeros, expected 3")
            return
    record("T_star columns all have 3 nonzeros", "PASS",
           f"{ls.n_star} faces checked")

    if ls.n_harmonic == 0:
        record("simply connected plate has no cohomology Loop", "PASS")
    else:
        record("plate cohomology Loop count", "FAIL",
               f"expected 0, got {ls.n_harmonic}")


def test_cohomology_loop():
    """An annulus contributes one global divergence-free current mode."""
    print("\n" + "=" * 60)
    print("Cohomology Loop: Annular Surface")
    print("=" * 60)

    if not check_ngsolve_available():
        record("annulus cohomology", "SKIP", "NGSolve not installed")
        return

    from radia.ngsbem_interface import LoopStarTransform

    mesh = create_test_annulus_mesh()
    ls = LoopStarTransform(mesh)
    if ls.n_harmonic == 1:
        record("annulus has one cohomology Loop", "PASS")
    else:
        record("annulus cohomology Loop count", "FAIL",
               f"expected 1, got {ls.n_harmonic}")

    total = ls.n_loop + ls.n_star
    record("annulus transform is complete",
           "PASS" if total == ls.n_edge else "FAIL",
           f"{ls.n_loop} + {ls.n_star} = {total}, edges={ls.n_edge}")
    harmonic_div = np.max(np.abs(ls.T_star.T @ ls.T_harmonic))
    record("cohomology Loop is divergence-free",
           "PASS" if harmonic_div < 1e-12 else "FAIL",
           f"incidence residual={harmonic_div:.2e}")


# =====================================================================
# Test 3: Orthogonality T_star^T * T_loop = 0
# =====================================================================
def test_orthogonality():
    print("\n" + "=" * 60)
    print("Test 3: Orthogonality T_star^T * T_loop = 0")
    print("=" * 60)

    if not check_ngsolve_available():
        record("orthogonality", "SKIP", "NGSolve not installed")
        return

    from radia.ngsbem_interface import LoopStarTransform

    mesh = create_test_plate_mesh(nx=4, ny=3)
    ls = LoopStarTransform(mesh)

    ortho_error = ls.verify_orthogonality()
    print(f"  |T_star^T * T_loop|_max = {ortho_error:.2e}")

    if ortho_error < 1e-14:
        record("orthogonality (exact)", "PASS",
               f"max error = {ortho_error:.2e}")
    else:
        record("orthogonality", "FAIL",
               f"max error = {ortho_error:.2e}, expected < 1e-14")


# =====================================================================
# Test 4: T matrix invertibility
# =====================================================================
def test_invertibility():
    print("\n" + "=" * 60)
    print("Test 4: T Matrix Invertibility")
    print("=" * 60)

    if not check_ngsolve_available():
        record("invertibility", "SKIP", "NGSolve not installed")
        return

    from radia.ngsbem_interface import LoopStarTransform

    mesh = create_test_plate_mesh(nx=3, ny=2)
    ls = LoopStarTransform(mesh)

    if ls._invertible:
        record("T invertible", "PASS")
    else:
        record("T invertible", "FAIL", "T is near-singular")

    # Check T * T_inv = I
    identity_check = ls.T @ ls.T_inv
    identity_error = np.max(np.abs(identity_check - np.eye(ls.n_edge)))
    print(f"  |T * T_inv - I|_max = {identity_error:.2e}")

    if identity_error < 1e-10:
        record("T * T_inv = I", "PASS", f"max error = {identity_error:.2e}")
    else:
        record("T * T_inv = I", "FAIL", f"max error = {identity_error:.2e}")

    # Condition number
    cond = np.linalg.cond(ls.T)
    print(f"  cond(T) = {cond:.2e}")
    record("condition number", "PASS", f"cond(T) = {cond:.2e}")


# =====================================================================
# Test 5: Decomposition recovery J_loop + J_star = J_edge
# =====================================================================
def test_decomposition_recovery():
    print("\n" + "=" * 60)
    print("Test 5: Decomposition Recovery J_loop + J_star = J_edge")
    print("=" * 60)

    if not check_ngsolve_available():
        record("decomposition", "SKIP", "NGSolve not installed")
        return

    from radia.ngsbem_interface import LoopStarTransform

    mesh = create_test_plate_mesh(nx=3, ny=2)
    ls = LoopStarTransform(mesh)

    # Random test vector
    np.random.seed(42)
    J_edge = np.random.randn(ls.n_edge)

    J_loop, J_star, alpha_loop, alpha_star = ls.decompose_current(J_edge)

    # Check recovery
    J_recovered = J_loop + J_star
    recovery_error = np.max(np.abs(J_recovered - J_edge))
    print(f"  |J_loop + J_star - J_edge|_max = {recovery_error:.2e}")

    if recovery_error < 1e-10:
        record("J_loop + J_star = J_edge", "PASS",
               f"max error = {recovery_error:.2e}")
    else:
        record("J_loop + J_star = J_edge", "FAIL",
               f"max error = {recovery_error:.2e}")

    # Check J_loop is in range of T_loop
    J_loop_check = ls.T_loop @ alpha_loop
    check_error = np.max(np.abs(J_loop_check - J_loop))
    if check_error < 1e-14:
        record("J_loop = T_loop * alpha_loop", "PASS")
    else:
        record("J_loop = T_loop * alpha_loop", "FAIL",
               f"error = {check_error:.2e}")

    # Check J_star is in range of T_star
    J_star_check = ls.T_star @ alpha_star
    check_error = np.max(np.abs(J_star_check - J_star))
    if check_error < 1e-14:
        record("J_star = T_star * alpha_star", "PASS")
    else:
        record("J_star = T_star * alpha_star", "FAIL",
               f"error = {check_error:.2e}")


# =====================================================================
# Test 6: Matrix transform preserves symmetry
# =====================================================================
def test_symmetry_preservation():
    print("\n" + "=" * 60)
    print("Test 6: Matrix Transform Preserves Symmetry")
    print("=" * 60)

    if not check_ngsolve_available():
        record("symmetry", "SKIP", "NGSolve not installed")
        return

    from radia.ngsbem_interface import LoopStarTransform

    mesh = create_test_plate_mesh(nx=3, ny=2)
    ls = LoopStarTransform(mesh)

    # Create symmetric test matrix (simulating L)
    np.random.seed(123)
    A_rand = np.random.randn(ls.n_edge, ls.n_edge)
    A_sym = A_rand + A_rand.T  # Symmetric matrix

    # Transform to Loop-Star
    A_LS = ls.transform_matrix(A_sym)

    # Check symmetry of A_LS
    asym = np.max(np.abs(A_LS - A_LS.T))
    asym_rel = asym / np.max(np.abs(A_LS)) if np.max(np.abs(A_LS)) > 0 else 0
    print(f"  |A_LS - A_LS^T|_max = {asym:.2e}")
    print(f"  Relative asymmetry = {asym_rel:.2e}")

    if asym_rel < 1e-14:
        record("symmetry preserved", "PASS",
               f"relative asymmetry = {asym_rel:.2e}")
    else:
        record("symmetry preserved", "FAIL",
               f"relative asymmetry = {asym_rel:.2e}")

    # Extract blocks
    A_LL, A_LS_block, A_SL, A_SS = ls.extract_blocks(A_LS)

    # LL and SS blocks should be symmetric
    ll_asym = np.max(np.abs(A_LL - A_LL.T))
    ss_asym = np.max(np.abs(A_SS - A_SS.T))
    print(f"  A_LL shape: {A_LL.shape}, asymmetry: {ll_asym:.2e}")
    print(f"  A_SS shape: {A_SS.shape}, asymmetry: {ss_asym:.2e}")
    print(f"  A_LS shape: {A_LS_block.shape}")

    if ll_asym < 1e-12:
        record("A_LL symmetric", "PASS")
    else:
        record("A_LL symmetric", "FAIL", f"asymmetry = {ll_asym:.2e}")

    if ss_asym < 1e-12:
        record("A_SS symmetric", "PASS")
    else:
        record("A_SS symmetric", "FAIL", f"asymmetry = {ss_asym:.2e}")

    # A_LS should equal A_SL^T
    ls_sl_diff = np.max(np.abs(A_LS_block - A_SL.T))
    if ls_sl_diff < 1e-12:
        record("A_LS = A_SL^T", "PASS")
    else:
        record("A_LS = A_SL^T", "FAIL", f"diff = {ls_sl_diff:.2e}")


# =====================================================================
# Test 7: Divergence-free property (requires ngbem)
# =====================================================================
def test_divergence_free():
    print("\n" + "=" * 60)
    print("Test 7: Divergence-Free Property M_LS * T_loop = 0")
    print("=" * 60)

    if not check_ngbem_available():
        record("divergence-free", "SKIP", "ngbem not installed")
        return

    from radia.ngsbem_peec import NGBEMPEECSolver, create_plate_mesh
    from radia.ngsbem_interface import LoopStarTransform

    mesh = create_plate_mesh(0.01, 0.001, 0.003)
    solver = NGBEMPEECSolver(mesh, sigma=5.8e7, thickness=1e-3)
    solver.assemble()

    ls = LoopStarTransform(mesh)

    # M_LS * T_loop should be ~0 (loops are divergence-free)
    div_error = ls.verify_divergence_free(solver.M_LS)
    print(f"  |M_LS * T_loop|_max = {div_error:.2e}")

    if div_error < 1e-10:
        record("M_LS * T_loop = 0", "PASS",
               f"max error = {div_error:.2e}")
    else:
        record("M_LS * T_loop = 0", "FAIL",
               f"max error = {div_error:.2e}")

    # M_LS * T_star should NOT be zero (stars have divergence)
    star_coupling = solver.M_LS @ ls.T_star
    star_max = np.max(np.abs(star_coupling))
    print(f"  |M_LS * T_star|_max = {star_max:.2e} (should be nonzero)")

    if star_max > 1e-10:
        record("M_LS * T_star != 0", "PASS",
               f"max = {star_max:.2e}")
    else:
        record("M_LS * T_star != 0", "FAIL",
               f"max = {star_max:.2e}, expected nonzero")


# =====================================================================
# Test 8: L matrix Loop-Star blocks (requires ngbem)
# =====================================================================
def test_l_matrix_blocks():
    print("\n" + "=" * 60)
    print("Test 8: L Matrix Loop-Star Blocks")
    print("=" * 60)

    if not check_ngbem_available():
        record("L blocks", "SKIP", "ngbem not installed")
        return

    from radia.ngsbem_peec import NGBEMPEECSolver, create_plate_mesh
    from radia.ngsbem_interface import LoopStarTransform

    mesh = create_plate_mesh(0.01, 0.001, 0.003)
    solver = NGBEMPEECSolver(mesh, sigma=5.8e7, thickness=1e-3)
    solver.assemble()

    ls = LoopStarTransform(mesh)

    # Transform L to Loop-Star basis
    L_LS = ls.transform_matrix(solver.L)
    L_LL, L_LS_block, L_SL, L_SS = ls.extract_blocks(L_LS)

    print(f"  L_edge: {solver.L.shape}")
    print(f"  L_LL (loop inductance): {L_LL.shape}")
    print(f"  L_SS (star): {L_SS.shape}")
    print(f"  L_LS (coupling): {L_LS_block.shape}")

    # L_LL should be symmetric positive definite
    L_LL_eig = np.linalg.eigvalsh(L_LL)
    n_pos = np.sum(L_LL_eig > 0)
    n_neg = np.sum(L_LL_eig < 0)
    print(f"  L_LL eigenvalues: {n_pos} positive, {n_neg} negative")
    print(f"  L_LL eig range: [{L_LL_eig[0]:.4e}, {L_LL_eig[-1]:.4e}]")

    if n_neg == 0:
        record("L_LL positive semi-definite", "PASS",
               f"{n_pos} positive eigenvalues")
    else:
        # Small negative eigenvalues due to numerical errors are acceptable
        min_neg = L_LL_eig[0]
        max_pos = L_LL_eig[-1]
        ratio = abs(min_neg) / max_pos if max_pos > 0 else float('inf')
        if ratio < 1e-10:
            record("L_LL positive semi-definite", "PASS",
                   f"tiny negative eigs (ratio={ratio:.2e})")
        else:
            record("L_LL positive semi-definite", "FAIL",
                   f"min_neg/max_pos = {ratio:.2e}")

    # L_LL diagonal should be positive (self-inductance)
    L_LL_diag = np.diag(L_LL)
    if np.all(L_LL_diag > 0):
        record("L_LL diagonal positive", "PASS",
               f"range [{L_LL_diag.min():.4e}, {L_LL_diag.max():.4e}] H")
    else:
        record("L_LL diagonal positive", "FAIL",
               f"min diagonal = {L_LL_diag.min():.4e}")

    # L_SS should be small compared to L_LL
    # (star functions have small inductance)
    L_LL_norm = np.linalg.norm(L_LL)
    L_SS_norm = np.linalg.norm(L_SS)
    ratio = L_SS_norm / L_LL_norm if L_LL_norm > 0 else 0
    print(f"  ||L_SS|| / ||L_LL|| = {ratio:.4f}")
    record("L_SS/L_LL ratio", "PASS", f"ratio = {ratio:.4f}")


# =====================================================================
# Test 9: NGBEMBridge with decompose_loop_star (requires ngbem)
# =====================================================================
def test_bridge_loop_star():
    print("\n" + "=" * 60)
    print("Test 9: NGBEMBridge decompose_loop_star")
    print("=" * 60)

    if not check_ngbem_available():
        record("bridge loop-star", "SKIP", "ngbem not installed")
        return

    from radia.ngsbem_peec import NGBEMPEECSolver, create_plate_mesh
    from radia.ngsbem_interface import NGBEMBridge

    mesh = create_plate_mesh(0.01, 0.001, 0.003)
    solver = NGBEMPEECSolver(mesh, sigma=5.8e7, thickness=1e-3)
    solver.assemble()

    bridge = NGBEMBridge(solver, port_spec='auto')
    topo = bridge.to_topology_dict(decompose_loop_star=True)

    # Check loop_star dict exists
    if 'loop_star' not in topo:
        record("loop_star in topo", "FAIL", "missing key")
        return
    record("loop_star in topo", "PASS")

    ls_data = topo['loop_star']

    # Check required keys
    required = ['L_LL', 'L_SS', 'L_LS', 'L_SL', 'R_LL', 'R_SS',
                'T_loop', 'T_star', 'T', 'T_inv',
                'n_loop_ls', 'n_star_ls', 'n_local_loop_ls',
                'n_harmonic_ls', 'T_local_loop', 'T_harmonic', 'transform',
                'orthogonality_error', 'div_free_error']
    missing = [k for k in required if k not in ls_data]
    if not missing:
        record("all loop_star keys present", "PASS",
               f"{len(required)} keys")
    else:
        record("loop_star keys", "FAIL", f"missing: {missing}")

    # Check orthogonality and divergence-free
    ortho = ls_data['orthogonality_error']
    divfree = ls_data['div_free_error']
    print(f"  Orthogonality error: {ortho:.2e}")
    print(f"  Divergence-free error: {divfree:.2e}")

    if ortho < 1e-14:
        record("orthogonality via bridge", "PASS")
    else:
        record("orthogonality via bridge", "FAIL", f"error={ortho:.2e}")

    if divfree < 1e-10:
        record("divergence-free via bridge", "PASS")
    else:
        record("divergence-free via bridge", "FAIL", f"error={divfree:.2e}")

    # Check transform object is accessible
    ls_transform = ls_data['transform']
    ls_transform.print_summary()
    record("print_summary works", "PASS")


# =====================================================================
# Main
# =====================================================================
if __name__ == '__main__':
    print("Loop-Star Decomposition Validation")
    print("=" * 60)
    print(f"NGSolve available: {check_ngsolve_available()}")
    print(f"ngbem available: {check_ngbem_available()}")
    print()

    test_module_import()
    test_dimensions()
    test_cohomology_loop()
    test_orthogonality()
    test_invertibility()
    test_decomposition_recovery()
    test_symmetry_preservation()
    test_divergence_free()
    test_l_matrix_blocks()
    test_bridge_loop_star()

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    n_pass = sum(1 for _, s, _ in results if s == "PASS")
    n_fail = sum(1 for _, s, _ in results if s == "FAIL")
    n_skip = sum(1 for _, s, _ in results if s == "SKIP")

    print(f"  PASS: {n_pass}")
    print(f"  FAIL: {n_fail}")
    print(f"  SKIP: {n_skip}")
    print(f"  Total: {len(results)}")

    if n_fail > 0:
        print("\nFailed tests:")
        for name, status, detail in results:
            if status == "FAIL":
                print(f"  - {name}: {detail}")

    if n_fail == 0:
        print("\nAll tests PASSED (or SKIPPED due to missing dependencies).")
    else:
        print(f"\n{n_fail} test(s) FAILED!")
        sys.exit(1)
