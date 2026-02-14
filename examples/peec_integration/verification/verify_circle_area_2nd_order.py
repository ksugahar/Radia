"""
verify_circle_area_2nd_order.py

Compare circle area accuracy between:
  1. GMSH 2nd order elements (Tri6) with CAD geometry
  2. Netgen 1st order mesh + mesh.Curve(2) with OCC geometry

Both methods should achieve 2nd order geometric accuracy for curved boundaries.
The exact area is pi * a^2.

Part of Radia project
"""

import sys
import os
import numpy as np

# ============================================================
# Method 1: GMSH 2nd order
# ============================================================

def gmsh_circle_area(radius, mesh_size, order=2):
    """
    Create circle mesh with GMSH, compute area from element integration.

    Args:
        radius: Circle radius
        mesh_size: Characteristic mesh size
        order: Element order (1 or 2)

    Returns:
        dict with area, n_elements, error
    """
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.model.add('circle')

    # Create circle disk
    gmsh.model.occ.addDisk(0, 0, 0, radius, radius)
    gmsh.model.occ.synchronize()

    # Set mesh size
    gmsh.option.setNumber('Mesh.CharacteristicLengthMax', mesh_size)
    gmsh.option.setNumber('Mesh.CharacteristicLengthMin', mesh_size * 0.5)

    # Generate 2D mesh
    gmsh.model.mesh.generate(2)

    # Set element order
    if order > 1:
        gmsh.model.mesh.setOrder(order)

    # Get element data
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)

    # Get node coordinates
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_map = {}
    for i, tag in enumerate(node_tags):
        node_map[tag] = np.array([node_coords[3*i], node_coords[3*i+1], node_coords[3*i+2]])

    total_area = 0.0
    n_elements = 0

    for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
        # Get element properties
        name, dim, eorder, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(etype)

        if dim != 2:
            continue

        n_elem = len(etags)
        n_elements += n_elem

        for i in range(n_elem):
            nodes_idx = enodes[i * n_nodes:(i + 1) * n_nodes]

            if order == 1:
                # Linear triangle: simple area formula
                if n_nodes == 3:
                    v0 = node_map[nodes_idx[0]]
                    v1 = node_map[nodes_idx[1]]
                    v2 = node_map[nodes_idx[2]]
                    area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
                    total_area += area
            else:
                # 2nd order: use Jacobian-based integration
                if n_nodes == 6:
                    # Tri6: 6-node quadratic triangle
                    # Use 3-point Gauss quadrature on reference triangle
                    total_area += _tri6_area(
                        [node_map[nodes_idx[j]] for j in range(6)])

    gmsh.finalize()

    exact_area = np.pi * radius**2
    error = abs(total_area - exact_area) / exact_area

    return {
        'area': total_area,
        'exact': exact_area,
        'error': error,
        'n_elements': n_elements,
        'order': order,
        'method': f'GMSH order={order}',
    }


def _tri6_area(verts):
    """
    Compute area of a 6-node quadratic triangle using Gauss quadrature.

    Nodes: 0,1,2 = corners; 3 = mid(0,1); 4 = mid(1,2); 5 = mid(2,0)
    Reference triangle: (xi, eta) in [0,1] x [0,1], xi+eta <= 1

    Shape functions (GMSH Tri6 ordering):
      N0 = (1-xi-eta)(1-2*xi-2*eta)  [= lambda0*(2*lambda0-1)]
      N1 = xi*(2*xi-1)
      N2 = eta*(2*eta-1)
      N3 = 4*xi*(1-xi-eta)
      N4 = 4*xi*eta
      N5 = 4*eta*(1-xi-eta)
    """
    # 4-point Gauss quadrature on triangle (degree 3, exact for quadratics)
    gauss_pts = [
        (1.0/3, 1.0/3, -27.0/96),
        (0.6, 0.2, 25.0/96),
        (0.2, 0.6, 25.0/96),
        (0.2, 0.2, 25.0/96),
    ]

    area = 0.0
    for xi, eta, w in gauss_pts:
        # Shape function derivatives w.r.t. xi and eta
        lam0 = 1.0 - xi - eta

        dN_dxi = np.array([
            -1.0*(2*lam0 - 1) + lam0*(-2),    # dN0/dxi
            2*xi - 1 + xi*2,                     # dN1/dxi  (= 4*xi - 1)
            0.0,                                  # dN2/dxi
            4*(lam0 - xi),                        # dN3/dxi  (= 4*(1-2*xi-eta))
            4*eta,                                # dN4/dxi
            -4*eta,                               # dN5/dxi
        ])

        dN_deta = np.array([
            -1.0*(2*lam0 - 1) + lam0*(-2),    # dN0/deta
            0.0,                                  # dN1/deta
            2*eta - 1 + eta*2,                    # dN2/deta (= 4*eta - 1)
            -4*xi,                                # dN3/deta
            4*xi,                                 # dN4/deta
            4*(lam0 - eta),                       # dN5/deta (= 4*(1-xi-2*eta))
        ])

        # Jacobian: dx/dxi, dx/deta, dy/dxi, dy/deta
        dx_dxi = sum(dN_dxi[j] * verts[j][0] for j in range(6))
        dy_dxi = sum(dN_dxi[j] * verts[j][1] for j in range(6))
        dx_deta = sum(dN_deta[j] * verts[j][0] for j in range(6))
        dy_deta = sum(dN_deta[j] * verts[j][1] for j in range(6))

        det_J = abs(dx_dxi * dy_deta - dx_deta * dy_dxi)
        area += w * det_J

    # Factor of 2 because reference triangle has area 0.5
    # but Gauss weights already account for this (sum of weights = 0.5)
    return area


# ============================================================
# Method 2: Netgen + mesh.Curve(2)
# ============================================================

def netgen_circle_area(radius, mesh_size, curve_order=2):
    """
    Create circle mesh with Netgen, apply mesh.Curve(), compute area.

    Args:
        radius: Circle radius
        mesh_size: Characteristic mesh size
        curve_order: Order for mesh.Curve()

    Returns:
        dict with area, n_elements, error
    """
    from netgen.occ import WorkPlane, OCCGeometry
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh, Integrate, CoefficientFunction

    # Create circle as 2D OCC geometry (two semicircles)
    wp = WorkPlane()
    wp.MoveTo(radius, 0).Arc(radius, 180).Arc(radius, 180)
    face = wp.Face()
    geo = OCCGeometry(face, dim=2)

    # Generate mesh
    mp = MeshingParameters(maxh=mesh_size)
    ngmesh = geo.GenerateMesh(mp)
    mesh = Mesh(ngmesh)

    # Apply curving
    if curve_order > 1:
        mesh.Curve(curve_order)

    # Compute area using NGSolve integration
    area = Integrate(CoefficientFunction(1.0), mesh)
    n_elements = mesh.ne

    exact_area = np.pi * radius**2
    error = abs(area - exact_area) / exact_area

    return {
        'area': area,
        'exact': exact_area,
        'error': error,
        'n_elements': n_elements,
        'order': curve_order,
        'method': f'Netgen Curve({curve_order})',
    }


# ============================================================
# Method 3: Netgen without Curve (1st order reference)
# ============================================================

def netgen_circle_area_linear(radius, mesh_size):
    """Netgen 1st order (no Curve) for reference."""
    return netgen_circle_area(radius, mesh_size, curve_order=1)


# ============================================================
# Main comparison
# ============================================================

def main():
    radius = 1.0
    exact_area = np.pi * radius**2

    print()
    print("Circle Area Accuracy: GMSH 2nd Order vs Netgen Curve(2)")
    print("=" * 70)
    print(f"  Radius = {radius}, Exact area = pi = {exact_area:.10f}")
    print()

    # Test various mesh sizes
    mesh_sizes = [0.8, 0.5, 0.3, 0.2, 0.15, 0.1]

    # Check which libraries are available
    has_gmsh = False
    has_netgen = False

    try:
        import gmsh
        has_gmsh = True
    except ImportError:
        print("  WARNING: gmsh not available, skipping GMSH tests")

    try:
        from ngsolve import Mesh
        has_netgen = True
    except ImportError:
        print("  WARNING: ngsolve not available, skipping Netgen tests")

    if not has_gmsh and not has_netgen:
        print("  ERROR: Neither gmsh nor ngsolve available. Cannot run test.")
        return

    print()

    # Collect results
    results_gmsh_1 = []
    results_gmsh_2 = []
    results_netgen_1 = []
    results_netgen_2 = []

    for h in mesh_sizes:
        if has_gmsh:
            try:
                r1 = gmsh_circle_area(radius, h, order=1)
                results_gmsh_1.append(r1)
            except Exception as e:
                print(f"  GMSH order=1 h={h}: FAILED ({e})")

            try:
                r2 = gmsh_circle_area(radius, h, order=2)
                results_gmsh_2.append(r2)
            except Exception as e:
                print(f"  GMSH order=2 h={h}: FAILED ({e})")

        if has_netgen:
            try:
                r3 = netgen_circle_area_linear(radius, h)
                results_netgen_1.append(r3)
            except Exception as e:
                print(f"  Netgen Curve(1) h={h}: FAILED ({e})")

            try:
                r4 = netgen_circle_area(radius, h, curve_order=2)
                results_netgen_2.append(r4)
            except Exception as e:
                print(f"  Netgen Curve(2) h={h}: FAILED ({e})")

    # Print results table
    print("-" * 70)
    print(f"  {'Method':<20s}  {'N_elem':>7s}  {'Area':>14s}  {'Error':>12s}")
    print("-" * 70)

    all_groups = [
        ("GMSH order=1", results_gmsh_1),
        ("GMSH order=2", results_gmsh_2),
        ("Netgen Curve(1)", results_netgen_1),
        ("Netgen Curve(2)", results_netgen_2),
    ]

    for label, results in all_groups:
        if not results:
            continue
        print(f"\n  {label}:")
        for r in results:
            print(f"    N={r['n_elements']:5d}  "
                  f"Area={r['area']:.10f}  "
                  f"Error={r['error']:.2e}")

    # Convergence comparison
    print()
    print("=" * 70)
    print("Convergence Summary")
    print("=" * 70)
    print()
    print(f"  {'Method':<20s}  {'N_min':>6s} {'Err_min':>10s}  "
          f"{'N_max':>6s} {'Err_max':>10s}  {'Rate':>6s}")
    print(f"  {'-'*20}  {'-'*6} {'-'*10}  {'-'*6} {'-'*10}  {'-'*6}")

    for label, results in all_groups:
        if len(results) < 2:
            continue
        r_first = results[0]
        r_last = results[-1]

        # Estimate convergence rate from first and last
        if r_first['error'] > 0 and r_last['error'] > 0:
            h_ratio = np.sqrt(r_first['n_elements'] / r_last['n_elements'])
            if h_ratio > 0 and h_ratio != 1:
                rate = np.log(r_first['error'] / r_last['error']) / np.log(1.0 / h_ratio)
            else:
                rate = 0
        else:
            rate = 0

        print(f"  {label:<20s}  {r_first['n_elements']:6d} {r_first['error']:10.2e}  "
              f"{r_last['n_elements']:6d} {r_last['error']:10.2e}  {rate:6.2f}")

    print()
    print("  Expected convergence rates:")
    print("    1st order geometry: O(h^2) -> rate ~= 2")
    print("    2nd order geometry: O(h^4) -> rate ~= 4  (or better)")
    print()


if __name__ == '__main__':
    main()
