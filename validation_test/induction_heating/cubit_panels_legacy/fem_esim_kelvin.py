"""
FEM prototype: A-formulation + Kelvin + ESIM Robin BC (2D axisymmetric).

Cross-check against BEM + ESIM (impedance_esim.py).

Two modes:
  A) FEM-full:  workpiece volume meshed with sigma (reference)
  B) FEM-ESIM:  workpiece replaced by ESIM Robin BC on its surface

Both use Kelvin transformation for open boundary (no PML).

2D axisymmetric A-formulation:
  Variable: u = r * A_theta
  Equation: -div(nu/r * grad(u)) + j*omega*sigma*u/r = J_theta  (in workpiece)
  Kelvin: nu_K = nu0 * (rho'/a)^2 in exterior domain

Usage:
    python fem_esim_kelvin.py                  # Both modes, compare
    python fem_esim_kelvin.py --mode full      # FEM-full only
    python fem_esim_kelvin.py --mode esim      # FEM-ESIM only
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from ngsolve import *
from ngsolve import TaskManager
from netgen.occ import *

MU_0 = 4e-7 * np.pi
NU_0 = 1.0 / MU_0


# ============================================================
# Geometry: air + coil + workpiece + Kelvin exterior
# ============================================================
def build_geometry(R_coil, a_coil, R_wp, H_wp, a_kelvin, z_offset,
                   include_workpiece_volume=True, maxh=0.003,
                   coil_shape='square', curve_order=2,
                   maxh_wp_bnd=None, maxh_coil=None):
    """Build 2D axisymmetric geometry.

    Args:
        R_coil: Coil center radius [m]
        a_coil: Coil wire radius [m]
        R_wp: Workpiece cylinder radius [m]
        H_wp: Workpiece cylinder height [m]
        a_kelvin: Kelvin boundary radius [m]
        z_offset: Z-offset for exterior domain [m]
        include_workpiece_volume: If True, mesh workpiece interior (FEM-full)
        maxh: Mesh size [m]
        coil_shape: 'square' (2a x 2a, cross-section area = 4*a^2) or
            'circular' (radius a, cross-section area = pi*a^2). Use
            'circular' to match a 3D circular-torus coil exactly.

    Returns:
        NGSolve Mesh
    """
    # Interior domain: half-circle
    wp_inner = WorkPlane()
    inner_full = wp_inner.Circle(a_kelvin).Face()

    # Cut to half (x >= 0)
    cutter = MoveTo(-a_kelvin - 0.1, -a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    inner_half = inner_full - cutter

    # Coil cross-section
    if coil_shape == 'circular':
        # Circle of radius a_coil centered at (R_coil, 0)
        wp_coil = WorkPlane()
        wp_coil.MoveTo(R_coil + a_coil, 0)
        wp_coil.Direction(0, 1)  # heading UP so arc center = (R_coil, 0)
        wp_coil.Arc(a_coil, 180).Arc(a_coil, 180)
        coil = wp_coil.Face()
    else:
        # Square 2a x 2a centered at (R_coil, 0)
        w = a_coil * 2
        coil = MoveTo(R_coil - w / 2, -w / 2).Rectangle(w, w).Face()
    coil.name = "coil"
    coil.maxh = maxh_coil if maxh_coil is not None else maxh / 2

    # Workpiece cross-section (rectangle in r-z plane)
    wp_rect = MoveTo(0, -H_wp / 2).Rectangle(R_wp, H_wp).Face()
    # Name workpiece edges BEFORE boolean (names survive subtraction)
    for edge in wp_rect.edges:
        cx, cy = edge.center.x, edge.center.y
        if cx < 0.001:
            edge.name = "axis"  # left edge on axis
        else:
            edge.name = "workpiece_bnd"
            if maxh_wp_bnd is not None:
                edge.maxh = maxh_wp_bnd

    if include_workpiece_volume:
        # FEM-full: workpiece is a meshed region
        wp_rect.name = "workpiece"
        wp_rect.maxh = maxh / 3  # finer mesh for skin effect

        # Air = inner_half - coil - workpiece
        air = inner_half - coil - wp_rect
        air.name = "air_inner"
        shapes = [air, coil, wp_rect]
    else:
        # FEM-ESIM: workpiece is a hole (not meshed)
        air = inner_half - coil - wp_rect
        air.name = "air_inner"
        shapes = [air, coil]

    # Exterior domain (Kelvin)
    wp_ext = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
    outer_full = wp_ext.Circle(a_kelvin).Face()
    cutter_ext = MoveTo(-a_kelvin - 0.1, z_offset - a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    outer_half = outer_full - cutter_ext
    outer_half.name = "air_outer"

    # GND point
    gnd = Vertex(Pnt(0, z_offset, 0))
    gnd.name = "GND"

    # Name edges - detect workpiece boundary by vertex positions
    for edge in air.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = math.sqrt(v0.p.x**2 + v0.p.y**2)
            d1 = math.sqrt(v1.p.x**2 + v1.p.y**2)
            is_arc = abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01

            # Workpiece boundary: edge on the workpiece rectangle with at least
            # one vertex at x > 0 (exclude pure axis edges)
            on_wp_rect = False
            if not include_workpiece_volume and max(v0.p.x, v1.p.x) > 1e-6:
                # Check if vertices lie on workpiece rectangle boundary
                def on_rect(px, py):
                    on_right = abs(px - R_wp) < 1e-4 and -H_wp/2 - 1e-4 <= py <= H_wp/2 + 1e-4
                    on_top = abs(py - H_wp/2) < 1e-4 and -1e-4 <= px <= R_wp + 1e-4
                    on_bot = abs(py + H_wp/2) < 1e-4 and -1e-4 <= px <= R_wp + 1e-4
                    return on_right or on_top or on_bot
                on_wp_rect = on_rect(v0.p.x, v0.p.y) and on_rect(v1.p.x, v1.p.y)
        except:
            is_arc = False
            on_wp_rect = False

        if on_wp_rect:
            edge.name = "workpiece_bnd"
        elif is_arc:
            edge.name = "kelvin_int"
        elif cx < 1e-4:
            edge.name = "axis"
        else:
            edge.name = "default"

    for edge in coil.edges:
        edge.name = "coil_bnd" if edge.center.x > 0.01 else "axis"

    if include_workpiece_volume:
        for edge in wp_rect.edges:
            edge.name = "wp_bnd" if edge.center.x > 0.001 else "axis"

    # Kelvin exterior edges
    kelvin_int_edges = []
    kelvin_ext_edges = []

    for edge in air.edges:
        if hasattr(edge, 'name') and edge.name == "kelvin_int":
            kelvin_int_edges.append(edge)

    for edge in outer_half.edges:
        cx = edge.center.x
        cy = edge.center.y - z_offset
        try:
            v0, v1 = edge.vertices
            d0 = math.sqrt(v0.p.x**2 + (v0.p.y - z_offset)**2)
            d1 = math.sqrt(v1.p.x**2 + (v1.p.y - z_offset)**2)
            is_arc = abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01
        except:
            is_arc = False

        if cx < 0.01:
            edge.name = "axis_ext"
        elif is_arc:
            edge.name = "kelvin_ext"
            kelvin_ext_edges.append(edge)
        else:
            edge.name = "default"

    # Periodic identification
    for int_e in kelvin_int_edges:
        for ext_e in kelvin_ext_edges:
            iy = int_e.center.y
            ey = ext_e.center.y - z_offset
            if (iy > 0 and ey > 0) or (iy < 0 and ey < 0):
                int_e.Identify(ext_e, "kelvin", IdentificationType.PERIODIC)
                break

    shape = Glue(shapes + [outer_half, gnd])
    geo = OCCGeometry(shape, dim=2)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh, grading=0.4))
    # Curve the mesh so circular coil / Kelvin arcs are represented
    # correctly at high FES order. Without this, circles are polygonal
    # and L is systematically under-predicted by several percent.
    if curve_order >= 2:
        mesh.Curve(curve_order)
    return mesh


# ============================================================
# FEM-full: workpiece with sigma (direct eddy current solve)
# ============================================================
def solve_fem_full(R_coil, a_coil, R_wp, H_wp, sigma, frequency,
                   I_total=1.0, a_kelvin=0.15, z_offset=0.4,
                   maxh=0.002, order=3, coil_shape='square'):
    """Solve frequency-domain eddy current with meshed workpiece.

    Args:
        coil_shape: 'square' (2a x 2a) or 'circular' (radius a).
            Use 'circular' to match 3D circular-torus coil.

    Returns dict with P_total, Q_total, L, mesh info.
    """
    omega = 2 * np.pi * frequency
    if coil_shape == 'circular':
        J0 = I_total / (np.pi * a_coil ** 2)
    else:
        J0 = I_total / (2 * a_coil) ** 2

    print(f"  Building geometry (workpiece meshed, coil={coil_shape})...")
    t0 = time.perf_counter()
    mesh = build_geometry(R_coil, a_coil, R_wp, H_wp, a_kelvin, z_offset,
                          include_workpiece_volume=True, maxh=maxh,
                          coil_shape=coil_shape, curve_order=order)
    t_mesh = time.perf_counter() - t0

    materials = mesh.GetMaterials()
    boundaries = mesh.GetBoundaries()
    print(f"  Mesh: {mesh.ne} elements, materials={materials}")
    print(f"  Boundaries: {boundaries}")

    # FE space
    fes_base = H1(mesh, order=order, complex=True,
                  dirichlet="axis|axis_ext", dirichlet_bbnd="GND")
    fes = Periodic(fes_base)
    u, v = fes.TnT()
    print(f"  DOFs: {fes.ndof}")

    # Coefficient functions
    r_cf = IfPos(x - 1e-10, x, 1e-10)

    # Kelvin factor for exterior
    y_local = y - z_offset
    rho_prime = sqrt(x**2 + y_local**2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    kelvin_fac = (rho_safe / a_kelvin)**2

    # Reluctivity per material
    nu_dict = {}
    for mat in materials:
        if "outer" in mat.lower():
            nu_dict[mat] = NU_0 * kelvin_fac
        else:
            nu_dict[mat] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # Conductivity (only in workpiece)
    sigma_dict = {}
    for mat in materials:
        if "workpiece" in mat.lower():
            sigma_dict[mat] = sigma
        else:
            sigma_dict[mat] = 0.0
    sigma_cf = mesh.MaterialCF(sigma_dict, default=0.0)

    # Bilinear form
    a_bf = BilinearForm(fes)
    a_bf += nu_cf / r_cf * grad(u) * grad(v) * dx
    # Eddy current term: j*omega*sigma * (u/r)*(v/r) * r = j*omega*sigma * u*v/r
    a_bf += 1j * omega * sigma_cf * u * v / r_cf * dx

    # Source
    f_lf = LinearForm(fes)
    f_lf += J0 * v * dx("coil")

    print("  Assembling...")
    t0 = time.perf_counter()
    with TaskManager():
        a_bf.Assemble()
        f_lf.Assemble()
        t_asm = time.perf_counter() - t0

        print("  Solving...")
        t0 = time.perf_counter()
        gfu = GridFunction(fes)
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f_lf.vec
        t_solve = time.perf_counter() - t0

        # Post-process: power loss in workpiece
        # P = 0.5 * sigma * omega^2 * |A_theta|^2 * volume
        # In 2D axi: P = 2*pi * integral( 0.5 * sigma * omega^2 * |u/r|^2 * r dr dz )
        #           = 2*pi * integral( 0.5 * sigma * omega^2 * |u|^2 / r dr dz )
        P_density = 0.5 * sigma_cf * omega**2 * Conj(gfu) * gfu / r_cf
        P_per_radian = Integrate(P_density, mesh, definedon=mesh.Materials("workpiece")).real
        P_total = 2 * np.pi * P_per_radian

        # Reactive power (stored magnetic energy rate)
        # Q = 2*pi * integral( 0.5 * omega * nu * |grad(u)|^2 / r )
        # Only meaningful for the workpiece region
        Q_density = 0.5 * omega * sigma_cf * Conj(gfu) * gfu / r_cf
        Q_per_radian = Integrate(Q_density, mesh, definedon=mesh.Materials("workpiece")).real
        Q_total = 2 * np.pi * Q_per_radian

        # Inductance from stored energy
        # W_mag = 0.5 * L * I^2 = 2*pi * integral( 0.5 * nu * |grad(u)|^2 / r )
        W_per_radian = Integrate(
            0.5 * nu_cf / r_cf * (grad(gfu) * Conj(grad(gfu))),
            mesh).real
        W_mag = 2 * np.pi * W_per_radian
        L = 2 * W_mag / (I_total**2)

        return {
            'mode': 'FEM-full',
            'P_total': float(P_total),
            'Q_total': float(Q_total),
            'L': float(L),
            'ndof': fes.ndof,
            'ne': mesh.ne,
            't_mesh': t_mesh,
            't_asm': t_asm,
            't_solve': t_solve,
        }


# ============================================================
# FEM-ESIM one-way: static FEM (no hole) + ESIM per panel
# ============================================================
def solve_fem_esim_oneway(R_coil, a_coil, R_wp, H_wp, sigma, frequency,
                          bh_curve=None, mu_r=1.0,
                          I_total=1.0, a_kelvin=0.15, z_offset=0.4,
                          maxh=0.003, order=3):
    """One-way coupled FEM-ESIM: static coil field + ESIM power.

    Equivalent to BEM-ESIM but using FEM for the coil field.
    No workpiece hole in the mesh (avoids PMC artifact).

    Pipeline:
      1. Build mesh WITHOUT workpiece hole (air + coil + Kelvin)
      2. Solve static A-formulation (no eddy currents)
      3. Sample H_incident at workpiece surface location
      4. ESIM per panel -> P

    Returns dict with P_total, L, H_t distribution.
    """
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * np.pi * frequency
    J0 = I_total / (2 * a_coil)**2

    print("  Building geometry (no workpiece hole)...")
    t0 = time.perf_counter()
    # Build geometry WITHOUT workpiece (include_workpiece_volume=True
    # but we skip the workpiece subtraction by building a simpler mesh)
    inner = WorkPlane().Circle(a_kelvin).Face()
    cutter = MoveTo(-a_kelvin - 0.1, -a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    inner = inner - cutter
    w = 2 * a_coil
    coil = MoveTo(R_coil - w / 2, -w / 2).Rectangle(w, w).Face()
    coil.name = "coil"
    coil.maxh = maxh / 2
    air = inner - coil
    air.name = "air_inner"

    # Kelvin exterior
    from netgen.occ import Axes, Z, X, Vertex, Pnt, IdentificationType
    outer = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X)).Circle(a_kelvin).Face()
    cut_ext = MoveTo(-a_kelvin - 0.1, z_offset - a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    outer = outer - cut_ext
    outer.name = "air_outer"
    gnd = Vertex(Pnt(0, z_offset, 0))
    gnd.name = "GND"

    ki, ke = [], []
    for e in air.edges:
        cx = e.center.x
        try:
            v0, v1 = e.vertices
            d0 = math.sqrt(v0.p.x**2 + v0.p.y**2)
            d1 = math.sqrt(v1.p.x**2 + v1.p.y**2)
            is_arc = abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01
        except:
            is_arc = False
        if is_arc:
            e.name = "kelvin_int"; ki.append(e)
        elif cx < 1e-4:
            e.name = "axis"
        else:
            e.name = "default"
    for e in coil.edges:
        e.name = "coil_bnd" if e.center.x > 0.01 else "axis"
    for e in outer.edges:
        cx = e.center.x
        try:
            v0, v1 = e.vertices
            d0 = math.sqrt(v0.p.x**2 + (v0.p.y - z_offset)**2)
            d1 = math.sqrt(v1.p.x**2 + (v1.p.y - z_offset)**2)
            is_arc = abs(d0 - a_kelvin) < 0.01 and abs(d1 - a_kelvin) < 0.01 and cx > 0.01
        except:
            is_arc = False
        if cx < 0.01:
            e.name = "axis_ext"
        elif is_arc:
            e.name = "kelvin_ext"; ke.append(e)
        else:
            e.name = "default"
    for ie in ki:
        for ee in ke:
            iy = ie.center.y
            ey = ee.center.y - z_offset
            if (iy > 0 and ey > 0) or (iy < 0 and ey < 0):
                ie.Identify(ee, "kelvin", IdentificationType.PERIODIC)
                break

    shape = Glue([air, coil, outer, gnd])
    mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh, grading=0.4))
    t_mesh = time.perf_counter() - t0

    materials = mesh.GetMaterials()
    print(f"  Mesh: {mesh.ne} elements")

    r_cf = IfPos(x - 1e-10, x, 1e-10)
    y_local = y - z_offset
    rho_prime = sqrt(x**2 + y_local**2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    kelvin_fac = (rho_safe / a_kelvin)**2
    nu_dict = {m: NU_0 * kelvin_fac if "outer" in m.lower() else NU_0
               for m in materials}
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # Static solve (real, no Robin BC, no eddy currents)
    fes = Periodic(H1(mesh, order=order, complex=False,
                      dirichlet="axis|axis_ext", dirichlet_bbnd="GND"))
    u, v = fes.TnT()
    bf = BilinearForm(fes)
    bf += nu_cf / r_cf * grad(u) * grad(v) * dx
    bf.Assemble()
    lf = LinearForm(fes)
    lf += J0 * v * dx("coil")
    lf.Assemble()

    t0_s = time.perf_counter()
    gfu = GridFunction(fes)
    gfu.vec.data = bf.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * lf.vec
    t_solve = time.perf_counter() - t0_s

    # Inductance
    W = 2 * np.pi * Integrate(
        0.5 * nu_cf / r_cf * grad(gfu) * grad(gfu), mesh).real
    L = 2 * W / I_total**2

    # ESIM per panel (sample H_incident at workpiece surface location)
    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None, n_nodes=200,
        geometry='cylinder')

    grad_u = grad(gfu)
    n_z_panel = 24
    n_r_cap = 6
    P_total = 0.0
    Q_total = 0.0
    H_t_vals = []

    # Lateral panels
    for j in range(n_z_panel):
        z_j = -H_wp / 2 + (j + 0.5) * H_wp / n_z_panel
        dA = 2 * np.pi * R_wp * (H_wp / n_z_panel)
        try:
            mip = mesh(R_wp, z_j)  # ON the surface (no hole, no offset needed)
            H_z = abs(float(grad_u[0](mip)) / (MU_0 * R_wp))
        except:
            H_z = 0.0
        H_t_vals.append(H_z)
        sol_p = esim_solver.solve(max(H_z, 1e-3))
        P_total += sol_p['P_prime'] * dA
        Q_total += sol_p['Q_prime'] * dA

    # Cap panels
    for sign in [+1, -1]:
        for k in range(n_r_cap):
            r_k = (k + 0.5) * R_wp / n_r_cap
            if r_k < 1e-6:
                continue
            dA = 2 * np.pi * r_k * (R_wp / n_r_cap)
            try:
                mip = mesh(r_k, sign * H_wp / 2)
                H_r = abs(float(grad_u[1](mip)) / (MU_0 * r_k))
            except:
                H_r = 0.0
            H_t_vals.append(H_r)
            sol_p = esim_solver.solve(max(H_r, 1e-3))
            P_total += sol_p['P_prime'] * dA
            Q_total += sol_p['Q_prime'] * dA

    H_t_arr = np.array(H_t_vals)
    H_t_avg = float(np.mean(H_t_arr[H_t_arr > 0])) if np.any(H_t_arr > 0) else 0

    print(f"  P (one-way ESIM) = {P_total:.6e} W")
    print(f"  L (static)       = {L*1e9:.2f} nH")
    print(f"  H_t range        = [{H_t_arr.min():.2f}, {H_t_arr.max():.2f}] A/m")

    return {
        'mode': 'FEM-ESIM one-way',
        'P_total': float(P_total),
        'Q_total': float(Q_total),
        'L': float(L),
        'H_t_avg': float(H_t_avg),
        'ndof': fes.ndof,
        'ne': mesh.ne,
        't_mesh': t_mesh,
        't_solve': t_solve,
    }


# ============================================================
# FEM-ESIM: workpiece as Robin BC (two-way coupled)
# ============================================================
def solve_fem_sibc(R_coil, a_coil, R_wp, H_wp, sigma, frequency,
                   bh_curve=None, mu_r=1.0,
                   I_total=1.0, a_kelvin=0.15, z_offset=0.4,
                   maxh=0.003, order=3, max_iter=15, tol=1e-4,
                   coil_shape='square', maxh_wp_bnd=None, maxh_coil=None):
    """Solve with SIBC Robin BC + Karl Hollaus iteration.

    Karl iteration (Hollaus, IEEE Trans. Mag., 2025):
      1. Solve FEM with current Z_s as Robin BC on workpiece surface
      2. Compute H_t from solution (gradient sampling near boundary)
      3. Update Z_s from ESIM cell problem at average H_t
      4. Repeat until Z_s converges

    Robin BC weak form (2D axi, u = r*A_theta):
      int (nu/r) grad(u).grad(v) dx - (jw/Zs) int u*v/r ds = int J*v dx

    Power (Poynting flux, correct factor = pi not 2*pi):
      P = pi * w^2 * Re(Zs)/|Zs|^2 * int |u|^2/r ds

    Args:
        coil_shape: 'square' (2a x 2a) or 'circular' (radius a).
            Use 'circular' to match 3D circular-torus coil.

    Returns dict with P_total, L, iteration info.
    """
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * np.pi * frequency
    if coil_shape == 'circular':
        J0 = I_total / (np.pi * a_coil ** 2)
    else:
        J0 = I_total / (2 * a_coil) ** 2

    print(f"  Building geometry (workpiece = hole, coil={coil_shape})...")
    t0 = time.perf_counter()
    mesh = build_geometry(R_coil, a_coil, R_wp, H_wp, a_kelvin, z_offset,
                          include_workpiece_volume=False, maxh=maxh,
                          coil_shape=coil_shape, curve_order=order,
                          maxh_wp_bnd=maxh_wp_bnd, maxh_coil=maxh_coil)
    t_mesh = time.perf_counter() - t0

    materials = mesh.GetMaterials()
    print(f"  Mesh: {mesh.ne} elements")

    r_cf = IfPos(x - 1e-10, x, 1e-10)

    # Kelvin
    y_local = y - z_offset
    rho_prime = sqrt(x**2 + y_local**2)
    rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
    kelvin_fac = (rho_safe / a_kelvin)**2
    nu_dict = {m: NU_0 * kelvin_fac if "outer" in m.lower() else NU_0
               for m in materials}
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # ESIM solver
    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None, n_nodes=200,
        geometry='cylinder')

    # Initial Z_s
    sol0 = esim_solver.solve(5.0)
    Z_s = sol0['Z']
    print(f"  Initial Z_s = {Z_s:.4e}, delta = {sol0['delta']*1e3:.3f} mm")

    # Sample points for H_t evaluation (just outside workpiece boundary)
    eps = maxh * 0.2
    pts_lat = [(R_wp + eps, z) for z in
               np.linspace(-H_wp / 2 + 0.001, H_wp / 2 - 0.001, 20)]
    pts_cap = [(r, H_wp / 2 + eps) for r in np.linspace(0.002, R_wp - 0.001, 8)]
    pts_cap += [(r, -H_wp / 2 - eps) for r in np.linspace(0.002, R_wp - 0.001, 8)]
    sample_pts = pts_lat + pts_cap

    # --- Karl iteration ---
    print(f"  {'Iter':>4s} {'|Z_s|':>12s} {'P [W]':>12s} {'H_t_avg':>10s} {'dZ/Z':>10s}")

    t0_solve = time.perf_counter()
    for iteration in range(max_iter):
        # Robin BC: -(jw/Zs) * int u*v/r ds
        robin = -1j * omega / Z_s

        fesc = Periodic(H1(mesh, order=order, complex=True,
                           dirichlet="axis|axis_ext", dirichlet_bbnd="GND"))
        uc, vc = fesc.TnT()
        bf = BilinearForm(fesc)
        bf += nu_cf / r_cf * grad(uc) * grad(vc) * dx
        bf += robin * uc * vc / r_cf * ds("workpiece_bnd")
        lfc = LinearForm(fesc)
        lfc += J0 * vc * dx("coil")
        bf.Assemble()
        lfc.Assemble()
        gfu = GridFunction(fesc)
        gfu.vec.data = bf.mat.Inverse(fesc.FreeDofs(), inverse="pardiso") * lfc.vec

        # Power (Poynting): P = pi * w^2 * Re(Zs)/|Zs|^2 * int|u|^2/r ds
        int_u2r = Integrate(gfu * Conj(gfu) / r_cf, mesh,
                            definedon=mesh.Boundaries("workpiece_bnd")).real
        P = np.pi * omega**2 * Z_s.real / abs(Z_s)**2 * int_u2r

        # H_t from SIBC relation on the boundary (NOT offset sampling)
        # SIBC: H_t = (-jw/Z_s) * A_t, so |H_t| = w*|A_t|/|Z_s| = w*|u|/(r*|Z_s|)
        # 3D surface element dA = 2*pi*r*ds, so
        #   int |A_phi|^2 dA = 2*pi * int |u|^2/r ds = 2*pi * int_u2r
        #   <|A_phi|^2>_area = 2*pi * int_u2r / A_wp
        #   H_t_rms^2 = w^2/|Z_s|^2 * 2*pi * int_u2r / A_wp
        wp_area = 2 * np.pi * (R_wp * H_wp + 2 * 0.5 * R_wp**2)  # lateral + caps
        H_t_avg = omega / abs(Z_s) * np.sqrt(
            max(2 * np.pi * int_u2r / wp_area, 0)) if int_u2r > 0 else 1e-3

        # Update Z_s (under-relaxation)
        Z_s_old = Z_s
        sol_new = esim_solver.solve(max(float(H_t_avg), 1e-3), relaxation=0.3)
        Z_s = 0.5 * Z_s_old + 0.5 * sol_new['Z']

        dZ = abs(Z_s - Z_s_old) / abs(Z_s_old)
        print(f"  {iteration:4d} {abs(Z_s):12.4e} {P:12.4e} {H_t_avg:10.2f} {dZ:10.4e}")

        if dZ < tol and iteration > 0:
            break

    t_solve = time.perf_counter() - t0_solve

    # Inductance: volume magnetic energy + skin-layer reactive energy.
    # SIBC excludes the conductor interior from the mesh, so the skin
    # layer's stored magnetic energy is missing from the volume integral.
    # Add it back via the reactive part of Z_s (Im(Z_s) = Y > 0 for metals):
    #   <W_m_skin> = Im(Z_s)*|H_t|^2/(4*omega)  per 3D surface area
    # With |H_t|^2 = omega^2*|u|^2/(r^2*|Z_s|^2) and dA_3D = 2*pi*r*ds:
    #   <W_m_skin> = pi*omega*Y/(2|Z_s|^2) * int |u|^2/r ds
    # L = 4*<W_m_total>/|I|^2 = 2*W_volume_peak/|I|^2 + 2*pi*omega*Y/|Z_s|^2 * B
    # Verified: matches 2D full-res within ~1% across R/delta = 12..85
    # (without this term, L is systematically low by 5-10%).
    W_vol = 2 * np.pi * Integrate(
        0.5 * nu_cf / r_cf * grad(gfu) * Conj(grad(gfu)), mesh).real
    L_vol = 2 * W_vol / I_total**2
    L_skin = (2 * np.pi * omega * Z_s.imag / abs(Z_s)**2 * int_u2r
              / I_total**2)
    L = L_vol + L_skin

    print(f"  P (Poynting Robin)  = {P:.6e} W")
    print(f"  L (volume only)     = {L_vol*1e9:.3f} nH")
    print(f"  L (+skin surface)   = {L*1e9:.3f} nH  (skin adds {L_skin*1e9:+.3f} nH)")

    return {
        'mode': 'FEM-SIBC (Karl)',
        'P_total': float(P),
        'L': float(L),
        'Z_s': Z_s,
        'H_t_avg': float(H_t_avg),
        'iterations': iteration + 1,
        'ndof': fesc.ndof,
        'ne': mesh.ne,
        't_mesh': t_mesh,
        't_solve': t_solve,
    }


# ============================================================
# Main
# ============================================================
def main(mode='both', material='copper'):
    R_coil = 0.030
    a_coil = 0.003
    R_wp = 0.010
    H_wp = 0.020
    I_total = 1.0

    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]

    if material == 'steel':
        sigma = 2e6
        frequency = 7000
        mu_r = None
        bh_curve = STEEL_BH
        mu_r_display = "~1600 (BH curve)"
    else:
        sigma = 5.8e7
        frequency = 1000
        mu_r = 1.0
        bh_curve = None
        mu_r_display = "1.0"

    rho = 1.0 / sigma
    mu_eff = MU_0 * (mu_r if mu_r else 1.0)
    delta = math.sqrt(2 * rho / (2 * np.pi * frequency * mu_eff))

    print("=" * 65)
    print(f"FEM-full vs FEM-ESIM (2D axisymmetric + Kelvin) [{material}]")
    print("=" * 65)
    print(f"Coil:      R={R_coil*1e3:.0f} mm, a={a_coil*1e3:.1f} mm")
    print(f"Workpiece: R={R_wp*1e3:.0f} mm, H={H_wp*1e3:.0f} mm")
    print(f"Material:  {material} (sigma={sigma:.1e} S/m, mu_r={mu_r_display})")
    print(f"Frequency: {frequency} Hz")
    print(f"Skin depth: {delta*1e3:.3f} mm (xi = R_wp/delta = {R_wp/delta:.1f})")
    print()

    results = {}

    maxh_full = max(delta / 3, 0.0002)  # resolve skin layer for FEM-full

    if mode in ('both', 'full'):
        print("[FEM-full] Workpiece volume meshed (direct eddy current)")
        print("-" * 65)
        r_full = solve_fem_full(
            R_coil, a_coil, R_wp, H_wp, sigma, frequency,
            I_total=I_total, maxh=maxh_full, order=3,
            a_kelvin=0.15, z_offset=0.4)
        results['full'] = r_full
        print(f"  P = {r_full['P_total']:.6e} W")
        print(f"  L = {r_full['L']*1e9:.2f} nH")
        print(f"  Time: mesh={r_full['t_mesh']:.1f}s, "
              f"asm={r_full['t_asm']:.1f}s, solve={r_full['t_solve']:.1f}s")
        print()

    if mode in ('both', 'esim', 'oneway'):
        print("[FEM-ESIM one-way] Static FEM + ESIM per panel (BEM-equivalent)")
        print("-" * 65)
        r_oneway = solve_fem_esim_oneway(
            R_coil, a_coil, R_wp, H_wp, sigma, frequency,
            bh_curve=bh_curve, mu_r=mu_r if mu_r else 1.0,
            I_total=I_total, maxh=0.003, order=3,
            a_kelvin=0.15, z_offset=0.4)
        results['oneway'] = r_oneway
        print(f"  P = {r_oneway['P_total']:.6e} W")
        print(f"  L = {r_oneway['L']*1e9:.2f} nH")
        print()

    # Comparison
    print("=" * 65)
    print("Comparison")
    print("=" * 65)
    header = f"  {'':>20s}"
    if 'full' in results:
        header += f"  {'FEM-full':>12s}"
    if 'oneway' in results:
        header += f"  {'one-way':>12s}"
    print(header)

    if 'full' in results:
        P_f = results['full']['P_total']
        L_f = results['full']['L']
    if 'oneway' in results:
        P_o = results['oneway']['P_total']
        L_o = results['oneway']['L']

    row_P = f"  {'P [W]':>20s}"
    row_L = f"  {'L [nH]':>20s}"
    if 'full' in results:
        row_P += f"  {P_f:12.4e}"
        row_L += f"  {L_f*1e9:12.4f}"
    if 'oneway' in results:
        row_P += f"  {P_o:12.4e}"
        row_L += f"  {L_o*1e9:12.4f}"
    if 'full' in results and 'oneway' in results:
        row_P += f"  {(P_o-P_f)/max(abs(P_f),1e-30)*100:+8.1f}%"
        row_L += f"  {(L_o-L_f)/max(abs(L_f),1e-30)*100:+8.1f}%"
    print(row_P)
    print(row_L)
    print()
    print(f"  delta = {delta*1e3:.3f} mm, xi = {R_wp/delta:.1f}")
    print(f"  FEM-ESIM one-way: static FEM (no hole) + cylinder ESIM")

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="FEM A-formulation + Kelvin + ESIM prototype")
    parser.add_argument("--mode", default="both",
                        choices=["both", "full", "esim", "oneway"])
    parser.add_argument("--material", default="copper",
                        choices=["copper", "steel"])
    args = parser.parse_args()
    main(mode=args.mode, material=args.material)
