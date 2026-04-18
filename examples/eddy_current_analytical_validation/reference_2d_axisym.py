"""2D axisymmetric FEM reference solution for coil + cylinder workpiece.

Resolves eddy currents DIRECTLY inside the workpiece (no SIBC).
Mesh is refined to delta/5 inside the conductor for accurate skin
depth resolution.

Geometry (matches ih_bem_sample.jou):
  Coil:      circle cross-section at (r=30mm, z=0), radius 3mm
  Workpiece: solid cylinder r=0..25mm, z=-12.5..12.5mm
  Air:       half-circle (Kelvin) or rectangle (no Kelvin)

A-formulation in axisymmetric (r,z) with phi = r*A_phi substitution:
  Bilinear: nu/r * grad(phi).grad(psi) + jw*sigma*phi*psi/r
  B_r = -(1/r)*dphi/dz,  B_z = (1/r)*dphi/dr

Output: L, P, H_t_rms on workpiece surface.

Usage:
    python reference_2d_axisym.py --material copper --freq 7000
    python reference_2d_axisym.py --sweep
    python reference_2d_axisym.py --sweep --kelvin
"""

import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
from em_material import EMMaterial, MU_0


def solve_2d_axisym(r_coil=0.030, a_coil=0.003, r_wp=0.025, h_wp=0.025,
                    mat=None, frequency=7000, I_total=1.0, maxh_air=0.010,
                    kelvin=False, order=2,
                    skin_mode="uniform", maxh_wp=None,
                    maxh_wp_interior=None, skin_ratio=5.0,
                    grading=0.5,
                    sigma_coil=0.0, mu_r_coil=1.0, maxh_coil_override=None):
    """2D axisymmetric FEM with full eddy current resolution.

    Args:
        kelvin: If True, use z-offset Kelvin transformation for open
            boundary.  Eliminates the conducting-shell artifact from
            Dirichlet BC on a finite air domain.
        order: FE polynomial order (default 2).
        sigma_coil: coil conductivity [S/m]. 0 (default) = uniform-J
            DC mode (existing behavior, matches 2D SIBC). Set to
            5.8e7 (Cu) to capture coil skin/proximity self-consistently
            -- this matches the physics of a PEEC nwinc>=3 filament
            bundle and is the apples-to-apples reference for PEEC+FEM
            and PEEC+BEM solvers at AC frequencies.
        mu_r_coil: coil relative permeability (for magnetic coils).
        maxh_coil_override: mesh size in coil. When sigma_coil > 0 the
            coil needs h <= delta_coil/3 to resolve the skin layer;
            default a_coil/2 may be too coarse.

    Returns dict with L, P, H_t_rms, diagnostics.
    """
    if mat is None:
        mat = EMMaterial.from_name("steel")

    from ngsolve import (Mesh, H1, Periodic, BilinearForm, LinearForm,
                         GridFunction, Integrate, Conj, grad, sqrt as ngsqrt,
                         dx, ds, CF, CoefficientFunction, BND, VOL,
                         TaskManager, IfPos)
    from ngsolve import x, y
    from netgen.occ import (WorkPlane, OCCGeometry, Glue, Axes, Pnt,
                            Vertex, MoveTo, X, Y, Z)
    from netgen.meshing import IdentificationType

    omega = 2 * math.pi * frequency
    sigma = mat.sigma
    mu_r_val = mat.mu_r
    nu0 = 1.0 / MU_0
    nu_wp = 1.0 / (mu_r_val * MU_0)
    delta = mat.skin_depth(frequency)

    if maxh_wp is None:
        maxh_wp = delta / skin_ratio
    if maxh_wp_interior is None:
        maxh_wp_interior = max(delta * 3, maxh_wp)
    if maxh_coil_override is not None:
        maxh_coil = maxh_coil_override
    elif sigma_coil > 0:
        # Resolve skin depth inside conducting coil
        delta_coil = math.sqrt(2.0 / (omega * mu_r_coil * MU_0 * sigma_coil))
        maxh_coil = min(a_coil / 2, delta_coil / 3)
    else:
        maxh_coil = a_coil / 2
    nu_coil = 1.0 / (mu_r_coil * MU_0)

    t0 = time.perf_counter()

    r_safe = IfPos(x - 1e-10, x, 1e-10)

    if kelvin:
        # --- Kelvin geometry: half-circle interior + exterior ---
        a_kelvin = max(r_coil * 3, h_wp + r_coil * 2, 0.10)
        z_offset = a_kelvin * 2.5

        # Interior half-circle centered at (0, 0)
        wp1 = WorkPlane()
        inner_full = wp1.Circle(a_kelvin).Face()
        inner_full.name = "air"

        # Coil
        coil_full = MoveTo(r_coil - a_coil, -a_coil).Rectangle(
            2 * a_coil, 2 * a_coil).Face()
        # Better: use circle for coil
        wp_coil = WorkPlane()
        wp_coil.MoveTo(r_coil + a_coil, 0)
        wp_coil.Direction(0, 1)  # heading UP so arc center = (r_coil, 0)
        wp_coil.Arc(a_coil, 180).Arc(a_coil, 180)
        coil_full = wp_coil.Face()
        coil_full.name = "coil"
        coil_full.maxh = maxh_coil

        # Workpiece
        wp_cond = WorkPlane()
        wp_cond.MoveTo(0, -h_wp / 2).LineTo(r_wp, -h_wp / 2).LineTo(
            r_wp, h_wp / 2).LineTo(0, h_wp / 2).LineTo(0, -h_wp / 2)
        cond_full = wp_cond.Face()
        cond_full.name = "workpiece"
        if skin_mode == "skin":
            cond_full.maxh = maxh_wp_interior
            for e in cond_full.edges:
                if e.center.x > 1e-6:
                    e.maxh = maxh_wp
        else:
            cond_full.maxh = maxh_wp

        # Exterior half-circle centered at (0, z_offset)
        wp2 = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
        outer_full = wp2.Circle(a_kelvin).Face()
        outer_full.name = "air_kelvin"

        # GND point at center of exterior (represents infinity)
        gnd_point = Vertex(Pnt(0, z_offset, 0))
        gnd_point.name = "GND"

        # Cut to half-circles (x >= 0)
        cut_left = MoveTo(-a_kelvin - 0.1, -a_kelvin - 0.1).Rectangle(
            a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
        cut_left_ext = MoveTo(-a_kelvin - 0.1,
                              z_offset - a_kelvin - 0.1).Rectangle(
            a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()

        inner_half = inner_full - cut_left
        coil_half = coil_full - cut_left
        cond_half = cond_full - cut_left
        outer_half = outer_full - cut_left_ext

        inner_air = inner_half - coil_half - cond_half
        inner_air.name = "air"

        # Name edges: axis, kelvin arcs
        for edge in inner_air.edges:
            cx = edge.center.x
            try:
                v0, v1 = edge.vertices
                d0 = math.sqrt(v0.p.x ** 2 + v0.p.y ** 2)
                d1 = math.sqrt(v1.p.x ** 2 + v1.p.y ** 2)
                is_arc = (abs(d0 - a_kelvin) < 0.01
                          and abs(d1 - a_kelvin) < 0.01 and cx > 0.01)
            except Exception:
                is_arc = False
            if cx < 0.01:
                edge.name = "axis"
            elif is_arc:
                edge.name = "kelvin_int"

        for edge in coil_half.edges:
            if edge.center.x < 0.01:
                edge.name = "axis"

        for edge in cond_half.edges:
            if edge.center.x < 0.01:
                edge.name = "axis"

        kelvin_outer_edges = []
        for edge in outer_half.edges:
            cx = edge.center.x
            try:
                v0, v1 = edge.vertices
                d0 = math.sqrt(v0.p.x ** 2 + (v0.p.y - z_offset) ** 2)
                d1 = math.sqrt(v1.p.x ** 2 + (v1.p.y - z_offset) ** 2)
                is_arc = (abs(d0 - a_kelvin) < 0.01
                          and abs(d1 - a_kelvin) < 0.01 and cx > 0.01)
            except Exception:
                is_arc = False
            if cx < 0.01:
                edge.name = "axis_ext"
            elif is_arc:
                edge.name = "kelvin_ext"
                kelvin_outer_edges.append(edge)

        # Periodic identification: match arcs by y-sign
        kelvin_inner_edges = [e for e in inner_air.edges
                              if getattr(e, 'name', '') == "kelvin_int"]
        for int_edge in kelvin_inner_edges:
            int_y = int_edge.center.y
            for ext_edge in kelvin_outer_edges:
                ext_y = ext_edge.center.y - z_offset
                if (int_y > 0 and ext_y > 0) or (int_y < 0 and ext_y < 0):
                    int_edge.Identify(ext_edge, "kelvin",
                                      IdentificationType.PERIODIC)
                    break

        shape = Glue([inner_air, coil_half, cond_half, outer_half, gnd_point])
        geo = OCCGeometry(shape, dim=2)
        with TaskManager():
            ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=grading)
        mesh = Mesh(ngmesh)
        mesh.Curve(order)

        # FES with Periodic BC
        fes_pre = H1(mesh, order=order, complex=True,
                     dirichlet="axis|axis_ext", dirichlet_bbnd="GND")
        fes = Periodic(fes_pre)

        # Reluctivity: Kelvin factor (rho'/a)^2 in exterior domain
        y_local = y - z_offset
        rho_prime = ngsqrt(x ** 2 + y_local ** 2)
        rho_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)
        kelvin_factor = (rho_safe / a_kelvin) ** 2

        mat_list = list(mesh.GetMaterials())
        nu_vals = []
        sigma_vals = []
        for m in mat_list:
            if "kelvin" in m:
                nu_vals.append(nu0 * kelvin_factor)
                sigma_vals.append(0)
            elif m == "workpiece":
                nu_vals.append(nu_wp)
                sigma_vals.append(sigma)
            elif m == "coil":
                nu_vals.append(nu_coil)
                sigma_vals.append(sigma_coil)
            else:  # air
                nu_vals.append(nu0)
                sigma_vals.append(0)
        nu_cf = CoefficientFunction(nu_vals)
        sigma_cf = CoefficientFunction(sigma_vals)

    else:
        # --- Rectangular air domain (original) ---
        r_air = max(r_coil * 5, 0.15)
        z_air = max(r_air, h_wp + r_coil * 3)

        wp_air = WorkPlane()
        wp_air.MoveTo(0, -z_air).LineTo(r_air, -z_air).LineTo(
            r_air, z_air).LineTo(0, z_air).LineTo(0, -z_air)
        air_face = wp_air.Face()
        air_face.name = "air"
        for e in air_face.edges:
            e.name = "axis" if abs(e.center.x) < 1e-6 else "outer"

        wp_cond = WorkPlane()
        wp_cond.MoveTo(0, -h_wp / 2).LineTo(r_wp, -h_wp / 2).LineTo(
            r_wp, h_wp / 2).LineTo(0, h_wp / 2).LineTo(0, -h_wp / 2)
        cond_face = wp_cond.Face()
        cond_face.name = "workpiece"
        if skin_mode == "skin":
            cond_face.maxh = maxh_wp_interior
            for e in cond_face.edges:
                if e.center.x > 1e-6:
                    e.maxh = maxh_wp
        else:
            cond_face.maxh = maxh_wp
        for e in cond_face.edges:
            if abs(e.center.x) < 1e-6:
                e.name = "axis"

        wp_coil = WorkPlane()
        wp_coil.MoveTo(r_coil + a_coil, 0)
        wp_coil.Direction(0, 1)  # heading UP so arc center = (r_coil, 0)
        wp_coil.Arc(a_coil, 180).Arc(a_coil, 180)
        coil_face = wp_coil.Face()
        coil_face.name = "coil"
        coil_face.maxh = maxh_coil

        air_cut = air_face - cond_face - coil_face
        air_cut.name = "air"
        shape = Glue([air_cut, cond_face, coil_face])
        geo = OCCGeometry(shape, dim=2)
        with TaskManager():
            ngmesh = geo.GenerateMesh(maxh=maxh_air, grading=grading)
        mesh = Mesh(ngmesh)
        mesh.Curve(order)

        fes = H1(mesh, order=order, complex=True,
                 dirichlet="outer|axis")
        nu_cf = mesh.MaterialCF({"air": nu0, "coil": nu_coil,
                                  "workpiece": nu_wp}, default=nu0)
        sigma_cf = mesh.MaterialCF({"workpiece": sigma, "coil": sigma_coil},
                                     default=0)

    # ---- Common solve path ----
    t_mesh = time.perf_counter() - t0
    ne = mesh.ne
    ndof = fes.ndof
    u, v = fes.TnT()

    def Bop(phi):
        """B from phi = r*A_phi."""
        return CF((-grad(phi)[1] / r_safe, grad(phi)[0] / r_safe))

    # Bilinear form: nu/r * grad . grad  +  jw*sigma/r * phi*psi
    # (equivalent to nu*r*Bop.Bop + jw*sigma*u*v/r)
    a_bf = BilinearForm(fes)
    a_bf += nu_cf / r_safe * grad(u) * grad(v) * dx
    a_bf += 1j * omega * sigma_cf * u * v / r_safe * dx
    a_bf.Assemble()

    J0 = I_total / (math.pi * a_coil ** 2)
    f_lf = LinearForm(fes)
    f_lf += J0 * v * dx("coil")
    f_lf.Assemble()

    # Solve
    t0 = time.perf_counter()
    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(),
                                         inverse="pardiso") * f_lf.vec
    t_solve = time.perf_counter() - t0

    # Current normalization for sigma_coil > 0 (conductive coil):
    # The imposed J0 alone gives total I = J0 * pi a^2 = I_target, but
    # with sigma > 0 the coil eddy -jw sigma A_phi subtracts from J0,
    # so the actual net current through the coil cross-section differs
    # from I_target. Compute I_actual and rescale gfu linearly so the
    # re-scaled solution corresponds to net I = I_target.
    #   J_phi_actual = J0 - jw sigma A_phi = J0 - jw sigma * phi / r
    #   I_actual = integral over coil meridional = int(J_phi_actual dr dz)
    #     (no 2pi, no r-weight -- this is the net azimuthal current).
    if sigma_coil > 0:
        J_phi_cf = J0 - 1j * omega * sigma_coil * gfu / r_safe
        I_actual = complex(Integrate(
            J_phi_cf, mesh, definedon=mesh.Materials("coil")))
        scale = complex(I_total) / I_actual
        gfu.vec.data = complex(scale) * gfu.vec
    else:
        scale = 1.0 + 0j

    # L and P from impedance (no 1/r singularity at axis)
    # NOTE: ftu uses raw gfu; if we scaled gfu, also scale f_lf-weighted
    # inner product by the same factor (it's in gfu now already).
    ftu = sum(f_lf.vec[i] * gfu.vec[i].conjugate()
              for i in range(fes.ndof)).real
    L = 2 * math.pi * ftu / I_total ** 2

    ftu_im = sum(f_lf.vec[i].real * gfu.vec[i].imag
                 for i in range(fes.ndof))
    P_total = -math.pi * omega * ftu_im

    # H_t on workpiece cylindrical surface
    B_cf = Bop(gfu)
    n_sample = 50
    z_pts = np.linspace(-h_wp / 2 * 0.95, h_wp / 2 * 0.95, n_sample)
    H_t_sq_sum = 0.0
    n_valid = 0
    for zi in z_pts:
        try:
            mip = mesh(r_wp * 1.001, zi)
            Bz_val = B_cf[1](mip)
            Hz_val = Bz_val / MU_0
            H_t_sq_sum += abs(Hz_val) ** 2
            n_valid += 1
        except Exception:
            pass
    H_t_rms_side = math.sqrt(H_t_sq_sum / max(n_valid, 1))

    return {
        'L': float(L),
        'P_total': float(P_total),
        'H_t_rms': float(H_t_rms_side),
        'W_mag': 0.5 * float(L) * I_total ** 2,
        'ndof': ndof,
        'ne': ne,
        'delta': float(delta),
        'maxh_wp': float(maxh_wp),
        't_mesh': round(t_mesh, 2),
        't_solve': round(t_solve, 2),
        'frequency': frequency,
        'sigma': sigma,
        'mu_r': mu_r_val,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="2D axisymmetric FEM reference (full eddy current)")
    parser.add_argument("--material", default="copper",
                        choices=["copper", "steel", "aluminum"])
    parser.add_argument("--freq", type=float, default=7000)
    parser.add_argument("--sweep", action="store_true")
    parser.add_argument("--kelvin", action="store_true",
                        help="Use z-offset Kelvin for open boundary")
    parser.add_argument("--sigma-coil", type=float, default=0.0,
                        help="Coil conductivity [S/m]. 0 (default) = "
                             "uniform-J DC coil. 5.8e7 = Cu full-eddy "
                             "coil (apples-to-apples with PEEC nwinc>=3).")
    args = parser.parse_args()

    bc_label = "Kelvin" if args.kelvin else "Dirichlet"
    print("=" * 70)
    print(f"2D Axisymmetric FEM Reference (Full Eddy Current, {bc_label} BC)")
    print("Coil R=30mm a=3mm, WP R=25mm H=25mm")
    print("=" * 70)

    if args.sweep:
        print(f"\n{'Mat':>6s} {'f[Hz]':>8s} {'R/d':>6s} {'maxh_wp':>8s} "
              f"{'L[nH]':>8s} {'P[W]':>10s} {'H_t[A/m]':>10s} "
              f"{'ndof':>6s} {'t[s]':>5s}")
        print("-" * 80)
        for mname in ['copper', 'steel']:
            mat = EMMaterial.from_name(mname)
            for freq in [1000, 7000, 50000]:
                try:
                    r = solve_2d_axisym(mat=mat, frequency=freq,
                                        kelvin=args.kelvin)
                    delta = mat.skin_depth(freq)
                    rd = 0.025 / delta
                    print(f"{mname:>6s} {freq:>8.0f} {rd:>6.1f} "
                          f"{r['maxh_wp']*1e3:>7.3f}mm "
                          f"{r['L']*1e9:>8.2f} {r['P_total']:>10.2e} "
                          f"{r['H_t_rms']:>10.1f} "
                          f"{r['ndof']:>6d} {r['t_solve']:>5.1f}")
                except Exception as e:
                    print(f"{mname:>6s} {freq:>8.0f}  FAILED: {e}")
    else:
        mat = EMMaterial.from_name(args.material)
        print(f"\nMaterial: {mat.name} (sigma={mat.sigma:.2e}, "
              f"mu_r={mat.mu_r})")
        print(f"Frequency: {args.freq:.0f} Hz")
        delta = mat.skin_depth(args.freq)
        print(f"Skin depth: {delta*1e3:.3f} mm")
        print()

        r = solve_2d_axisym(mat=mat, frequency=args.freq,
                            kelvin=args.kelvin,
                            sigma_coil=args.sigma_coil)
        print(f"L       = {r['L']*1e9:.4f} nH")
        print(f"P_total = {r['P_total']:.6e} W")
        print(f"H_t_rms = {r['H_t_rms']:.4f} A/m (cylindrical surface)")
        print(f"ndof    = {r['ndof']}, ne = {r['ne']}")
        print(f"t_mesh  = {r['t_mesh']:.1f}s, t_solve = {r['t_solve']:.1f}s")
