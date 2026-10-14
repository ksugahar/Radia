"""Clean test for scattered-field Robin RHS bug investigation.

HISTORICAL (2026-04-24): the production scattered path
(``calc_fem_kelvin.solve_fem_biot_savart`` + ``--source-mode scattered``)
was REMOVED from the main codebase on 2026-04-24 — the investigation
this script drove concluded that the Kelvin-free clean test cannot
reproduce the alpha-dependence observed in the full PEEC+Kelvin case,
so the "missing (n x B_s) surface term" fix could not be validated here.
This script is retained as a research artifact only; the Biot-Savart
helpers it imports (kelvin_source.biot_savart_A_cf / _B_cf) remain
available for other research uses (line-integral back-reaction etc.).

Isolates the physics from PEEC / Kelvin / multi-FD / proximity effects:

  - Single straight filament along z, I = 1A real (no PEEC complex phases)
  - Flat small Cu slab at y = d (wp "hole" in air mesh, sibc sideset)
  - Dirichlet outer air box (no Kelvin transformation)
  - Frequency 7 kHz (Cu skin << wp thickness, SIBC valid)

Analytical H_t at the center of the wp surface is exactly  I / (2*pi*d)
modulo a small finite-length correction (L=2m -> <0.2% shift at d=5cm).

For a small wp face a x a with a << d, H_t_rms is also ~= I/(2*pi*d):

    H_t_rms**2 = (I/2pi)**2 * [1/(2*(a^2/4+d^2)) + arctan(a/2d)/(a*d)]

Sweeps the "missing surface term" factor alpha in the scattered Robin RHS:

    f_lf = -robin * A_s_t * v_t * ds(sibc)                  # current bug
         + alpha * NU_0 * (n x B_s) . v_t * ds(sibc)        # candidate fix

If alpha=0.78 reproduces the analytical H_t, the 14% overshoot observed
in the full PEEC+Kelvin problem is caused by PEEC phase / Kelvin coupling
and NOT by the pure scattered formulation.

If alpha=1.0 reproduces the analytical H_t, the scattered formulation
derivation is correct and the 14% overshoot originates in the Karl
iteration / sum|I|=1.21 issue.
"""

import os
import sys
import math
import json
import argparse
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'src', 'radia'))
sys.path.insert(0, os.path.join(REPO, 'src', 'radia', 'panels'))


# ---------------------------------------------------------------------------
# Analytical H_t_rms over a square face a x a at distance d from an infinite
# straight wire (filament along z, face normal in +/-y):
#   |H_t|(x, z) = I d / (2 pi (x^2 + d^2))
#   H_t_rms^2 = (1/a^2) * int_{-a/2}^{a/2} int_{-a/2}^{a/2} |H_t|^2 dx dz
# Closed form (see header).
# ---------------------------------------------------------------------------
def analytical_H_t_rms(I, d, A):
    """H_t_rms over bottom face (x,z) in [-A, A]^2 at y = -d,
    for infinite wire along z-axis carrying current I.

    Derivation: |H_t|(x, z) = I d / (2 pi (x^2 + d^2))  (pure x-component)
    H_t_rms^2 = (I^2 / 8 pi^2) * [1/(A^2 + d^2) + arctan(A/d)/(A*d)]
    """
    term = 1.0 / (A * A + d * d) + math.atan(A / d) / (A * d)
    return abs(I) / (2.0 * math.pi * math.sqrt(2.0)) * math.sqrt(term)


def analytical_H_t_center(I, d, L_wire):
    """Finite-length wire: H = I/(4 pi d) * (sin a2 - sin a1).
    For wire from z=-L/2 to +L/2 at x=0, observation at (0, d, 0):
      a1 = -atan(L/2d), a2 = +atan(L/2d), sin a2 - sin a1 = 2 sin atan(L/2d).
    """
    s = math.sin(math.atan(L_wire / (2.0 * d)))
    return abs(I) / (2.0 * math.pi * d) * s


# ---------------------------------------------------------------------------
# Geometry: air box (Dirichlet "GND") minus flat Cu slab "hole" with "sibc"
# boundary. Uses netgen.occ directly. No Cubit required.
# ---------------------------------------------------------------------------
def build_vol(vol_file, A=0.05, d=0.05, h_gap=0.01,
              maxh_air=0.01, maxh_sibc=0.005):
    """Air box placed BELOW the wire so the filament never enters the domain.

    Wire runs along z at (x=0, y=0), z in [-L/2, L/2].
    Domain: (x, y, z) in [-A, A] x [-d, -h_gap] x [-A, A]
      - Bottom face (y = -d)             -> "sibc"
      - Other 5 faces                    -> "GND" (Dirichlet)

    h_gap > 0 is essential: it keeps the filament axis strictly outside the
    FEM domain, otherwise biot_savart_A_cf logs a zero denominator on the
    wire axis (at mesh vertices that happen to land on the z-piercing faces).

    sibc face size: 2A x 2A,  distance wire -> sibc: d.
    """
    from netgen.occ import Box, Pnt, OCCGeometry

    air = Box(Pnt(-A, -d, -A), Pnt(A, -h_gap, A))

    for f in air.faces:
        c = f.center
        if abs(c.y - (-d)) < 1e-8:
            f.name = "sibc"
            f.maxh = maxh_sibc
        else:
            f.name = "GND"
    air.mat("air")

    geo = OCCGeometry(air)
    ng_mesh = geo.GenerateMesh(maxh=maxh_air)
    ng_mesh.Save(vol_file)

    print(f"[mesh] {vol_file}  ne={ng_mesh.ne}")
    return vol_file


# ---------------------------------------------------------------------------
# Scattered FEM-SIBC solve with injectable RHS variant.
# Mirrors solve_fem_biot_savart in calc_fem_kelvin.py but (a) single filament
# instead of PEEC, (b) no Kelvin, (c) configurable rhs_factor for the
# (n x B_s) surface term.
# ---------------------------------------------------------------------------
def solve(vol_file, L_wire, I, d, A_half, frequency,
          rhs_factor, fes_order=1, solver="pardiso", mat_name="copper"):
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
                          Integrate, curl, dx, ds, CF, BND, VOL, TaskManager,
                          specialcf)
    global MU_0, NU_0
    from calc_common import MU_0 as _MU_0, NU_0 as _NU_0, EMMaterial
    MU_0 = _MU_0
    NU_0 = _NU_0
    from kelvin_source import biot_savart_A_cf, biot_savart_B_cf

    mesh = Mesh(vol_file)
    materials = mesh.GetMaterials()
    boundaries = set(mesh.GetBoundaries())
    assert 'sibc' in boundaries, "mesh missing 'sibc' boundary"
    dirichlet = 'GND' if 'GND' in boundaries else ''

    # --- Single straight filament: z in [-L_wire/2, +L_wire/2] at (0, 0) ---
    p1 = np.array([0.0, 0.0, -L_wire / 2.0])
    p2 = np.array([0.0, 0.0, +L_wire / 2.0])
    fil_paths = [[(p1, p2)]]
    fil_currents = np.array([complex(I)])

    # --- Biot-Savart A_s and B_s as symbolic CFs ---
    A_s_cf = biot_savart_A_cf(fil_paths, fil_currents)
    B_s_cf = biot_savart_B_cf(fil_paths, fil_currents)

    # --- SIBC ---
    mat = EMMaterial.from_name(mat_name)
    t_wp_for_zs = 0.005  # slab thickness; skin delta << 5mm at 7 kHz Cu
    Z_s = mat.dowell_Zs(frequency, t_wp_for_zs)
    omega = 2.0 * math.pi * frequency
    robin = 1j * omega / Z_s

    # --- FES (no Kelvin -> nu_cf = NU_0 everywhere) ---
    fes = HCurl(mesh, order=fes_order, nograds=True, complex=True,
                dirichlet=dirichlet)
    u, v_ = fes.TnT()

    a_bf = BilinearForm(fes)
    a_bf += NU_0 * curl(u) * curl(v_) * dx(bonus_intorder=4)
    a_bf += 1e-6 * NU_0 * u * v_ * dx  # gauge regularization
    a_bf += robin * u.Trace() * v_.Trace() * ds('sibc')

    n_bnd = specialcf.normal(3)

    f_lf = LinearForm(fes)
    # Current (buggy) RHS: -robin * A_s_t . v_t  on sibc
    f_lf += -robin * A_s_cf * v_.Trace() * ds('sibc', bonus_intorder=4)
    # Candidate missing surface term: +alpha * NU_0 * (n x B_s) . v_t
    if rhs_factor != 0.0:
        nxB = CF((n_bnd[1] * B_s_cf[2] - n_bnd[2] * B_s_cf[1],
                  n_bnd[2] * B_s_cf[0] - n_bnd[0] * B_s_cf[2],
                  n_bnd[0] * B_s_cf[1] - n_bnd[1] * B_s_cf[0]))
        f_lf += rhs_factor * NU_0 * nxB * v_.Trace() * ds('sibc',
                                                          bonus_intorder=4)

    t0 = time.perf_counter()
    a_bf.Assemble()
    f_lf.Assemble()
    t_asm = time.perf_counter() - t0

    gfu = GridFunction(fes)
    t0 = time.perf_counter()
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse=solver) * f_lf.vec
    t_solve = time.perf_counter() - t0

    # --- Measurements ---
    def _surface_rms(vec_cf, has_complex=True):
        """sqrt(int |vec_t|^2 dS / area) for tangential part on sibc."""
        sq = sum(vec_cf[i].real * vec_cf[i].real
                 + vec_cf[i].imag * vec_cf[i].imag for i in range(3))
        dn_re = sum(vec_cf[i].real * n_bnd[i] for i in range(3))
        dn_im = sum(vec_cf[i].imag * n_bnd[i] for i in range(3))
        tsq = sq - (dn_re * dn_re + dn_im * dn_im)
        wp_region = mesh.Boundaries('sibc')
        area = Integrate(CF(1), mesh, BND, definedon=wp_region).real
        val = Integrate(tsq, mesh, BND, definedon=wp_region, order=8).real
        return math.sqrt(max(val, 0) / max(area, 1e-30)), area

    from ngsolve import curl as ng_curl

    # A_total = A_s (CF) + A_r (gfu)
    A_total = CF(tuple(A_s_cf[i] + gfu[i] for i in range(3)))
    # Also build A_s only and A_r only (for diagnostics)
    A_s_only = A_s_cf
    A_r_only = CF(tuple(gfu[i] for i in range(3)))
    # curl(A_total) gives B_total via curl of GridFunction + CF
    B_total = CF(tuple((B_s_cf[i] + ng_curl(gfu)[i]) for i in range(3)))

    At_rms_total, A_wp = _surface_rms(A_total)
    At_rms_s, _ = _surface_rms(A_s_only)
    At_rms_r, _ = _surface_rms(A_r_only)
    Bt_rms_total, _ = _surface_rms(B_total)

    # Robin relation: |H_t| = omega * |A_t| / |Z_s|  (valid if BC satisfied)
    H_t_rms_robin = abs(1j * omega / Z_s) * At_rms_total
    # Physical: |H_t| from B_total / MU_0
    H_t_rms_phys = Bt_rms_total / MU_0
    # Physical from A_s only (free-space)
    Bt_rms_s, _ = _surface_rms(B_s_cf)
    H_t_rms_free = Bt_rms_s / MU_0

    return {
        "rhs_factor": rhs_factor,
        "H_t_rms_robin": H_t_rms_robin,
        "H_t_rms_phys": H_t_rms_phys,
        "H_t_rms_free": H_t_rms_free,
        "At_rms_total": At_rms_total,
        "At_rms_s": At_rms_s,
        "At_rms_r": At_rms_r,
        "A_wp": A_wp,
        "ndof": fes.ndof,
        "t_asm": round(t_asm, 2),
        "t_solve": round(t_solve, 2),
        "Z_s_real": Z_s.real,
        "Z_s_imag": Z_s.imag,
        "omega_over_Zs": float(abs(1j * omega / Z_s)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", default=os.path.join(HERE, "clean_filament.vol"))
    parser.add_argument("--rebuild-mesh", action="store_true")
    parser.add_argument("--L-wire", type=float, default=2.0,
                        help="Filament length along z [m] (longer -> closer to infinite)")
    parser.add_argument("--I", type=float, default=1.0)
    parser.add_argument("--d", type=float, default=0.05,
                        help="Distance from filament to sibc face [m]")
    parser.add_argument("--A-half", type=float, default=0.05,
                        help="Air box half-width in x and z [m]")
    parser.add_argument("--h-gap", type=float, default=0.01,
                        help="Gap between box top and wire [m]; wire is OUTSIDE the domain")
    parser.add_argument("--maxh-air", type=float, default=0.01)
    parser.add_argument("--maxh-sibc", type=float, default=0.005)
    parser.add_argument("--frequency", type=float, default=7000)
    parser.add_argument("--fes-order", type=int, default=1)
    parser.add_argument("--alphas", default="0.0,0.5,0.78,1.0,1.2",
                        help="Comma-separated rhs_factor values to sweep")
    parser.add_argument("--output", default=os.path.join(HERE, "results.json"))
    args = parser.parse_args()

    if args.rebuild_mesh or not os.path.exists(args.vol):
        build_vol(args.vol, A=args.A_half, d=args.d, h_gap=args.h_gap,
                  maxh_air=args.maxh_air, maxh_sibc=args.maxh_sibc)

    # Analytical targets
    H_t_rms_ana = analytical_H_t_rms(args.I, args.d, args.A_half)
    H_t_center_ana_inf = abs(args.I) / (2.0 * math.pi * args.d)
    H_t_center_ana_finite = analytical_H_t_center(args.I, args.d, args.L_wire)

    # Target for a Cu half-space SIBC (close to PEC): H_t ~= 2 * H_free
    H_t_rms_pec = 2.0 * H_t_rms_ana

    alphas = [float(s) for s in args.alphas.split(",")]
    results = []
    for alpha in alphas:
        print(f"\n=== alpha = {alpha} ===")
        r = solve(args.vol, args.L_wire, args.I, args.d, args.A_half,
                   args.frequency, rhs_factor=alpha, fes_order=args.fes_order)
        r["H_t_rms_free_ana"] = H_t_rms_ana
        r["H_t_rms_pec_ana"] = H_t_rms_pec
        r["ratio_robin_to_pec"] = r["H_t_rms_robin"] / H_t_rms_pec
        r["ratio_phys_to_pec"] = r["H_t_rms_phys"] / H_t_rms_pec
        r["ratio_phys_to_free"] = r["H_t_rms_phys"] / r["H_t_rms_free"]
        results.append(r)
        print(f"   H_t_rms [Robin] = {r['H_t_rms_robin']:.4f}  "
              f"[curl(A_tot)] = {r['H_t_rms_phys']:.4f}  "
              f"[free from B_s] = {r['H_t_rms_free']:.4f}")
        print(f"   A_t_rms  total = {r['At_rms_total']:.3e}  "
              f"s = {r['At_rms_s']:.3e}  r = {r['At_rms_r']:.3e}")
        print(f"   ndof = {r['ndof']}, |omega/Z_s| = {r['omega_over_Zs']:.3e}")

    out = {
        "args": vars(args),
        "analytical": {
            "H_t_rms_free": H_t_rms_ana,
            "H_t_rms_pec": H_t_rms_pec,
            "H_t_center_infinite": H_t_center_ana_inf,
            "H_t_center_finite": H_t_center_ana_finite,
        },
        "sweep": results,
    }
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\n=== Summary ===")
    print(f"analytical H_t_rms free-space      = {H_t_rms_ana:.4f} A/m")
    print(f"analytical H_t_rms PEC (2 x free)  = {H_t_rms_pec:.4f} A/m")
    print(f"{'alpha':>6}  {'H_Robin':>10}  {'H_phys':>10}  "
          f"{'H_free':>10}  {'ph/pec':>7}  {'ph/free':>8}")
    for r in results:
        print(f"{r['rhs_factor']:>6.2f}  "
              f"{r['H_t_rms_robin']:>10.4f}  "
              f"{r['H_t_rms_phys']:>10.4f}  "
              f"{r['H_t_rms_free']:>10.4f}  "
              f"{r['ratio_phys_to_pec']:>7.3f}  "
              f"{r['ratio_phys_to_free']:>8.3f}")
    print(f"\nResults: {args.output}")


if __name__ == "__main__":
    main()
