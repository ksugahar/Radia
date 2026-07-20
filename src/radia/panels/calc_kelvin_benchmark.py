"""
Kelvin Benchmark: verify a Cubit-meshed Kelvin pipeline against a known
analytical solution.

Loads a `.vol` (with `magnetic`, `air`, `kelvin` materials and
`sphere`, `kelvin_int`, `kelvin_ext` sidesets), applies a uniform
external field along one of the principal axes, solves the
Omega-Reduced Omega + Kelvin scalar-potential equation and probes the
field at the probe point (default origin).  Compares to the analytical
interior field for a magnetic sphere in uniform field:

    H_inside = 3 / (mu_r + 2) * H0  along the field-axis direction

JSON output keys:
    Hi_origin       float   probed component along the field-axis
    Hi_analytical   float   3 / (mu_r + 2) * H0
    error_pct       float   100 * (Hi_origin - Hi_analytical) / Hi_analytical
    ne, ndof, slaved_dofs   FES diagnostics
    converged       bool
    n_iterations    int     (linear solve = 1)

Required sidesets / blocks in the .vol:
    blocks:    magnetic, air, kelvin
    sidesets:  sphere       (mag-air interface)
               kelvin_int   (air outer cap, paired with kelvin_ext)
               kelvin_ext   (kelvin outer cap)
    bbnodeset: GND          (Kelvin sphere centre, Dirichlet anchor)

Usage (called by the Electromagnet Simulink block / `EMDesignSpec` Kelvin Benchmark mode):
    python calc_kelvin_benchmark.py --vol model.vol --mu-r 100 \\
        --H0 1.0 --field-axis z --fes-order 2 --R-kelvin 0.20

Layer 4 contract (CLAUDE.md): NO `import cubit`, NO PySide6.
"""

import argparse
import json
import math
import os
import sys

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (
    setup_paths, progress, calc_main,
    detect_kelvin_offset, add_periodic_kelvin,
)


_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def _log(msg):
    progress("KELVIN_BENCH", msg)


def solve_kelvin_benchmark(vol_path, mu_r=100.0, H0=1.0,
                            field_axis="z", fes_order=2,
                            R_kelvin=0.20,
                            probe=(1e-10, 1e-10, 1e-10),
                            output_msh=None):
    """Run Omega-Reduced Omega + Kelvin and return field probe vs analytical.

    Args:
        vol_path: Netgen .vol with magnetic / air / kelvin blocks and
            sphere / kelvin_int / kelvin_ext sidesets.
        mu_r: Magnetic sphere relative permeability.
        H0: Uniform external field amplitude [A/m].
        field_axis: "x" | "y" | "z" -- direction of H0.
        fes_order: H1 polynomial order.  Use >= 2 for Kelvin (curved
            geometry) per memory/feedback_kelvin_high_order_required.md.
        R_kelvin: Kelvin sphere radius [m].  Must match the .vol mesh.
        probe: (x, y, z) probe point in the magnetic interior.
        output_msh: Optional path to write a GMSH .msh post-export.
    """
    if field_axis not in _AXIS_INDEX:
        raise ValueError(f"--field-axis must be x|y|z, got {field_axis!r}")
    field_idx = _AXIS_INDEX[field_axis]

    setup_paths()

    import ngsolve
    from ngsolve import (Mesh, H1, Periodic, BilinearForm, LinearForm,
                         GridFunction, CoefficientFunction, BND, ds, dx,
                         specialcf, x, y, z, grad, sqrt, TaskManager)

    # Pull the Kelvin reluctivity helper from radia.kelvin_source.
    sys.path.insert(0, os.path.join(_this_dir, ".."))
    from kelvin_source import (kelvin_mu_factor_3d_cf, build_material_cf)

    mu0 = 4.0 * math.pi * 1.0e-7

    # ---- Load mesh + Kelvin Periodic identification ----
    if not vol_path or not os.path.exists(vol_path):
        return {"error": f".vol not found: {vol_path}"}

    _log(f"loading {vol_path}")
    mesh = Mesh(vol_path)
    materials = mesh.GetMaterials()
    boundaries = mesh.GetBoundaries()
    _log(f"materials={list(materials)}, BNDs={sorted(set(boundaries))}")

    for required in ("magnetic", "air", "kelvin"):
        if required not in materials:
            return {"error": f"required material {required!r} missing; "
                             f"got {list(materials)}"}
    for required in ("sphere", "kelvin_int", "kelvin_ext"):
        if required not in boundaries:
            return {"error": f"required boundary {required!r} missing; "
                             f"got {sorted(set(boundaries))}"}

    offset = detect_kelvin_offset(mesh)
    _log(f"kelvin_offset = {offset}")
    if add_periodic_kelvin(mesh, offset):
        mesh = ngsolve.Mesh(mesh.ngmesh)

    if fes_order >= 2:
        mesh.Curve(fes_order)
        _log(f"mesh.Curve({fes_order}) applied")

    # ---- Build Periodic FES with Dirichlet BCs ----
    #
    # Note (2026-04-26): `kelvin_far` is intentionally NOT included in
    # the Dirichlet set when it is present on a 1/8 reduction mesh.
    # For H1 + Kelvin transform, the "infinity" anchor only needs ONE
    # Dirichlet point -- the `GND` vertex/nodeset placed at the Kelvin
    # centre (= image of physical infinity).  Constraining the whole
    # kelvin_far face overconstrains the system; it doesn't match the
    # physical sym_bn=0 BC that the offset axis has in the un-mapped
    # domain, and it removes ~250 DOFs that should be free.  This fix
    # alone does NOT make the 1/8 sample produce the analytical
    # answer (a separate corner-discretisation / sym_ht=0_z issue
    # gives Hz ~ 0 there), but the change is theoretically correct
    # for any future 1/8 path that gets the corner right and is
    # transparent for 1/2 and 1/4 (they have no kelvin_far face).
    bnd = list(mesh.GetBoundaries())
    bbb = list(mesh.GetBBBoundaries()) if hasattr(mesh, "GetBBBoundaries") else []
    dir_parts = []
    for label in ("GND", "sym_ht=0_x", "sym_ht=0_y", "sym_ht=0_z"):
        if label in bnd or label in bbb:
            dir_parts.append(label)
    if "GND" not in dir_parts:
        return {"error": "GND vertex (Kelvin centre) Dirichlet missing -- "
                         "add a `GND` nodeset / vertex name in the .jou."}
    dirichlet = "|".join(dir_parts)
    _log(f"dirichlet = {dirichlet!r}")

    fes_plain = H1(mesh, order=fes_order, dirichlet=dirichlet)
    fes = Periodic(fes_plain)
    free_total = sum(1 for d in fes.FreeDofs() if d)
    free_plain = sum(1 for d in fes_plain.FreeDofs() if d)
    slaved = free_plain - free_total
    _log(f"FES: ndof={fes.ndof}, FreeDofs {free_plain}->{free_total}, "
         f"slaved={slaved}")

    # ---- Material reluctivity (Mu) ----
    mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=tuple(offset),
                                              R=R_kelvin)
    Mu = build_material_cf(
        mesh, mu0, mu_kelvin_factor,
        outer_keyword="kelvin",
        overrides={"magnetic": mu_r * mu0},
    )

    # ---- Source (uniform H0 along field_axis) ----
    coords = (x, y, z)
    Omega_s = H0 * coords[field_idx]
    Bs_components = [0.0, 0.0, 0.0]
    Bs_components[field_idx] = mu0 * H0
    Bs = CoefficientFunction(tuple(Bs_components))

    # ---- Bilinear form: Mu * grad(u) * grad(v) over all materials ----
    u = fes.TrialFunction()
    v = fes.TestFunction()
    a = BilinearForm(fes, symmetric=True)
    a += Mu * grad(u) * grad(v) * dx

    # ---- Linear form: source lift in air + Neumann correction on sphere ----
    gfOmega = GridFunction(fes)
    gfOmega.vec[:] = 0.0
    gfOmega.Set(Omega_s, BND, mesh.Boundaries("sphere"))

    f = LinearForm(fes)
    f += Mu * grad(gfOmega) * grad(v) * dx("air")
    f.Assemble()

    # Remove Dirichlet components from f.vec (matches the archived
    # classic Omega_ReducedOmega.py pattern).
    import numpy as _np
    fcut = _np.array(f.vec.FV())[fes.FreeDofs()]
    _np.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

    # Cubit-meshed FaceDescriptor for "sphere" sideset has the
    # (domin=air, domout=mag) orientation; the Neumann correction
    # therefore uses -specialcf.normal (verified 2026-04-25 against
    # Cubit_1_4_p_convergence at p=2 -> +0.71% error).
    normal = -specialcf.normal(mesh.dim)
    f += (normal * Bs) * v * ds("sphere")
    f.Assemble()

    # ---- Solve ----
    _log("assembling + solving (pardiso)")
    with TaskManager():
        a.Assemble()
        gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                          inverse="pardiso") * f.vec

    # ---- Probe field at probe point (in magnetic region: H = grad(Omega)) ----
    grad_Omega = grad(gfOmega)
    mip = mesh(probe[0], probe[1], probe[2])
    Hi_num = float(grad_Omega[field_idx](mip))
    Hi_ana = 3.0 / (mu_r + 2.0) * H0
    err_pct = (Hi_num - Hi_ana) / Hi_ana * 100.0 if Hi_ana != 0 else 0.0
    _log(f"H{field_axis}_origin = {Hi_num:.6e} (analytical {Hi_ana:.6e}, "
         f"err {err_pct:+.2f}%)")

    # ---- Optional GMSH post-export ----
    msh_file = ""
    if output_msh:
        try:
            from gmsh_post_export import GmshPostExport
            post = GmshPostExport(mesh)
            B_field = Mu * grad_Omega
            post.add_vector_field("B", B_field)
            post.add_vector_field("H", grad_Omega)
            post.write(output_msh)
            msh_file = output_msh
            _log(f"wrote {output_msh}")
        except Exception as e:
            _log(f"GMSH post-export skipped: {e}")

    return {
        "vol_path": vol_path,
        "field_axis": field_axis,
        "mu_r": mu_r,
        "H0": H0,
        "fes_order": fes_order,
        "R_kelvin": R_kelvin,
        "probe": list(probe),
        "Hi_origin": Hi_num,
        "Hi_analytical": Hi_ana,
        "error_pct": err_pct,
        "ne": mesh.ne,
        "ndof": fes.ndof,
        "slaved_dofs": slaved,
        "kelvin_offset": list(offset),
        "dirichlet": dirichlet,
        "converged": True,
        "iterations": 1,
        "msh_file": msh_file,
    }


def build_argparser():
    """argparse factory shared by main() and notebook DesignSpec callers."""
    parser = argparse.ArgumentParser(
        description="Kelvin Benchmark: verify Cubit-Kelvin pipeline "
                    "against analytical magnetic-sphere-in-uniform-field.")
    parser.add_argument("--vol", required=True,
                        help="Netgen .vol mesh file (must have "
                             "magnetic/air/kelvin blocks and "
                             "sphere/kelvin_int/kelvin_ext sidesets)")
    parser.add_argument("--fes-order", type=int, default=2,
                        help="FES polynomial order (>=2 recommended for "
                             "Kelvin; default 2)")
    parser.add_argument("--mu-r", type=float, default=100.0,
                        help="Magnetic sphere relative permeability (default 100)")
    parser.add_argument("--H0", type=float, default=1.0,
                        help="Uniform external field amplitude [A/m] (default 1)")
    parser.add_argument("--field-axis", choices=("x", "y", "z"), default="z",
                        help="External field direction (default z)")
    parser.add_argument("--R-kelvin", type=float, default=0.20,
                        help="Kelvin sphere radius [m] (default 0.20)")
    parser.add_argument("--probe-x", type=float, default=1e-10)
    parser.add_argument("--probe-y", type=float, default=1e-10)
    parser.add_argument("--probe-z", type=float, default=1e-10)
    parser.add_argument("--msh-output", default=None,
                        help="Optional GMSH .msh post-export path")
    # calc_main auto-adds --output at runtime when missing; declare it
    # here too so the panel generator (which only sees this factory)
    # knows the calc accepts it (kubota 2026-05-30).
    parser.add_argument("--output", default="",
                        help="JSON output file (optional)")
    return parser


def main():
    parser = build_argparser()

    def run(args):
        return solve_kelvin_benchmark(
            vol_path=args.vol,
            mu_r=args.mu_r,
            H0=args.H0,
            field_axis=args.field_axis,
            fes_order=args.fes_order,
            R_kelvin=args.R_kelvin,
            probe=(args.probe_x, args.probe_y, args.probe_z),
            output_msh=args.msh_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
