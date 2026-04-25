"""calc_peec_bem.py -- PEEC coil + BEM-SIBC workpiece (1-way, P-focus).

Layer 4 subprocess calc for the IH panel "PEEC+BEM SIBC (1-way)" method.

Pipeline (forward-only, no workpiece back-reaction on coil L):
  1. STEP -> filaments (CoilBuilder, profile-aware placement)
  2. PEEC bundle solve -> I_fil, V_port -> L_coil, R_coil (vacuum)
  3. Extract workpiece surface mesh from .vol (keep_label filter)
  4. BEM-SIBC solve on wp surface:
       phi_inc = Biot-Savart from filaments with I_fil
       SIBC BIE -> phi_s -> H_t, J_s, P
  5. Report L_coil + R_coil (PEEC vacuum) + P_wp + H_t_rms

Rationale: The back-reaction Delta_L from BEM is mesh-dependent
(point-source panel approximation + midpoint-rule line integral +
near-field log contamination).  For IH heating analysis P is the
primary output (5% accurate in this pipeline).  Use calc_fem_kelvin
with --peec-step for L to 1%.

Usage (from the IH panel):
    python calc_peec_bem.py --peec-step coil.step \\
        --peec-nwinc 3 --peec-nhinc 3 \\
        --frequency 7000 --current 1.0 --coil-sigma 5.8e7 \\
        --vol workpiece.vol --wp-label wp_surface \\
        --sigma 5.8e7 --half-thickness 0.0125 --mu-r 1.0 \\
        --impedance-model sibc

Output: JSON to stdout (calc_main contract).
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, ".."))
# bem_sibc_solver lives under examples/bem_reference/ (demoted from
# production per memory project_ih_bem_to_examples.md).  We import it
# from there; this keeps only one canonical copy.
BEM_REF = os.path.abspath(os.path.join(
    HERE, "..", "..", "..", "examples", "induction_heating", "bem_reference"))
for p in (SRC_RADIA, HERE, BEM_REF):
    if p not in sys.path:
        sys.path.insert(0, p)

from calc_common import calc_main, progress


def _extract_bnd_only(vol_mesh, bnd_label):
    """Extract a surface mesh containing only BND elements with the
    given boundary label.  Used for the workpiece-as-hole + sibc
    sideset geometry (no workpiece volume, closed-torus style).

    The same bc name often appears across MULTIPLE FaceDescriptors
    (e.g. the closed-torus sample has 4 separate FDs all named "sibc"
    because the wp cylinder is split by webcut + side + cap).  We
    therefore match by NAME and collect every FD index whose label
    equals bnd_label.  ``el.index`` for BND elements in NGSolve is
    0-based into the GetBoundaries() list (see
    calc_heating_bem._extract_surface_mesh_filtered).
    """
    from ngsolve import Mesh, BND
    import netgen.meshing as ngm

    bnd_labels = list(vol_mesh.GetBoundaries())
    target_indices = {i for i, n in enumerate(bnd_labels)
                       if n == bnd_label}
    if not target_indices:
        raise ValueError(
            f"Boundary label {bnd_label!r} not found in .vol. "
            f"Available: {sorted(set(bnd_labels))}")

    ngmesh_new = ngm.Mesh(dim=3)
    fd = ngmesh_new.Add(ngm.FaceDescriptor(bc=1, domin=1))
    ngmesh_new.SetBCName(0, bnd_label)

    used_vtx = set()
    for el in vol_mesh.Elements(BND):
        if el.index in target_indices:
            for v in el.vertices:
                used_vtx.add(v.nr)

    if not used_vtx:
        raise ValueError(
            f"No BND elements found with label {bnd_label!r} "
            f"(matched FD indices {sorted(target_indices)})")

    old_to_new = {}
    for old in sorted(used_vtx):
        p = vol_mesh.vertices[old].point
        old_to_new[old] = ngmesh_new.Add(
            ngm.MeshPoint(ngm.Pnt(p[0], p[1], p[2])))

    # When wp is a HOLE in the air volume (sibc sideset on air's bnd),
    # different FaceDescriptors may have their BND elements wound with
    # OPPOSITE orientations (e.g. air_top vs air_bot halves).  Mixing
    # inward and outward normals on the same extracted surface ruins
    # the BEM solve.  Robust fix: compute the wp centroid from the
    # collected vertex cloud, then for each triangle orient it so its
    # normal points AWAY from the centroid (= outward from the hole,
    # matching BEM-SIBC's wp-exterior-normal convention).
    coords_np = np.array([vol_mesh.vertices[v].point for v in sorted(used_vtx)])
    wp_centroid = coords_np.mean(axis=0)

    n_flipped = 0
    for el in vol_mesh.Elements(BND):
        if el.index not in target_indices:
            continue
        verts = [old_to_new[v.nr] for v in el.vertices]
        pts = np.array([vol_mesh.vertices[v.nr].point for v in el.vertices])
        if len(verts) >= 3:
            e1 = pts[1] - pts[0]
            e2 = pts[2] - pts[0]
            n = np.cross(e1, e2)
            outward = pts.mean(axis=0) - wp_centroid
            if np.dot(n, outward) < 0:
                # flip by reversing the non-first vertices
                if len(verts) == 3:
                    verts = [verts[0], verts[2], verts[1]]
                elif len(verts) == 4:
                    verts = [verts[0], verts[3], verts[2], verts[1]]
                n_flipped += 1
        ngmesh_new.Add(ngm.Element2D(fd, verts))

    progress("BEM", f"oriented {n_flipped} triangles outward from wp "
                    f"centroid {wp_centroid.tolist()}")
    return Mesh(ngmesh_new)


def solve_peec_bem_forward(peec_step, peec_nwinc, peec_nhinc,
                            frequency, current, coil_sigma,
                            vol, wp_label, sigma, half_thickness, mu_r,
                            impedance_model, h1_order):
    """One-way forward PEEC->BEM pipeline.  No Delta_L back-reaction.

    wp_label is interpreted as:
      - a material name first (extract BND of that volume), else
      - a boundary/sideset label (extract those BND elements directly).
    The latter is the 'workpiece-as-hole + sibc sideset' convention
    used by the closed-torus sample and calc_fem_coilmesh.
    """
    # Trigger Radia's MKL DLL path setup before peec_matrices loads.
    # Required because this script is launched as a subprocess by the
    # IH panel and by regression tests; neither context imports radia
    # ambiently.
    import radia  # noqa: F401

    from coil_from_cad import filaments_from_step
    from peec_bundle import (build_loop_bundle_impedance,
                              solve_loop_bundle)
    from ngsolve import Mesh, BND
    from em_material import EMMaterial
    from surface_mesh_extract import _extract_surface_mesh_filtered
    from bem_sibc_solver import (ScalarBIESIBCSolver,
                                  compute_phi_inc_from_filaments)

    omega = 2 * math.pi * frequency

    # 1. STEP -> filaments
    progress("PEEC", f"STEP -> filaments (nw={peec_nwinc}, nh={peec_nhinc})")
    t0 = time.perf_counter()
    topo = filaments_from_step(peec_step, sigma=coil_sigma,
                                nwinc=peec_nwinc, nhinc=peec_nhinc,
                                use_coil_builder=True)
    paths = topo["filament_paths"]
    seg_of_fil = topo["seg_of_filament"]
    solver = topo["solver"]
    t_topo = time.perf_counter() - t0
    progress("PEEC", f"{len(paths)} filaments, {t_topo:.1f}s")

    # 2. PEEC bundle solve (Loop form)
    progress("PEEC", f"Loop-bundle solve @ {frequency:.0f} Hz")
    t0 = time.perf_counter()
    R_f, L_f = build_loop_bundle_impedance(solver, seg_of_fil)
    I_fil, V_port = solve_loop_bundle(R_f, L_f, frequency, I_port=current)
    t_peec = time.perf_counter() - t0
    Z_coil = V_port / current
    L_coil = Z_coil.imag / omega if omega > 0 else 0.0
    R_coil = Z_coil.real
    progress("PEEC", f"L_coil={L_coil*1e9:.2f} nH, R_coil={R_coil*1e3:.4f} mOhm"
                      f" ({t_peec:.1f}s)")

    # 3. Workpiece surface mesh from .vol
    progress("BEM", f"load wp surface from {os.path.basename(vol)}")
    t0 = time.perf_counter()
    vol_mesh = Mesh(vol)
    mats = set(vol_mesh.GetMaterials())
    bnds = set(vol_mesh.GetBoundaries())
    if wp_label in mats:
        progress("BEM", f"wp_label {wp_label!r} found as material")
        wp_mesh = _extract_surface_mesh_filtered(vol_mesh, keep_label=wp_label)
    elif wp_label in bnds:
        progress("BEM", f"wp_label {wp_label!r} found as boundary sideset")
        wp_mesh = _extract_bnd_only(vol_mesh, wp_label)
    else:
        raise ValueError(
            f"wp_label {wp_label!r} is neither a material "
            f"({sorted(mats)}) nor a boundary "
            f"({sorted(bnds)}) of {vol}")
    t_mesh = time.perf_counter() - t0
    progress("BEM", f"wp nv={wp_mesh.nv} ne(BND)={wp_mesh.GetNE(BND)}"
                    f" ({t_mesh:.1f}s)")

    # 4. BEM-SIBC assembly + solve
    progress("BEM", f"assembly (order={h1_order})")
    t0 = time.perf_counter()
    bem = ScalarBIESIBCSolver(wp_mesh, order=h1_order)
    t_asm = time.perf_counter() - t0
    progress("BEM", f"ndof={bem.ndof} ({t_asm:.1f}s)")

    # phi_inc from filaments with I_fil
    obs = np.array([[wp_mesh.vertices[i].point[j] for j in range(3)]
                     for i in range(wp_mesh.nv)])
    progress("BEM", "phi_inc from filaments")
    t0 = time.perf_counter()
    phi_inc = compute_phi_inc_from_filaments(obs, paths, I_fil)
    t_phi = time.perf_counter() - t0
    progress("BEM", f"phi_inc ({t_phi:.1f}s)")

    # Impedance model: linear (Z_s = (1+j)*rho/delta) or ESIM (nonlinear).
    # 1-way forward = single BEM solve.  ESIM Karl iteration is NOT
    # run here (that is a workpiece-internal self-consistency, still
    # 1-way wrt coil).  For the 1-way-P-only scope we use the linear
    # SIBC only.  ESIM users should pick the FEM method.
    mat = EMMaterial(name="custom", sigma=sigma, mu_r=mu_r)
    delta = mat.skin_depth(frequency)
    rho = 1.0 / sigma
    Z_s = (1.0 + 1j) * rho / delta * math.sqrt(mu_r)

    progress("BEM", f"BIE solve (Z_s model={impedance_model})")
    t0 = time.perf_counter()
    res = bem.solve(phi_inc, Z_s=Z_s, omega=omega)
    t_bie = time.perf_counter() - t0
    progress("BEM", f"BIE ({t_bie:.1f}s)")

    # 5. Report P over wp surface
    from ngsolve import Integrate, CF
    A_wp = float(Integrate(CF(1), wp_mesh, VOL_or_BND=BND).real)
    P_wp = float(res['P_density'] * A_wp)
    H_t_rms = float(res['H_t_rms'])

    return {
        "status": "ok",
        "method": "PEEC+BEM (1-way, forward)",
        "frequency_hz": float(frequency),
        "current_A": float(current),
        "n_filaments": int(len(paths)),
        # Coil (PEEC vacuum, no back-reaction)
        "L_coil_nH": float(L_coil * 1e9),
        "R_coil_mOhm": float(R_coil * 1e3),
        "L_coil_note": "vacuum coil, no workpiece back-reaction "
                        "(use PEEC+FEM for L with back-reaction)",
        # Workpiece (BEM-SIBC, 1-way forward)
        "P_wp_W": P_wp,
        "H_t_rms_Am": H_t_rms,
        "wp_area_m2": float(A_wp),
        "skin_depth_mm": float(delta * 1e3),
        "Z_s_real": float(Z_s.real),
        "Z_s_imag": float(Z_s.imag),
        # Diagnostics
        "bem_ndof": int(bem.ndof),
        "bem_nv": int(wp_mesh.nv),
        "bem_ne": int(wp_mesh.GetNE(BND)),
        "t_topology_s": float(t_topo),
        "t_peec_solve_s": float(t_peec),
        "t_bem_assembly_s": float(t_asm),
        "t_bem_solve_s": float(t_bie),
    }


def main():
    parser = argparse.ArgumentParser(
        description="PEEC+BEM SIBC 1-way forward coupling (P-focus)")
    # PEEC coil
    parser.add_argument("--peec-step", required=True,
                        help="STEP file for PEEC coil")
    parser.add_argument("--peec-nwinc", type=int, default=3,
                        help="Filament sub-cells width direction")
    parser.add_argument("--peec-nhinc", type=int, default=3,
                        help="Filament sub-cells height direction")
    parser.add_argument("--frequency", type=float, required=True,
                        help="Frequency [Hz]")
    parser.add_argument("--current", type=float, default=1.0,
                        help="Port current [A]")
    parser.add_argument("--coil-sigma", type=float, default=5.8e7,
                        help="Coil conductivity [S/m]")
    # Workpiece
    parser.add_argument("--vol", required=True,
                        help=".vol file containing workpiece mesh")
    parser.add_argument("--wp-label", default="sibc",
                        help="Workpiece label: volume material name OR "
                             "boundary sideset label (wp-as-hole + sibc)")
    parser.add_argument("--sigma", type=float, required=True,
                        help="Workpiece conductivity [S/m]")
    parser.add_argument("--half-thickness", type=float, default=0.01,
                        help="Half-thickness for ESIM (unused in linear)")
    parser.add_argument("--mu-r", type=float, default=1.0,
                        help="Workpiece relative permeability")
    parser.add_argument("--impedance-model", default="sibc",
                        choices=["sibc", "esim"],
                        help="sibc: linear Dowell (production). "
                             "esim: Karl iteration (WIP, raises NotImplemented).")
    parser.add_argument("--peec-solver", default="dense",
                        choices=["dense", "hacapk"],
                        help="PEEC bundle linear solver.")
    parser.add_argument("--bh-file", default="",
                        help="BH table for ESIM (WIP).")
    parser.add_argument("--esim-max-iter", type=int, default=15,
                        help="ESIM Karl iteration max (WIP).")
    parser.add_argument("--esim-tol", type=float, default=1e-3,
                        help="ESIM Karl iteration tol (WIP).")
    parser.add_argument("--h1-order", type=int, default=1,
                        help="H1 polynomial order for BEM")
    parser.add_argument("--msh-output", default="",
                        help="Optional GMSH .msh output path. When set, "
                             "the workpiece .vol mesh is converted to "
                             ".msh after solve so the panel's OpenGmsh "
                             "button can view the geometry.")

    def run(args):
        if args.impedance_model == "esim":
            return {
                "error": "ESIM is not implemented in calc_peec_bem yet. "
                         "Use --impedance-model sibc for now. For "
                         "nonlinear BH, use the FEM A-V panel method "
                         "(also WIP for ESIM)."
            }
        result = solve_peec_bem_forward(
            peec_step=args.peec_step,
            peec_nwinc=args.peec_nwinc,
            peec_nhinc=args.peec_nhinc,
            frequency=args.frequency,
            current=args.current,
            coil_sigma=args.coil_sigma,
            vol=args.vol,
            wp_label=args.wp_label,
            sigma=args.sigma,
            half_thickness=args.half_thickness,
            mu_r=args.mu_r,
            impedance_model=args.impedance_model,
            h1_order=args.h1_order,
        )
        # GMSH export (mesh geometry only — minimum to activate panel's
        # OpenGmsh button). Field export (H_t, P_density on wp surface)
        # is a future enhancement.
        if args.msh_output and isinstance(result, dict) and "error" not in result:
            try:
                import sys, os as _os
                radia_src = _os.path.dirname(_os.path.abspath(__file__)) + "/.."
                if _os.path.abspath(radia_src) not in sys.path:
                    sys.path.insert(0, _os.path.abspath(radia_src))
                from gmsh_post_export import vol2msh
                vol2msh(args.msh_output, args.vol, [])
                result["msh_file"] = args.msh_output
            except Exception as e:
                result["msh_export_error"] = str(e)
        return result

    calc_main(run, parser)


if __name__ == "__main__":
    main()
