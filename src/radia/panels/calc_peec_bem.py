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
        --peec-n-peri 16 \\
        --frequency 7000 --current 1.0 --coil-sigma 5.8e7 \\
        --vol workpiece.vol --wp-label wp_surface \\
        --sigma 5.8e7 --half-thickness 0.0125 --mu-r 1.0 \\
        --impedance-model sibc

Note on inductance: the L_coil in the output is the VACUUM-COIL
inductance only.  PEEC-BEM is a 1-way coupling (no Delta_L back-
reaction), so L_coil here is NOT the physically correct in-system
L when a magnetic / lossy workpiece is present.  Use
calc_fem_coilmesh.py for accurate in-system L.

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


def solve_peec_bem_forward(peec_step,
                            frequency, current, coil_sigma,
                            vol, wp_label, sigma, half_thickness, mu_r,
                            impedance_model, h1_order,
                            n_peri=16,
                            msh_output=""):
    """One-way forward PEEC->BEM pipeline.  No Delta_L back-reaction.

    Inductance reliability: the PEEC vacuum L_coil reported here is
    the COIL-ONLY inductance from the Loop-bundle solve.  Because PEEC-
    BEM is a 1-way coupling (BEM workpiece does NOT feed back into the
    coil via Delta_L), the L_coil from this script is **not the
    physically correct in-system inductance** when a magnetic / lossy
    workpiece is present.  Use ``calc_fem_coilmesh.py`` (FEM A-V)
    when an accurate in-system L is required; this PEEC-BEM script is
    optimised for P_wp.

    PEEC parameters:
      - ``n_peri``: number of perimeter-only filaments around the
        cross-section (thin-skin regime).  16 by default.
      - ``nwinc`` / ``nhinc`` (volume-grid filament counts) are NOT
        supported here -- PEEC-BEM uses perimeter placement only.

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

    # 1. STEP -> filaments (perimeter-only, n_peri filaments)
    progress("PEEC", f"STEP -> perimeter filaments (n_peri={n_peri})")
    t0 = time.perf_counter()
    topo = filaments_from_step(peec_step, sigma=coil_sigma,
                                n_peri=n_peri,
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
    # Backend selection (kwarg `bem_backend`):
    #   "ngsolve" : legacy lazy ngsolve.bem + GMRES (use_fmm-free)
    #   "intree"  : in-tree C++ Sauter-Schwab assembler + dense LU
    #   "hacapk"  : in-tree C++ + HACApK ACA + GMRES  (default; same
    #               accuracy, ~4x less memory, ~2x faster solve at N>2K)
    bem_backend = kwargs.get("bem_backend", "hacapk") if "kwargs" in dir() else "hacapk"
    # (calc_peec_bem doesn't accept **kwargs in its signature -- expose
    #  via global locals if added in future; for now hardcode "hacapk".)
    # Override with environment variable for benchmarking
    bem_backend = os.environ.get("RADIA_BEM_BACKEND", "hacapk")

    progress("BEM", f"assembly (order={h1_order}, backend={bem_backend})")
    t0 = time.perf_counter()
    if bem_backend == "ngsolve":
        bem = ScalarBIESIBCSolver(wp_mesh, order=h1_order,
                                    assemble_dense=False)
    elif bem_backend == "intree":
        bem = ScalarBIESIBCSolver(wp_mesh, order=h1_order,
                                    assemble_dense=True,
                                    use_intree_bem=True,
                                    intree_singular_n_q=6,
                                    intree_regular_quad_degree=7)
    elif bem_backend == "hacapk":
        bem = ScalarBIESIBCSolver(wp_mesh, order=h1_order,
                                    assemble_dense=True,
                                    use_intree_bem=True,
                                    intree_singular_n_q=6,
                                    intree_regular_quad_degree=7,
                                    use_intree_hacapk=True,
                                    hacapk_aca_eps=1e-10,
                                    hacapk_leaf=64,
                                    hacapk_eta=2.0)
    else:
        raise ValueError(f"unknown bem_backend={bem_backend!r}")
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
    mat = EMMaterial(name="custom", sigma=sigma, mu_r=mu_r)
    delta = mat.skin_depth(frequency)
    rho = 1.0 / sigma
    Z_s = (1.0 + 1j) * rho / delta * math.sqrt(mu_r)

    progress("BEM", f"BIE solve (Z_s model={impedance_model}, backend={bem_backend})")
    t0 = time.perf_counter()
    if bem_backend == "ngsolve":
        res = bem.solve_iterative(phi_inc, Z_s=Z_s, omega=omega,
                                    tol=1e-8, maxiter=300, restart=50)
    elif bem_backend == "intree":
        res = bem.solve(phi_inc, Z_s=Z_s, omega=omega)
    elif bem_backend == "hacapk":
        res = bem.solve_hacapk(phi_inc, Z_s=Z_s, omega=omega,
                                tol=1e-8, maxiter=500, restart=80)
    t_bie = time.perf_counter() - t0
    info = res.get('gmres_info', 0)
    progress("BEM", f"BIE ({t_bie:.1f}s, gmres info={info})")

    # 5. Report P over wp surface
    from ngsolve import (Integrate, CF, GridFunction, InnerProduct, grad)
    A_wp = float(Integrate(CF(1), wp_mesh, VOL_or_BND=BND).real)
    P_wp = float(res['P_density'] * A_wp)
    H_t_rms = float(res['H_t_rms'])

    # 5b. Surface heat-flux distribution q_surf(x,y,z) [W/m^2].
    # Mirrors calc_fem_kelvin.py Step 6b for optimization filtering:
    # the scalar P_wp answers "how much", the q_surf_max/p95/mean
    # answer "where".  |H_t|^2 = |grad_s phi_re|^2 + |grad_s phi_im|^2
    # is the same CF the bem.solve() integrand uses, just kept as a
    # CoefficientFunction so we can sample it pointwise.
    q_surf_max = None
    q_surf_mean = None
    q_surf_p95 = None
    P_total_check = None
    try:
        gf_phi_re = GridFunction(bem.fes)
        gf_phi_im = GridFunction(bem.fes)
        gf_phi_re.vec.FV().NumPy()[:] = res['phi_vec'].real
        gf_phi_im.vec.FV().NumPy()[:] = res['phi_vec'].imag
        H_t_sq_cf = (InnerProduct(grad(gf_phi_re), grad(gf_phi_re))
                     + InnerProduct(grad(gf_phi_im), grad(gf_phi_im)))
        Zs_re = float(Z_s.real)  # global SIBC: scalar
        q_surf_cf = 0.5 * Zs_re * H_t_sq_cf

        P_total_check = float(Integrate(q_surf_cf, wp_mesh,
                                        VOL_or_BND=BND).real)
        q_surf_mean = P_total_check / max(A_wp, 1e-30)

        # Stats via per-vertex sampling on the BEM H1 FES.
        gf_q = GridFunction(bem.fes)
        gf_q.Set(q_surf_cf, definedon=wp_mesh.Boundaries(".*"))
        q_arr = np.asarray(gf_q.vec.FV().NumPy())
        if q_arr.size > 0:
            q_surf_max = float(np.max(q_arr))
            q_surf_p95 = float(np.percentile(q_arr, 95))
        else:
            q_surf_max = 0.0
            q_surf_p95 = 0.0
        progress("BEM", f"Q_SURF max={q_surf_max:.3e} "
                        f"mean={q_surf_mean:.3e} p95={q_surf_p95:.3e} "
                        f"W/m^2 (P_check={P_total_check:.4e} vs "
                        f"P_wp={P_wp:.4e})")
    except Exception as e:
        gf_q = None
        progress("BEM", f"Q_SURF stats failed: {type(e).__name__}: {e}")

    # 5c. Surface current density on the wp face.
    #   - |J_s| [A/m] (scalar)  : intensity colormap.
    #   - Re(J_s) (3D vector)   : direction arrows (instantaneous-
    #     phase snapshot, real part of -grad_s(phi) in BEM scalar
    #     potential terms).
    gf_J = None
    gf_J_vec = None
    try:
        from ngsolve import sqrt as ng_sqrt, specialcf as _scf
        H_t_mag_cf = ng_sqrt(InnerProduct(grad(gf_phi_re), grad(gf_phi_re))
                              + InnerProduct(grad(gf_phi_im), grad(gf_phi_im)))
        gf_J = GridFunction(bem.fes)
        gf_J.Set(H_t_mag_cf, definedon=wp_mesh.Boundaries(".*"))
    except Exception as e:
        progress("BEM", f"J_SURF scalar gen failed: {type(e).__name__}: {e}")
    try:
        # Surface current density J_s for visualization (real-part
        # snapshot at instantaneous phase).
        #
        # Sign convention -- the hole approach matters here.  In this
        # .vol the wp is NOT meshed: only air + coil are present, and
        # the "sibc" boundary is the air-side cavity surface.  So
        # specialcf.normal(3) on the sibc bnd is the AIR-domain outward
        # normal -- it points FROM air INTO the (fictitious) wp, i.e.
        # OPPOSITE the standard "conductor outward normal" used in the
        # textbook formula J_s = n_cond x H_t.
        #
        # Let n = NGSolve normal (air-outward, into-wp).  Then
        #   n_cond = -n
        #   H_t   = -grad_s(phi)
        #   J_s   = n_cond x H_t  =  (-n) x (-grad_s phi)  =  +n x grad_s phi
        # so the correct expression with the NGSolve normal is
        #   J_s = +Cross(n, grad_s phi).
        # 2026-04-29 morning version had -Cross (got the textbook formula
        # right but missed the inverted normal), so all azimuthal arrows
        # pointed the wrong way.  Fixed evening 2026-04-29.
        from ngsolve import Cross
        n_bnd_v = _scf.normal(3)
        gphi = grad(gf_phi_re)
        gphi_dot_n = sum(gphi[i] * n_bnd_v[i] for i in range(3))
        gphi_t = CF(tuple(gphi[i] - gphi_dot_n * n_bnd_v[i]
                           for i in range(3)))
        J_s_re_vec = Cross(n_bnd_v, gphi_t)  # = n_air x grad_s(phi) = J_s
        from ngsolve import H1 as _H1
        fes_J_vec = _H1(wp_mesh, order=int(h1_order), dim=3)
        gf_J_vec = GridFunction(fes_J_vec)
        gf_J_vec.vec[:] = 0
        gf_J_vec.Set(J_s_re_vec, definedon=wp_mesh.Boundaries(".*"))
    except Exception as e:
        progress("BEM", f"J_SURF vector gen failed: {type(e).__name__}: {e}")
        gf_J_vec = None

    # 5d. GMSH .msh export with q_surf + J_surf as ElementData /
    # NodeData.  Saves the BEM workpiece SURFACE mesh (a 2D-in-3D
    # triangle mesh) so Kubota's Open GMSH button shows the
    # hotspot + current-intensity distribution coloured per surface
    # vertex.  Without this the exported .msh is just bare triangles
    # and the user reports "GMSH opens but I don't see anything".
    gmsh_file = ""
    qsurf_sol_path = ""
    if msh_output and gf_q is not None:
        try:
            from gmsh_post_export import save_vol_sol_pair, vol2msh
            base_dir = os.path.dirname(os.path.abspath(msh_output))
            stem = os.path.splitext(os.path.basename(msh_output))[0]
            vol_q = os.path.join(base_dir,
                                 f"{stem}_bem.vol").replace("\\", "/")
            sol_q = os.path.join(base_dir,
                                 f"{stem}_qsurf.sol").replace("\\", "/")
            save_vol_sol_pair(vol_q, sol_q, wp_mesh.ngmesh, gf_q)
            qsurf_sol_path = sol_q
            sol_entries = [
                {"sol": sol_q, "fes": "H1",
                 "fes_order": int(h1_order),
                 "fes_dim": 1,
                 "name": "q_surf", "ncomp": 1},
            ]
            if gf_J is not None:
                sol_J = os.path.join(base_dir,
                                     f"{stem}_Jsurf.sol").replace("\\", "/")
                gf_J.Save(sol_J)
                sol_entries.append(
                    {"sol": sol_J, "fes": "H1",
                     "fes_order": int(h1_order),
                     "fes_dim": 1,
                     "name": "J_surf_Am", "ncomp": 1})
            if gf_J_vec is not None:
                sol_Jv = os.path.join(base_dir,
                                      f"{stem}_Jvec.sol").replace("\\", "/")
                gf_J_vec.Save(sol_Jv)
                sol_entries.append(
                    {"sol": sol_Jv, "fes": "H1",
                     "fes_order": int(h1_order),
                     "fes_dim": 3,
                     "name": "J_surf_vec", "ncomp": 3})
            vol2msh(msh_output, vol_q, sol_entries)
            gmsh_file = msh_output
            progress("BEM", f"GMSH:wrote {os.path.basename(msh_output)} "
                             f"({len(sol_entries)} fields)")
        except Exception as e:
            progress("BEM", f"GMSH_ERROR:{type(e).__name__}: {e}")

    return {
        "status": "ok",
        "method": "PEEC+BEM (1-way, forward)",
        "frequency_hz": float(frequency),
        "current_A": float(current),
        "n_filaments": int(len(paths)),
        # Coil (PEEC vacuum, no back-reaction).
        # WARNING: this L_coil is the VACUUM-COIL inductance only --
        # PEEC-BEM is a 1-way coupling, so the BEM workpiece does NOT
        # feed back into the coil via Delta_L.  When a magnetic /
        # lossy workpiece is present this number is NOT the physically
        # correct in-system inductance.  Use ``calc_fem_coilmesh.py``
        # (FEM A-V) when an accurate in-system L is required.
        "L_coil_nH": float(L_coil * 1e9),
        "R_coil_mOhm": float(R_coil * 1e3),
        "L_coil_note": "VACUUM coil only -- PEEC-BEM is 1-way and L "
                        "ignores workpiece back-reaction.  For accurate "
                        "in-system L with a workpiece, use FEM A-V "
                        "(calc_fem_coilmesh.py).",
        "L_coil_reliability": "low",
        # Workpiece (BEM-SIBC, 1-way forward)
        "P_wp_W": P_wp,
        "H_t_rms_Am": H_t_rms,
        "wp_area_m2": float(A_wp),
        "skin_depth_mm": float(delta * 1e3),
        "Z_s_real": float(Z_s.real),
        "Z_s_imag": float(Z_s.imag),
        # Surface heat-flux distribution (optimization filter).
        # Same schema as calc_fem_kelvin.py for cross-comparison.
        "q_surf_max": q_surf_max,
        "q_surf_mean": q_surf_mean,
        "q_surf_p95": q_surf_p95,
        "P_total_check": P_total_check,
        "qsurf_sol": qsurf_sol_path,
        "msh_file": gmsh_file,
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
    parser.add_argument("--peec-n-peri", type=int, default=16,
                        help="Number of perimeter filaments around the "
                             "cross-section (thin-skin regime, d/delta>=3). "
                             "PEEC-BEM does NOT support volume-grid "
                             "filament counts (--peec-nwinc / --peec-nhinc) "
                             "-- perimeter placement only.")
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
            n_peri=args.peec_n_peri,
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
            msh_output=args.msh_output,
        )
        return result

    calc_main(run, parser)


if __name__ == "__main__":
    main()
