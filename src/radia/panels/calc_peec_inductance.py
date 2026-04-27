"""calc_peec_inductance.py -- PEEC coil inductance from STEP (Layer 4, no GUI).

Layer 4 subprocess calc for the IH panel "PEEC inductance (coil only, STEP)"
method.  No workpiece, no BEM, no FEM mesh.  Just STEP -> filaments -> L, R.

Pipeline:
  1. STEP -> filaments (CoilBuilder, perimeter-only placement)
  2. PEEC Loop-bundle solve at given frequency
  3. Report L_coil, R_coil (vacuum)

Filaments are placed on the cross-section PERIMETER only (thin-skin limit,
d / delta >= 3 — the typical IH operating regime).  Use n_peri filaments
spread around the arc-length perimeter of each cross-section.

Usage (from the IH panel):
    python calc_peec_inductance.py --peec-step coil.step \\
        --peec-n-peri 16 \\
        --frequency 7000 --current 1.0 --coil-sigma 5.8e7

Output: JSON to stdout (calc_main contract).
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, ".."))
for p in (SRC_RADIA, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from calc_common import calc_main, progress


def _build_field_domain_mesh(filament_paths, expand=2.5, n_cells=18):
    """Generate an OCC Box mesh enclosing the coil bbox * ``expand``.

    Used for the auto-field-domain visualisation in PEEC-inductance
    mode (no user-supplied .vol).  Sized so a typical coil fits in
    the centre with ~``n_cells`` divisions in each direction; the
    mesh is dense enough that GMSH renders the field smoothly but
    coarse enough that B-field evaluation runs in <10s.
    """
    import numpy as np
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh

    # filament_paths[k] is a list of segments [(p1, p2), ...]; each
    # endpoint is an (x, y, z) tuple.  Flatten to a (n_pts, 3) view
    # so we can take the AABB.
    pts_list = []
    for path in filament_paths:
        arr = np.asarray(path, dtype=float)  # (n_seg, 2, 3)
        if arr.ndim != 3 or arr.shape[-1] != 3:
            continue
        pts_list.append(arr.reshape(-1, 3))
    if not pts_list:
        raise ValueError(
            "_build_field_domain_mesh: no valid filament paths.")
    pts = np.concatenate(pts_list, axis=0)
    pmin = pts.min(axis=0)
    pmax = pts.max(axis=0)
    centre = 0.5 * (pmin + pmax)
    half = 0.5 * (pmax - pmin)
    # Expand by ``expand`` factor; keep a floor so a planar coil
    # (height ~ 0) still has some thickness in z.
    half_exp = np.maximum(half * expand, half.max() * 0.3)
    p_lo = centre - half_exp
    p_hi = centre + half_exp

    # netgen.occ.Pnt requires native Python floats; numpy scalars
    # raise "Invoked with: array(...)" on the C++ binding.
    pnt_lo = Pnt(float(p_lo[0]), float(p_lo[1]), float(p_lo[2]))
    pnt_hi = Pnt(float(p_hi[0]), float(p_hi[1]), float(p_hi[2]))
    box = Box(pnt_lo, pnt_hi)
    box.faces.name = "outer"
    box.solids.name = "field_domain"
    geo = OCCGeometry(box)
    maxh = float(np.max(p_hi - p_lo)) / max(n_cells, 4)
    ngm = geo.GenerateMesh(maxh=maxh)
    return Mesh(ngm)


def _biot_savart_B_on_mesh(mesh, filament_paths, fil_currents):
    """Vectorised Biot-Savart evaluation of B at every vertex of
    ``mesh``.  Returns an (n_vertices, 3) NumPy array of B [T].
    """
    import numpy as np
    from biot_savart import h_segments_batch, MU0

    # Vertex coordinates as (n, 3).
    verts = np.asarray([list(v.point) for v in mesh.vertices], dtype=float)
    if verts.shape[1] == 2:
        # 2D mesh would be a programming error here -- this helper is
        # called from PEEC-inductance which always builds a 3D Box.
        raise ValueError(
            "_biot_savart_B_on_mesh: 2D mesh not supported.")

    H_total = np.zeros((verts.shape[0], 3), dtype=float)
    for path, I in zip(filament_paths, fil_currents):
        # PEEC paths are ALREADY a list of segments (p1, p2 pairs);
        # h_segments_batch expects exactly that, so pass the path
        # straight through (no consecutive-pair zip).
        segments = np.asarray(path, dtype=float)
        if segments.ndim != 3 or segments.shape[1:] != (2, 3):
            continue
        # PEEC currents are complex phasors; visualise the magnitude
        # of the time-harmonic B (peak amplitude) by feeding |I| to
        # the Biot-Savart kernel.  Direction follows the geometry of
        # the filament path -- identical to the DC current pattern.
        Im = float(abs(complex(I)))
        if Im < 1e-30:
            continue
        H_total += h_segments_batch(segments, verts, current=Im)
    return MU0 * H_total  # B [T]


def solve_peec_inductance(peec_input, n_peri,
                          frequency, current, coil_sigma,
                          msh_output=""):
    """Build PEEC topology from STEP or JOU input, then solve for L, R.

    Routing:
      - `.jou`            -> `coil_from_jou.filaments_from_jou` (explicit
                              centerline from `move Surface` commands).
                              Use for multi-turn / 3turnCoil-class coils
                              where STEP walker heuristics fail.
      - `.step` / `.stp`  -> `coil_from_cad.filaments_from_step` (walker
                              centerline extraction).  Use for clean
                              single-loop torus / helix coils.
      - other             -> error.
    """
    # Trigger Radia's MKL DLL path setup before peec_matrices loads.
    import radia  # noqa: F401

    from peec_bundle import (build_loop_bundle_impedance,
                              solve_loop_bundle)

    omega = 2 * math.pi * frequency

    ext = os.path.splitext(peec_input)[1].lower()

    # Auto-prefer sibling .jou (exact stem match, case-insensitive).
    # When user happens to have both foo.step and foo.jou in the same
    # directory — the common Cubit export pattern, since the panel's
    # ensure_jou_path() saves the .jou before every STEP export — we
    # use the .jou for correct L on multi-turn lofts.  The STEP
    # longest-edge path currently mis-estimates cross-section area on
    # tight-pancake multi-turn lofts (Kubota's 3turncoil.stp: STEP
    # path L=4.8 nH WRONG vs sibling 3turnCoil.jou L=426 nH correct).
    # Fixing that is a cross-section-geometry bug (see TODO below);
    # meanwhile sibling-match gives the user the right answer when
    # they have the right file.
    if ext in (".step", ".stp"):
        peec_dir = os.path.dirname(peec_input) or "."
        base = os.path.splitext(os.path.basename(peec_input))[0]
        try:
            entries = os.listdir(peec_dir)
        except OSError:
            entries = []
        jou_sibling = None
        for name in entries:
            stem, e = os.path.splitext(name)
            if e.lower() == ".jou" and stem.lower() == base.lower():
                jou_sibling = os.path.join(peec_dir, name)
                break
        if jou_sibling is not None:
            # Only prefer the sibling if it actually contains the
            # PEEC explicit-centerline pattern (`move Surface N x Y y Y z Z`).
            # A generic Cubit .jou (rect sweep, brick mesh script, etc.)
            # has no centerline and would crash the .jou parser.  In that
            # case fall through to the STEP path.
            is_peec_jou = False
            try:
                with open(jou_sibling, encoding="utf-8",
                          errors="replace") as jf:
                    for line in jf:
                        # Avoid importing the parser here; check the regex
                        # inline.  Match `move Surface <int> x <n> y <n> z <n>`.
                        s = line.strip()
                        if (s.startswith("move Surface ")
                                and " x " in s and " y " in s
                                and " z " in s):
                            is_peec_jou = True
                            break
            except OSError:
                pass
            if is_peec_jou:
                progress("PEEC", f"found sibling .jou, preferring it: "
                                  f"{os.path.basename(jou_sibling)}")
                peec_input = jou_sibling
                ext = ".jou"
            else:
                progress("PEEC", f"sibling .jou {os.path.basename(jou_sibling)} "
                                  f"has no explicit centerline — using STEP")

    if ext == ".jou":
        from coil_from_jou import filaments_from_jou
        source_kind = "JOU"
        progress("PEEC", f"JOU -> explicit centerline (n_peri={n_peri})")
        t0 = time.perf_counter()
        topo = filaments_from_jou(peec_input, sigma=coil_sigma,
                                    n_peri=n_peri)
    elif ext in (".step", ".stp"):
        from coil_from_cad import filaments_from_step
        source_kind = "STEP"
        progress("PEEC", f"STEP -> perimeter filaments (n_peri={n_peri})")
        t0 = time.perf_counter()
        topo = filaments_from_step(peec_input, sigma=coil_sigma,
                                    n_peri=n_peri,
                                    use_coil_builder=True)
    else:
        raise ValueError(
            f"Unsupported input extension {ext!r}. "
            f"Expected .jou (Cubit journal with explicit centerline) or "
            f".step / .stp (geometry file for walker extraction).")

    paths = topo["filament_paths"]
    seg_of_fil = topo["seg_of_filament"]
    solver = topo["solver"]
    t_topo = time.perf_counter() - t0
    progress("PEEC", f"{len(paths)} filaments via {source_kind}, {t_topo:.1f}s")

    progress("PEEC", f"Loop-bundle solve @ {frequency:.0f} Hz")
    t0 = time.perf_counter()
    R_f, L_f = build_loop_bundle_impedance(solver, seg_of_fil)
    I_fil, V_port = solve_loop_bundle(R_f, L_f, frequency, I_port=current)
    t_peec = time.perf_counter() - t0
    Z_coil = V_port / current
    L_coil = Z_coil.imag / omega if omega > 0 else 0.0
    R_coil = Z_coil.real
    progress("PEEC", f"L_coil={L_coil*1e9:.3f} nH, "
                      f"R_coil={R_coil*1e3:.4f} mOhm ({t_peec:.1f}s)")

    # B-field visualisation on an auto-generated bounding-box mesh.
    # PEEC-inductance has no .vol input; we fabricate a 3D Box around
    # the coil so GMSH gets something to render.  Optional: only
    # runs when --msh-output is supplied.
    gmsh_file = ""
    t_field = 0.0
    field_n_vertices = 0
    if msh_output:
        try:
            from ngsolve import H1, GridFunction
            from gmsh_post_export import save_vol_sol_pair, vol2msh
            t_field_start = time.perf_counter()
            field_mesh = _build_field_domain_mesh(paths)
            field_n_vertices = field_mesh.nv
            B_arr = _biot_savart_B_on_mesh(field_mesh, paths, I_fil)

            fes_B = H1(field_mesh, order=1, dim=3)
            gf_B = GridFunction(fes_B)
            v = gf_B.vec.FV().NumPy()
            v[:] = B_arr.reshape(-1)

            base_dir = os.path.dirname(os.path.abspath(msh_output))
            stem = os.path.splitext(os.path.basename(msh_output))[0]
            sol_B = os.path.join(base_dir,
                                 f"{stem}_B.sol").replace("\\", "/")
            vol_B = os.path.join(base_dir,
                                 f"{stem}_field.vol").replace("\\", "/")
            save_vol_sol_pair(vol_B, sol_B, field_mesh.ngmesh, gf_B)
            sol_entries = [
                {"sol": sol_B, "fes": "H1",
                 "fes_order": 1, "fes_dim": 3,
                 "name": "B", "ncomp": 3},
            ]
            vol2msh(msh_output, vol_B, sol_entries)
            gmsh_file = msh_output
            t_field = time.perf_counter() - t_field_start
            progress("PEEC", f"GMSH:wrote {os.path.basename(msh_output)} "
                              f"(B field on {field_n_vertices} verts, "
                              f"{t_field:.1f}s)")
        except Exception as e:
            progress("PEEC", f"GMSH_ERROR:{type(e).__name__}: {e}")

    return {
        "status": "ok",
        "method": "PEEC inductance (coil only, STEP)",
        "input_kind": source_kind,
        "input_path": peec_input,
        "placement": "perimeter",
        "n_peri": int(n_peri),
        "frequency_hz": float(frequency),
        "current_A": float(current),
        "n_filaments": int(len(paths)),
        "coil_sigma_Sm": float(coil_sigma),
        "L_coil_H": float(L_coil),
        "L_coil_nH": float(L_coil * 1e9),
        "R_coil_Ohm": float(R_coil),
        "R_coil_mOhm": float(R_coil * 1e3),
        "impedance_real": float(Z_coil.real),
        "impedance_imag": float(Z_coil.imag),
        "abs_Z": float(abs(Z_coil)),
        "inductance_H": float(L_coil),
        "t_topology_s": float(t_topo),
        "t_peec_solve_s": float(t_peec),
        "t_field_s": float(t_field),
        "t_solve_s": float(t_topo + t_peec + t_field),
        "field_n_vertices": int(field_n_vertices),
        "msh_file": gmsh_file,
    }


def main():
    parser = argparse.ArgumentParser(
        description="PEEC coil inductance from STEP (vacuum, no workpiece)")
    parser.add_argument("--peec-step", required=True,
                        help="Input file for PEEC coil: .step / .stp "
                             "(walker extraction) or .jou (explicit "
                             "centerline from Cubit journal). The arg "
                             "name stays --peec-step for backward "
                             "compatibility; extension chooses the path.")
    parser.add_argument("--peec-n-peri", type=int, default=16,
                        help="Number of filaments on the cross-section "
                             "perimeter (thin-skin regime; d/delta>=3)")
    parser.add_argument("--frequency", type=float, required=True,
                        help="Frequency [Hz]")
    parser.add_argument("--current", type=float, default=1.0,
                        help="Port current [A]")
    parser.add_argument("--coil-sigma", type=float, default=5.8e7,
                        help="Coil conductivity [S/m]")
    parser.add_argument("--msh-output", default="",
                        help="Optional .msh path; when set, B field on "
                             "an auto-generated bounding-box mesh is "
                             "evaluated via Biot-Savart and exported "
                             "for the panel's Open GMSH button.")

    def run(args):
        return solve_peec_inductance(
            peec_input=args.peec_step,
            n_peri=args.peec_n_peri,
            frequency=args.frequency,
            current=args.current,
            coil_sigma=args.coil_sigma,
            msh_output=args.msh_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
