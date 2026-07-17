"""calc_inductance.py -- unified coil inductance (vacuum + weak-coupled).

Layer 4 subprocess for the IH "Inductance" panel.  Replaces three
scoped CLIs that overlapped on workflow but differed on coil solver:

  * calc_peec_inductance.py        -> --coil-solver peec, --vol omitted
  * calc_peec_bem.py               -> --coil-solver peec, --vol present
  * calc_coil_bem_a_workpiece.py   -> --coil-solver bem-a, --vol present

Coil-input format (one per solver):
  * PEEC   -> --coil-step (CAD STEP; perimeter filaments derived from
              centerline)
  * BEM-A  -> --coil-vol  (pre-meshed surface .vol from Cubit / Netgen;
              must contain 'source'/'sink' boundary labels on the cap
              faces)

Coil solver (``--coil-solver``):
  * ``peec`` : PEEC perimeter filaments + Loop-bundle solve.
               Fast (~10x BEM-A); n_peri-discretisation-limited (~3% under
               on rect cross-section at n_peri=16 vs converged BEM-A).
  * ``bem-a``: intree Weggler stabilized EFIE on HDivSurface RWG
               (Phase C.1-C.5).  Surface-current physics, n_peri-free.
               Reads a pre-meshed coil .vol; no on-the-fly OCC re-mesh.

Workpiece is OPTIONAL.  Without ``--vol`` the script returns only the
vacuum coil L_coil + R_coil.  With ``--vol`` it adds weak-coupled BEM-SIBC:
gauge-invariant Telegen φ·(n·B_inc) ΔL captures the workpiece's
back-reaction at the port level (coil current distribution is NOT
recomputed — that is "strong coupling", available via FEM A-V in
calc_fem_coilmesh.py).

Coupling terminology (corrected 2026-05-02):
  * **weak coupling** = back-reacted Telegen ΔL is computed; coil source
    distribution FIXED.  Default and only mode here.
  * **one-way (forward)** = only P_wp; no ΔL.  Not exposed: weak is a
    strict superset at no extra cost.
  * **strong coupling** = coil source recomputed iteratively in response
    to workpiece reflection.  Use calc_fem_coilmesh.py.

Output: JSON to stdout (calc_main contract).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np


def _detect_vol_curving_order(vol_path: str) -> int:
    """Read companion ``<vol_path>.json`` (written by Cubit's
    ``export netgen``) to find the curving order baked into the
    .vol.  Returns the integer ``order`` field, or 1 if the companion
    JSON is missing / malformed (= flat, safe default).

    This is the SOURCE OF TRUTH for the workpiece BEM basis order in
    ``calc_inductance.py``: the .vol's ``curvedelements`` is fixed at
    Cubit-export time, and trying to upgrade via ``mesh.Curve(p)``
    after-the-fact silently falls back to flat without a CAD callback.
    """
    json_path = vol_path + ".json"
    if not os.path.isfile(json_path):
        return 1
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        order = int(data.get("order", 1))
        return max(1, order)
    except (OSError, ValueError, KeyError):
        return 1


HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, ".."))
for p in (SRC_RADIA, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from calc_common import calc_main, progress  # noqa: E402


MU_0 = 4.0e-7 * math.pi
INV_4PI = 1.0 / (4.0 * math.pi)


# ======================================================================
# Coil solver: PEEC
# ======================================================================
def _equivalent_wire_radius_m(topo, n_peri):
    """Derive an equivalent-circle wire radius [m] from the PEEC topology.

    Priority:
      1. ``cross_section_radius_m_mean`` (circle-edge / round-wire path)
      2. ``cross_section_area_m2_mean`` -> sqrt(A / pi)
      3. ``cell_wh`` mean (last resort): each filament's (w, h) area
         times n_peri gives the bundle area, then sqrt(A_bundle / pi).
    Returns 0.0 if none of the above are available -- the caller falls
    back to the R_DC-only behaviour.
    """
    r = float(topo.get("cross_section_radius_m_mean") or 0.0)
    if r > 0:
        return r
    A = float(topo.get("cross_section_area_m2_mean") or 0.0)
    if A > 0:
        return math.sqrt(A / math.pi)
    cell_wh = topo.get("cell_wh") or []
    if cell_wh:
        areas = [float(w) * float(h)
                 for fil in cell_wh for (w, h) in fil]
        if areas:
            A_bundle = (sum(areas) / len(areas)) * int(n_peri)
            return math.sqrt(A_bundle / math.pi)
    return 0.0


def _peec_skin_impedance_per_filament(paths, wire_radius_m, sigma, omega,
                                       n_peri):
    """Per-filament complex skin-impedance ``Zs_fil`` for ``solve_loop_bundle``.

    Derivation (BEM-A analogue, see docs/esim/R_MISMATCH_PEEC_VS_BEMA.md).
    A round-wire bundle of n_peri parallel filaments has:

        Z_port(omega) = Z_cyl(omega) * L_filament

    where ``Z_cyl(omega)`` is the closed-form per-unit-length Bessel
    impedance of a SOLID cylinder (analytical_formulas.conductor_impedance.
    cylinder_ac_impedance).  In the parallel-filament solve each
    filament's self-impedance contributes as Z_fil / n_peri (parallel
    of n_peri identical filaments).  Combined with the PEEC builder's
    R_DC_per_filament = n_peri * rho * L / (pi a^2) the self diagonal
    should reach::

        Z_self_k = n_peri * Z_cyl(omega) * L_k

    so the additive correction over the existing R_DC term is::

        Zs_fil_k = n_peri * (Z_cyl(omega) - R_DC_per_unit_length) * L_k.

    This recovers the full Bessel impedance smoothly across frequencies:
    at omega=0 we get R_DC; at high omega we get (1+j)/(sigma*delta*P) per
    unit length (the BEM-A SIBC limit).
    """
    import numpy as np
    from radia.analytical_formulas.conductor_impedance import (
        cylinder_ac_impedance, cylinder_dc_resistance,
    )

    if wire_radius_m <= 0 or omega <= 0:
        return None
    Z_ac = cylinder_ac_impedance(a=wire_radius_m, sigma=sigma, omega=omega)
    R_dc = cylinder_dc_resistance(a=wire_radius_m, sigma=sigma)
    dZ_per_m = complex(Z_ac) - complex(R_dc, 0.0)  # AC-skin contribution / m
    Zs_fil = []
    for path_k in paths:
        # filament length = sum of |p2 - p1| over polyline tuples
        L_k = 0.0
        for (p1, p2) in path_k:
            d = np.asarray(p2, dtype=float) - np.asarray(p1, dtype=float)
            L_k += float(np.linalg.norm(d))
        Zs_fil.append(int(n_peri) * dZ_per_m * L_k)
    return np.asarray(Zs_fil, dtype=complex)


def _solve_coil_peec(args):
    """PEEC coil layer.

    Returns a dict with:
      * L_coil [H], R_coil [Ω] at unit terminal current.
        R is the full Bessel AC resistance for a round-wire bundle:
        per-filament ``Zs_fil = n_peri * (Z_cyl(omega) - R_DC_per_m) * L_k``
        is injected into ``solve_loop_bundle`` so the loop-bundle
        impedance approaches ``Z_cyl(omega) * L_filament`` at all
        frequencies (R_DC at omega=0, SIBC asymptote at high omega).
        Falls back to R_DC-only when the wire radius cannot be inferred
        from the PEEC topology (rectangular cross-sections with no
        cross_section_radius_m_mean / cross_section_area_m2_mean keys);
        in that case a WARNING is logged on the progress channel.
      * source_type = "filament"
      * paths : list of filament polylines (m)
      * I_fil : (K,) complex per-filament current at args.current
      * t_coil_s, n_filaments
    """
    progress("PEEC", f"STEP -> perimeter filaments (n_peri={args.peec_n_peri})")
    from radia.coil_from_cad import filaments_from_step
    from radia.peec_bundle import build_loop_bundle_impedance, solve_loop_bundle

    omega = 2.0 * math.pi * args.frequency
    t0 = time.perf_counter()
    topo = filaments_from_step(args.coil_step,
                               sigma=args.coil_sigma,
                               n_peri=args.peec_n_peri,
                               use_coil_builder=True)
    paths = topo["filament_paths"]
    seg_of_fil = topo["seg_of_filament"]
    solver = topo["solver"]
    t_topo = time.perf_counter() - t0
    progress("PEEC", f"{len(paths)} filaments ({t_topo:.1f}s)")

    # Compute per-filament Bessel skin impedance Zs_fil so the bundle
    # impedance reaches the proper round-wire AC value at args.frequency.
    a_eq = _equivalent_wire_radius_m(topo, args.peec_n_peri)
    Zs_fil = _peec_skin_impedance_per_filament(
        paths, a_eq, args.coil_sigma, omega, args.peec_n_peri)
    if Zs_fil is None:
        progress("PEEC", "WARNING: wire radius unknown -- returning R_DC "
                          "(no Bessel skin correction).")
    else:
        from radia.analytical_formulas.conductor_impedance import (
            skin_depth as _skin_depth,
        )
        delta_m = _skin_depth(args.coil_sigma, omega)
        progress("PEEC",
            f"Bessel skin: a_eq={a_eq * 1e3:.3f} mm, "
            f"delta={delta_m * 1e3:.3f} mm, "
            f"a/delta={a_eq / delta_m:.2f}")

    use_prox = bool(getattr(args, "peec_proximity", True)) and a_eq > 0 \
        and omega > 0
    progress("PEEC",
        f"Loop-bundle solve @ {args.frequency:.0f} Hz "
        f"({'proximity-iterative' if use_prox else 'self-only Bessel'})")
    t0 = time.perf_counter()
    R_f, L_f = build_loop_bundle_impedance(solver, seg_of_fil)

    proximity_info = None
    if use_prox:
        from radia.peec_proximity import solve_proximity_iterative
        prox = solve_proximity_iterative(
            R_f=R_f, L_f=L_f, filament_paths=paths,
            frequency=args.frequency,
            sigma=args.coil_sigma, wire_radius_m=a_eq,
            n_peri=int(args.peec_n_peri),
            I_port=args.current,
            max_iter=int(args.peec_proximity_max_iter),
            tol=float(args.peec_proximity_tol),
            relax=float(args.peec_proximity_relax),
        )
        I_fil = prox["I_fil"]; V_port = prox["V_port"]
        proximity_info = {
            "n_iter": int(prox["n_iter"]),
            "converged": bool(prox["converged"]),
        }
    else:
        I_fil, V_port = solve_loop_bundle(R_f, L_f, args.frequency,
                                           I_port=args.current,
                                           Zs_fil=Zs_fil)
    t_peec = time.perf_counter() - t0
    Z_coil = V_port / args.current
    L_coil = Z_coil.imag / omega if omega > 0 else 0.0
    R_coil = Z_coil.real
    if use_prox and proximity_info is not None:
        progress("PEEC",
            f"L_coil={L_coil * 1e9:.3f} nH, "
            f"R_coil={R_coil * 1e3:.4f} mΩ "
            f"(proximity {proximity_info['n_iter']} iter, "
            f"{'conv' if proximity_info['converged'] else 'NO-conv'}, "
            f"{t_peec:.1f}s)")
    else:
        progress("PEEC", f"L_coil={L_coil * 1e9:.3f} nH, "
                          f"R_coil={R_coil * 1e3:.4f} mΩ ({t_peec:.1f}s)")

    out = {
        "source_type": "filament",
        "L_coil": L_coil, "R_coil": R_coil,
        "paths": paths, "I_fil": I_fil,
        "n_filaments": len(paths),
        "t_coil_topology_s": t_topo,
        "t_coil_solve_s": t_peec,
    }
    if proximity_info is not None:
        out["proximity"] = proximity_info
    return out


# ======================================================================
# Coil solver: BEM-A (intree Weggler EFIE saddle on RWG)
# ======================================================================
def _build_bema_coil_mesh(args):
    """Load pre-meshed coil surface from --coil-vol.

    BEM-A operates on a Cubit-/Netgen-meshed .vol whose 2D boundary
    elements form the coil surface.  The .vol must already contain
    boundary labels matching args.coil_source_name / args.coil_sink_name
    on the two cap faces (Layer-1 Cubit responsibility); no on-the-fly
    OCC re-mesh and no auto-detection -- fail-fast per CLAUDE.md
    "No Fallbacks" policy if the labels are missing.
    """
    from ngsolve import Mesh

    if not args.coil_vol:
        raise ValueError(
            "--coil-vol is required for --coil-solver bem-a.  PEEC consumes "
            "STEP (--coil-step); BEM-A consumes a pre-meshed surface .vol.")
    if not os.path.isfile(args.coil_vol):
        raise FileNotFoundError(f"--coil-vol not found: {args.coil_vol}")

    mesh = Mesh(args.coil_vol)
    bnd_names = list(mesh.GetBoundaries())
    mat_names = list(mesh.GetMaterials())
    src_name = args.coil_source_name
    snk_name = args.coil_sink_name
    if src_name not in bnd_names or snk_name not in bnd_names:
        # Diagnostic hint for Block-vs-Sideset confusion:
        # in Cubit, source/sink may have been assigned to a BLOCK
        # (volume material) instead of a SIDESET (surface boundary),
        # which puts the names into mesh.GetMaterials() but NOT
        # GetBoundaries().  Surface the materials list too so the
        # user can spot the mistake without re-running anything.
        in_mats = [n for n in (src_name, snk_name) if n in mat_names]
        hint = ""
        if in_mats:
            hint = (f"  HINT: {in_mats!r} appear in the .vol's "
                    f"MATERIAL labels (mesh.GetMaterials() -> "
                    f"{mat_names}) but NOT the boundary labels.  "
                    f"In Cubit, assign source/sink as SIDESETS, not "
                    f"BLOCKS, and re-export.  See CLAUDE.md "
                    f"\"Cubit Block/Sideset Label Convention\".")
        # Case-insensitive near-match hint (e.g. user wrote "Source"
        # in Cubit but the panel default is "source"):
        bnd_lower = {n.lower(): n for n in bnd_names}
        near = [bnd_lower[n.lower()] for n in (src_name, snk_name)
                if n.lower() in bnd_lower and bnd_lower[n.lower()] != n]
        if near:
            hint += (f"  CASE HINT: case-insensitive matches found: "
                      f"{near!r}.  Either rename the Cubit sideset to "
                      f"match {src_name!r}/{snk_name!r}, or pass "
                      f"--coil-source-name / --coil-sink-name with the "
                      f"actual case as it appears in the .vol.")
        # Fuzzy / typo hint (e.g. user wrote "sorce" instead of
        # "source").  Use difflib.get_close_matches at cutoff=0.7 so
        # 1-2 char edits count but unrelated names don't.  This is the
        # actual kubota / mdx 2026-05-11 case (sideset name typo).
        try:
            import difflib
            typos = []
            for required in (src_name, snk_name):
                if required in bnd_names:
                    continue
                matches = difflib.get_close_matches(
                    required, bnd_names, n=3, cutoff=0.7)
                # Skip matches already covered by the case hint above
                # so we don't double-report.
                matches = [m for m in matches if m.lower() != required.lower()]
                if matches:
                    typos.append((required, matches))
            if typos:
                lines = [f"{r!r} ~~ {m!r}" for r, m in typos]
                hint += (f"  TYPO HINT: similar names in the .vol "
                          f"(possible typo): {'; '.join(lines)}.  Fix "
                          f"the sideset name in the Cubit .jou or pass "
                          f"--coil-source-name / --coil-sink-name with "
                          f"the actual spelling.")
        except Exception:
            pass
        raise ValueError(
            f"--coil-vol {args.coil_vol!r} is missing the required boundary "
            f"labels {src_name!r}/{snk_name!r}.  Available boundaries: "
            f"{bnd_names}.{hint}  Set the source/sink sidesets in the "
            f"Cubit .jou before exporting the coil .vol.")

    # compute_inductance_source_sink uses HDivSurface(mesh, order=0) +
    # LaplaceSL — these assume a PURE SURFACE MESH (no volume elements).
    # When the .vol has volume tetrahedra (the common Cubit export
    # case), the saddle-point K matrix gains a 0-pivot from extra null
    # modes that the existing D[:-1, :] deflation doesn't remove, and
    # LU fails with "singular matrix detected".
    # Convert volume → surface in-memory before returning.
    # (Keiko mdx + LAB gapped_torus_2port both reproduced 2026-05-12.)
    if mesh.ne > 0:
        from surface_mesh_extract import _extract_surface_mesh_filtered
        mesh = _extract_surface_mesh_filtered(mesh, keep_label="")
    return mesh


def _arrays_from_bema_coil_mesh(mesh, source_name="source", sink_name="sink"):
    """Extract verts / tris / source-mask / sink-mask arrays.

    The source/sink boundary names default to "source" / "sink" but
    callers MUST pass ``source_name`` / ``sink_name`` from
    ``args.coil_source_name`` / ``args.coil_sink_name`` so the array
    extraction stays in lockstep with the .vol-validation step in
    ``_build_bema_coil_mesh``.  Hardcoding the strings here was a
    silent failure mode: validation passed (using the user's name) but
    the masks were empty (using "source"/"sink"), giving a misleading
    "coil source/sink mesh empty after labelling" runtime error
    (kubota / mdx, 2026-05-11).
    """
    from ngsolve import BND
    n_v = mesh.nv
    verts = np.empty((n_v, 3), dtype=np.float64)
    for i in range(n_v):
        verts[i] = [mesh.vertices[i].point[k] for k in range(3)]
    bnd_names = list(mesh.GetBoundaries())
    src_idx = {i for i, n in enumerate(bnd_names) if n == source_name}
    snk_idx = {i for i, n in enumerate(bnd_names) if n == sink_name}
    tris, src_mask, snk_mask = [], [], []
    for el in mesh.Elements(BND):
        tris.append([int(v.nr) for v in el.vertices])
        src_mask.append(el.index in src_idx)
        snk_mask.append(el.index in snk_idx)
    return (verts, np.array(tris, dtype=np.int64),
            np.array(src_mask, dtype=bool),
            np.array(snk_mask, dtype=bool))


def _solve_coil_bem_a(args):
    """BEM-A coil layer via ngsolve.bem (Weggler EFIE saddle on HDivSurface RT0).

    Replaced the intree pure-Python assembler 2026-05-03 after benchmarks
    showed ngsolve.bem was 50-60x faster across all tested sizes.  See
    ``src/radia/bem/coil_inductance_ngsolve.py`` docstring for the full
    rationale.

    Returns a dict with:
      * L_coil [H] (external inductance), R_coil [Ω] at unit terminal
        current.  R comes from the impedance-EFIE (Z_s inside the
        saddle, J = finite-impedance current) -- the sole formulation
        since 2026-07-02; the PEC post-hoc R was removed (it
        over-estimated ~3x on tightly-wound coils).
      * source_type = "surface"
      * coil_centroids (n_t, 3), coil_areas (n_t,) -- real;
        coil_J_per_tri (n_t, 3) -- COMPLEX per-triangle surface current
        at unit terminal current (Im is all-zero for frequency=0).
        Multiply by complex args.current downstream as needed.
      * coil_residual, coil_mesh_nv, coil_mesh_n_tris
    """
    from radia.bem.coil_inductance_ngsolve import (
        compute_inductance_source_sink, compute_centroids_areas_J,
        check_source_sink_current_path)

    omega = 2.0 * math.pi * args.frequency
    progress("BEMA", f"loading pre-meshed coil .vol: "
                       f"{os.path.basename(args.coil_vol)}")
    t0 = time.perf_counter()
    coil_mesh = _build_bema_coil_mesh(args)
    # We still extract verts/tris/masks for diagnostics + N_tris reporting.
    # Pass user-supplied source/sink names so the array masks stay in
    # sync with the .vol-validation step inside _build_bema_coil_mesh.
    coil_verts, coil_tris, src_mask, snk_mask = \
        _arrays_from_bema_coil_mesh(
            coil_mesh,
            source_name=args.coil_source_name,
            sink_name=args.coil_sink_name)
    t_mesh = time.perf_counter() - t0
    progress("BEMA",
        f"coil nv={coil_mesh.nv} n_tris={len(coil_tris)} "
        f"src/snk={int(src_mask.sum())}/{int(snk_mask.sum())} "
        f"({t_mesh:.1f}s)")
    if int(src_mask.sum()) == 0 or int(snk_mask.sum()) == 0:
        raise RuntimeError(
            f"coil source/sink mesh empty after labelling "
            f"(src={int(src_mask.sum())}, snk={int(snk_mask.sum())}).  "
            f"Check that the {args.coil_source_name!r}/"
            f"{args.coil_sink_name!r} sidesets in --coil-vol contain "
            f"surface elements (not just labelled but empty).")

    # Fail loud BEFORE the expensive EFIE solve when the source/sink labels
    # let the current shortcut across the conductor instead of running
    # around the coil (silently returns a tiny, meaningless L -- see the
    # ih_fem_kelvin 0.28 nH vs 90 nH incident, 2026-07-15).
    path_info = check_source_sink_current_path(
        coil_verts, coil_tris, src_mask, snk_mask,
        source_label=args.coil_source_name,
        sink_label=args.coil_sink_name)
    if path_info:
        progress("BEMA",
            f"source->sink current path {path_info['geodesic_m'] * 1e3:.1f} mm "
            f"({path_info['ratio']:.1f}x coil extent) -- OK")

    # Complex Leontovich surface impedance for the impedance-EFIE
    # (the sole AC formulation).  omega == 0 -> DC vacuum-L solve.
    delta_skin = math.sqrt(2.0 / (omega * MU_0 * args.coil_sigma)) \
                  if omega > 0 else 1.0
    Z_s_coil_complex = (1.0 + 1.0j) / (args.coil_sigma * delta_skin) \
                        if omega > 0 else None

    # Saddle-point solver selection.
    #   "auto" (default):  LU below a tri-count threshold, Krylov above.
    #       The AC impedance-EFIE saddle is COMPLEX (16 B/entry, 2x the
    #       old real PEC saddle the historical 10k threshold was
    #       calibrated on), so the AC LU cutoff is halved to 5000 tris;
    #       the omega=0 real vacuum solve keeps 10000.  Above the
    #       cutoff: cocr (div-free loop reduction + complex-symmetric
    #       COCR, ~24 mesh-independent iters, exact vs LU) for both AC
    #       and DC -- it replaces the historical GMRES-for-AC /
    #       MINRES-for-DC branch (GMRES stalled unpreconditioned on the
    #       indefinite AC saddle at ~1e5 matvecs).
    #   "lu"          dense LAPACK direct (overwrite_a=True, smallest
    #                 memory for direct solve)
    #   "cocr"        div-free loop reduction + COCR, dense SL matvec --
    #                 SCALABLE large-N
    #   "hacapk_cocr" same, HACApKBEMManager-compressed O(N log N) SL
    #                 matvec (accuracy ~3e-7)
    #   "minres"      Krylov for REAL symmetric indefinite -- omega=0 ONLY
    #   "gmres"       Generic complex-capable Krylov -- unpreconditioned,
    #                 stalls on large saddles (kept for comparison)
    saddle = args.coil_saddle_solver
    if saddle == "auto":
        lu_cutoff = 5000 if omega > 0 else 10000
        # Large systems -> COCR (div-free loop reduction + complex-symmetric
        # COCR, ~24 mesh-independent iters, exact vs LU) instead of the O(N^3)
        # LU or the unpreconditioned GMRES that stalled (~1e5 matvecs) on the
        # indefinite AC saddle.  Dense SL matvec; hacapk_cocr is the HACApK
        # opt-in.  Valid for AC and DC.
        saddle = "lu" if len(coil_tris) < lu_cutoff else "cocr"
    # cocr / hacapk_cocr are passed straight through: the library solver
    # values are the same two names (dense vs HACApK-compressed SL matvec).
    progress("BEMA",
        f"ngsolve.bem impedance-EFIE solve (n_tris={len(coil_tris)}, "
        f"fes_order=0, "
        f"Z_s={Z_s_coil_complex if omega > 0 else 'DC (vacuum L)'}, "
        f"saddle={saddle})")
    t0 = time.perf_counter()
    res = compute_inductance_source_sink(
        coil_mesh,
        source_label=args.coil_source_name,
        sink_label=args.coil_sink_name,
        fes_order=0,
        solver=saddle,
        omega=omega,
        Z_s_complex=Z_s_coil_complex,
        log_fn=progress,
    )
    t_solve = time.perf_counter() - t0
    if "error" in res:
        raise RuntimeError(f"BEM-A solve failed: {res['error']}")
    L_coil = float(res["L"])
    R_coil = float(res["R"])
    residual = float(res["residual"])
    # DoF breakdown of the BEM-A saddle system (HDivSurface RT0 J + SurfaceL2
    # divergence multiplier f).  Total saddle size = n_J + n_f - 1 (the -1
    # comes from the Lagrange constraint reduction inside the solver, where
    # one row of the divergence multiplier is dropped to fix the gauge).
    coil_n_J = int(res.get("n_J", 0))
    coil_n_f = int(res.get("n_f", 0))
    coil_saddle_ndof = coil_n_J + coil_n_f - 1 if coil_n_J and coil_n_f else 0
    progress("BEMA",
        f"L_coil={L_coil * 1e9:.3f} nH, R_coil(SIBC)={R_coil * 1e3:.4f} mΩ "
        f"residual={residual:.2e} ({t_solve:.1f}s)")
    progress("BEMA",
        f"coil DoF: n_J={coil_n_J} (HDivSurface RT0), n_f={coil_n_f} "
        f"(SurfaceL2 P0), saddle_ndof={coil_saddle_ndof} "
        f"on {len(coil_tris)} triangles / {coil_mesh.nv} vertices")

    # Per-triangle COMPLEX J_s vectors at centroids (workpiece bridge
    # prep).  The impedance-EFIE current is complex; sample Re and Im
    # GridFunctions separately and combine so the phi_inc / Telegen /
    # qsurf bridges downstream see the full phasor (they already split
    # Re/Im before calling the real-only compute_phi_inc_from_surface_J).
    coil_centroids, coil_areas, coil_J_re = compute_centroids_areas_J(
        coil_mesh, res["gf_J"])
    _, _, coil_J_im = compute_centroids_areas_J(
        coil_mesh, res["gf_J_im"])
    coil_J_per_tri = coil_J_re + 1j * coil_J_im

    return {
        "source_type": "surface",
        "L_coil": L_coil, "R_coil": R_coil,
        "coil_centroids": coil_centroids,
        "coil_areas": coil_areas,
        "coil_J_per_tri": coil_J_per_tri,
        "coil_mesh_nv": int(coil_mesh.nv),
        "coil_mesh_n_tris": int(len(coil_tris)),
        "coil_mesh_n_src_tris": int(src_mask.sum()),
        "coil_mesh_n_snk_tris": int(snk_mask.sum()),
        "coil_n_J": coil_n_J,
        "coil_n_f": coil_n_f,
        "coil_saddle_ndof": coil_saddle_ndof,
        "bem_a_residual": residual,
        "t_coil_topology_s": t_mesh,
        "t_coil_solve_s": t_solve,
    }


# ======================================================================
# Surface-current Biot-Savart helpers (mirror compute_phi_inc_from_surface_J
# but accept complex J -- output complex H -- for Telegen ΔL evaluation
# in the BEM-A coil case)
# ======================================================================
def _H_from_surface_J_complex(eval_pts, src_centroids, src_areas, src_J_complex):
    """Vectorised Biot-Savart H from a complex surface-current source.

    H(r) = (1/4π) Σ_t [J_t × (r - c_t)] / |r - c_t|^3 * A_t
    """
    dx = eval_pts[:, None, :] - src_centroids[None, :, :]
    r2 = np.sum(dx * dx, axis=2)
    r2 = np.maximum(r2, 1e-60)
    r3_inv = src_areas[None, :] / (r2 * np.sqrt(r2))
    cross = np.cross(src_J_complex[None, :, :], dx)
    return INV_4PI * np.sum(cross * r3_inv[..., None], axis=1)


def _delta_L_telegen_phiB_from_surface_J(
        coil_centroids, coil_areas, coil_J_complex,
        wp_centroids, wp_areas, wp_norms, wp_phi_avg, I_port):
    """Gauge-invariant Δ L_port via ∫_S φ (n · B_inc) dS, surface-J source.

    Parallel to ``radia.workpiece_surface.delta_L_telegen_phiB`` but for
    BEM-A surface-current coil.  1-point centroid quadrature for B_inc.
    """
    H_inc = _H_from_surface_J_complex(
        wp_centroids, coil_centroids, coil_areas, coil_J_complex)
    B_inc = MU_0 * H_inc
    n_dot_B = np.sum(wp_norms * B_inc, axis=1)
    integ = np.sum(wp_phi_avg * n_dot_B * wp_areas)
    Ip = complex(I_port)
    if abs(Ip) < 1e-30:
        return 0.0 + 0.0j
    return integ / (Ip * Ip)


# ======================================================================
# Workpiece weak-coupled BEM-SIBC block (shared by both coil solvers)
# ======================================================================
def _wp_genus_check(wp_mesh, tag="BEM"):
    """Euler characteristic / genus of the extracted workpiece surface.

    Returns ``(chi, genus)`` and logs a WARNING for genus >= 1: the scalar
    potential BIE (J_s = n x -grad phi, single-valued phi) carries ZERO net
    current through any cut of the surface, so on a handle that links the
    coil flux the physical shorted-turn eddy current -- and its Lenz
    screening -- is unrepresentable.  The analytic genus-0 sphere locks
    the simply connected path; flux-linked genus-1 heating needs the
    explicit loop-DOF extension.
    """
    from radia.bem_sibc_solver import surface_euler_characteristic
    chi = surface_euler_characteristic(wp_mesh)
    genus = (2 - chi) // 2
    if chi != 2:
        progress(tag,
            f"WARNING: workpiece surface has genus={genus} (Euler chi={chi}, "
            f"e.g. a tube/ring).  The scalar-potential BIE cannot carry a net "
            f"circulating (shorted-turn) eddy current on the handle, so its "
            f"Lenz screening is LOST: H_t / P_wp are over-estimated when the "
            f"coil flux links the handle.  L / delta_L remain usable.  "
            f"Use --wp-loop-dof on the supported weak-coupling path for "
            f"absolute heating.")
    return chi, genus


_GENUS_P_WP_CAVEAT = (
    "workpiece surface genus >= 1 (Euler chi != 2): the scalar-potential "
    "BIE cannot represent the net circulating (shorted-turn) eddy current "
    "on the handle, so its Lenz screening is missing and H_t / P_wp are "
    "over-estimated when the coil flux links the handle.  delta_L is "
    "unaffected.  On the supported weak-coupling path, use --wp-loop-dof "
    "to add the missing cohomology mode; the extension is locked by the "
    "analytic shorted-ring golden.")


def _apply_wp_loop_dof(args, bem, phi_inc, Z_s_wp, omega, coil_data,
                       res_bem, wp_genus, wp_chi):
    """Apply the genus-1 loop-DOF extension to the linear-SIBC weak solve.

    Replaces the reported P_wp / H_t_rms with the loop-extended values
    (``radia.bem_loop_extension.solve_loop_extended``); ``phi_vec`` stays
    the plain solve's, so the Telegen dL and the qsurf spatial pattern
    keep their validated plain-phi convention (the loop DOF corrects the
    DISSIPATION, not the reactive coupling -- cf. the genus caveat
    "delta_L is unaffected").  Returns ``(res_bem_updated, loop_meta)``.
    """
    import cmath

    if wp_genus != 1:
        raise ValueError(
            f"--wp-loop-dof requires a genus-1 workpiece surface "
            f"(got genus={wp_genus}, Euler chi={wp_chi}).  A genus-0 "
            f"workpiece needs no loop DOF (the sphere benchmark validates "
            f"the plain solver to 0.3%); genus >= 2 is not implemented.")
    from radia.bem_loop_extension import (solve_loop_extended,
                                          A_from_filaments)
    if coil_data["source_type"] == "filament":
        def A_inc_fn(points):
            return A_from_filaments(points, coil_data["paths"],
                                    coil_data["I_fil"])
    elif coil_data["source_type"] == "surface":
        from radia.bem_sibc_solver import A_from_surface_J
        coil_Jc = (complex(args.current)
                   * coil_data["coil_J_per_tri"]).astype(complex)

        def A_inc_fn(points):
            A_re = A_from_surface_J(points, coil_data["coil_centroids"],
                                    coil_data["coil_areas"],
                                    np.real(coil_Jc))
            A_im = A_from_surface_J(points, coil_data["coil_centroids"],
                                    coil_data["coil_areas"],
                                    np.imag(coil_Jc))
            return A_re + 1j * A_im
    else:
        raise RuntimeError(
            f"unknown coil source_type {coil_data['source_type']!r}")

    t0 = time.perf_counter()
    loop_out = solve_loop_extended(bem, phi_inc, Z_s_wp, omega, A_inc_fn)
    t_loop = time.perf_counter() - t0

    # sanity: the frozen sub-solve must reproduce the plain solve
    P_plain = float(res_bem["P_density"]) * float(res_bem["area"])
    frz_rel = abs(loop_out["P_frozen"] - P_plain) / max(P_plain, 1e-30)
    if frz_rel > 1e-3:
        raise RuntimeError(
            f"loop-DOF frozen sub-solve disagrees with the plain BIE solve "
            f"by {frz_rel:.2e} (P {loop_out['P_frozen']:.4e} vs "
            f"{P_plain:.4e} W) -- operator mismatch, refusing to report "
            f"loop-extended numbers.")

    alpha = loop_out["alpha"]
    progress("BEM",
        f"loop-DOF: alpha={abs(alpha):.4g} A "
        f"(phase {math.degrees(cmath.phase(alpha)):+.1f} deg), "
        f"P_wp {P_plain:.4e} -> {loop_out['P_total']:.4e} W, "
        f"H_t {res_bem['H_t_rms']:.4g} -> {loop_out['H_t_rms']:.4g} A/m "
        f"({t_loop:.1f}s)")
    res_bem = dict(res_bem)
    res_bem["P_density"] = loop_out["P_total"] / float(res_bem["area"])
    res_bem["H_t_rms"] = loop_out["H_t_rms"]
    loop_meta = {
        "wp_loop_dof": True,
        "wp_loop_alpha_A": float(abs(alpha)),
        "wp_loop_alpha_deg": float(math.degrees(cmath.phase(alpha))),
        "wp_loop_theta_jump": float(loop_out["theta_jump"]),
        "wp_loop_cut_n_vertices": int(loop_out["cut_n_vertices"]),
        "t_loop_dof_s": float(t_loop),
    }
    return res_bem, loop_meta


# Fail-loud bound on the surface-Poisson tangential-gradient residual
# ||grad_S psi + H_t,inc|| / ||H_t,inc||.  A coil current piercing the
# workpiece surface (psi does not exist) or a broken incident-H evaluation
# produces an O(1) residual; 10% is the fail-loud contract.
_PHI_INC_POISSON_MAX_RESIDUAL = 0.10


def _phi_inc_poisson(args, coil_data, wp_mesh, obs):
    """phi_inc via surface-Poisson projection of the exact vertex H_inc.

    Evaluates the incident H at the workpiece vertices with the coil
    source's EXACT kernel (finite-segment Biot-Savart for filaments,
    centroid-quadrature panel Biot-Savart for the BEM-A surface J --
    the same helpers the Telegen bridge uses), then reconstructs the
    mean-zero psi with ``radia.bem_sibc_solver.
    compute_phi_inc_surface_poisson`` (fails loud when H_t,inc is not a
    surface gradient).  Wall-free replacement for the axis-ray +
    horizontal-ray path integration, selected by ``--wp-phi-inc
    poisson``.  Returns ``(psi, grad_residual)``.
    """
    from ngsolve import BND
    from radia.bem_sibc_solver import compute_phi_inc_surface_poisson

    if len(obs) != wp_mesh.nv:
        raise RuntimeError(
            f"--wp-phi-inc poisson expects the P1 vertex-nodal DOF path "
            f"(len(obs)={len(obs)} != wp_mesh.nv={wp_mesh.nv}).")

    if coil_data["source_type"] == "filament":
        from radia.biot_savart import h_segments_batch
        H = np.zeros((len(obs), 3), dtype=complex)
        for fil_segs, Ik in zip(coil_data["paths"], coil_data["I_fil"]):
            H += complex(Ik) * h_segments_batch(fil_segs, obs)
    elif coil_data["source_type"] == "surface":
        coil_Jc = (complex(args.current)
                   * coil_data["coil_J_per_tri"]).astype(complex)
        H = _H_from_surface_J_complex(
            obs, coil_data["coil_centroids"], coil_data["coil_areas"],
            coil_Jc)
    else:
        raise RuntimeError(
            f"unknown coil source_type {coil_data['source_type']!r}")

    tris = np.array([[v.nr for v in el.vertices]
                     for el in wp_mesh.Elements(BND)], dtype=np.int64)
    return compute_phi_inc_surface_poisson(
        obs, tris, H, max_grad_residual=_PHI_INC_POISSON_MAX_RESIDUAL)


def _solve_workpiece_weak_coupled(args, coil_data):
    """Workpiece BEM-SIBC + Telegen ΔL using whichever coil source exists.

    Branches on coil_data["source_type"] for both phi_inc bridge and
    Telegen evaluation:
      * "filament" -> compute_phi_inc_from_filaments + delta_L_telegen_phiB
      * "surface"  -> compute_phi_inc_from_surface_J + _delta_L_telegen_phiB_from_surface_J
    """
    from ngsolve import (Mesh, BND, Integrate, CF, GridFunction,
                          InnerProduct, grad)
    from ngsolve import (specialcf as _scf, BND as _BND,
                          Integrate as _Int, CF as _CF, GridFunction as _GF)
    from surface_mesh_extract import _extract_surface_mesh_filtered
    from radia.bem_sibc_solver import (
        ScalarBIESIBCSolver, compute_phi_inc_from_filaments,
        compute_phi_inc_from_surface_J)
    from em_material import EMMaterial
    # Hole-extractor for wp-as-hole geometry (closed-torus convention).
    # Inlined from the deleted calc_peec_bem.py.
    _extract_bnd_only = _extract_bnd_only_inline

    omega = 2.0 * math.pi * args.frequency

    # 1. Workpiece mesh
    progress("BEM", f"load wp surface from {os.path.basename(args.vol)}")
    t0 = time.perf_counter()
    vol_mesh = Mesh(args.vol)
    mats = set(vol_mesh.GetMaterials())
    bnds = set(vol_mesh.GetBoundaries())

    # Detect the .vol's baked-in curving order from companion .vol.json
    # (written by Cubit ``export netgen``).  The .vol's
    # ``curvedelements`` section is loaded automatically by Mesh(); a
    # post-load ``mesh.Curve(p)`` would try to UPGRADE and, without a
    # CAD callback, silently fall back to flat -- so we never call it.
    # The curve order is FIXED at Cubit-export time and is not a panel
    # knob.  Capped at 2 for the Lagrange-P2 in-tree assembler (the C++
    # Tri6 geometry kernel takes 6 nodes / tri, so a curve_order=3 .vol
    # is treated as Tri6 here -- the extra mid-face / interior nodes
    # are not currently used by the SS BIE assembler).
    vol_curve_order = min(_detect_vol_curving_order(args.vol), 2)
    # Basis order is a separate, user-selectable knob (--h1-order), not
    # tied to the geometry order.  P1 + P2 supported.
    basis_order = int(args.h1_order)
    if basis_order not in (1, 2):
        raise ValueError(
            f"--h1-order must be 1 or 2 (got {basis_order}).  P3+ "
            f"basis is not implemented in the in-tree Lagrange "
            f"assembler.")
    progress("BEM",
        f"vol curve_order={vol_curve_order} "
        f"(from {os.path.basename(args.vol)}.json), "
        f"basis_order={basis_order} (--h1-order)")
    # The Lagrange-P2 path (basis_order=2) builds M, K, SL, DL all in the
    # in-tree P2 basis directly off the parent vol_mesh's curved Tri6
    # nodes; it requires a sideset bnd_label so it can filter the parent
    # mesh's BND.
    if basis_order >= 2 and args.wp_label not in bnds:
        raise ValueError(
            f"--h1-order=2 requires --wp-label to name a boundary "
            f"sideset (got {args.wp_label!r}, available BND labels: "
            f"{sorted(bnds)}).  Material-adjacency wp_label is only "
            f"supported with basis_order=1 (flat-FES path).")

    if args.wp_label in mats:
        progress("BEM", f"wp_label {args.wp_label!r} as material")
        wp_mesh, _wp_new_to_old = _extract_surface_mesh_filtered(
            vol_mesh, keep_label=args.wp_label, return_vertex_map=True)
    elif args.wp_label in bnds:
        progress("BEM", f"wp_label {args.wp_label!r} as boundary sideset")
        wp_mesh = _extract_bnd_only(vol_mesh, args.wp_label)
        _wp_new_to_old = None
    else:
        raise ValueError(
            f"wp_label {args.wp_label!r} not found in materials "
            f"({sorted(mats)}) or boundaries ({sorted(bnds)})")
    t_wp_mesh = time.perf_counter() - t0
    progress("BEM",
        f"wp nv={wp_mesh.nv} ne(BND)={wp_mesh.GetNE(BND)} ({t_wp_mesh:.1f}s)")
    wp_chi, wp_genus = _wp_genus_check(wp_mesh)

    # Resolve the loop-DOF mode (see --wp-loop-dof).  "on" was already
    # early-guarded in run_inductance and _apply_wp_loop_dof enforces
    # genus-1; "auto" applies exactly when the supported mathematical and
    # solver prerequisites hold and otherwise records why it skipped.
    wp_loop_req = getattr(args, "wp_loop_dof", "auto")
    wp_loop_skip = None
    if wp_loop_req == "on":
        wp_loop_apply = True
    elif wp_genus != 1:
        wp_loop_apply = False
        wp_loop_skip = (f"genus-{wp_genus} surface (the loop DOF applies "
                        f"to genus-1 only; genus-0 needs none)")
    elif args.impedance_model != "sibc":
        wp_loop_apply = False
        wp_loop_skip = ("ESIM impedance model (the Karl loop is not "
                        "integrated with the loop DOF; pass "
                        "--impedance-model sibc for the loop-extended "
                        "solve)")
    elif args.wp_bem_backend != "intree-dense":
        wp_loop_apply = False
        wp_loop_skip = ("the HACApK backend exposes no dense SL/DL for "
                        "the loop column (pass --wp-bem-backend "
                        "intree-dense)")
    elif basis_order != 1:
        wp_loop_apply = False
        wp_loop_skip = "P1 nodal path only (pass --h1-order 1)"
    else:
        wp_loop_apply = True
    if wp_loop_req == "auto" and wp_loop_skip is not None and wp_genus >= 1:
        progress("BEM",
            f"loop-DOF auto: SKIPPED on a genus-{wp_genus} workpiece -- "
            f"{wp_loop_skip}.  The genus P_wp_caveat applies.")

    # 2. ESIM prerequisite check (Karl iteration needs a BH curve).
    if args.impedance_model == "esim" and not args.bh_file:
        raise ValueError(
            "--impedance-model esim requires --bh-file with a 2-column "
            "BH curve (H[A/m], B[T]).  Linear SIBC is the default.")

    # 3. ScalarBIESIBCSolver assembly.
    #
    # Two independent knobs feed the solver:
    #   * vol_curve_order : geometry order (auto-detected from .vol.json).
    #     Forwarded as ``intree_geom_order``.
    #   * basis_order     : Lagrange basis order (user-selectable, --h1-order).
    #     Forwarded as ``order``.
    #
    # When basis_order=2 the in-tree Lagrange-P2 path bypasses NGSolve's
    # H1 hierarchical FES entirely and routes through the parent
    # ``vol_mesh`` (whose ``curvedelements`` is already loaded) plus the
    # ``bnd_label`` filter; the extracted flat ``wp_mesh`` is used only
    # for post-processing (Telegen per-element ΔL averaging).
    # When basis_order=1 we keep the existing flat-extracted wp_mesh
    # path; H1+ngsolve.bem trace integration runs on flat triangles
    # consistent with the P1 hat (curving has no benefit at P1 since the
    # basis cannot resolve geometric curvature).
    use_curved_p2 = (basis_order >= 2)
    bem_input_mesh = vol_mesh if use_curved_p2 else wp_mesh
    bem_bnd_label = args.wp_label if use_curved_p2 else None
    progress("BEM",
        f"assembly (basis_order={basis_order}, "
        f"vol_curve_order={vol_curve_order}, "
        f"backend={args.wp_bem_backend}, "
        f"input_mesh={'vol_mesh+bnd_label' if use_curved_p2 else 'wp_mesh'})")
    t0 = time.perf_counter()
    if args.wp_bem_backend == "intree-dense":
        bem = ScalarBIESIBCSolver(
            bem_input_mesh, order=basis_order, assemble_dense=True,
            use_intree_bem=True, intree_geom_order=vol_curve_order,
            intree_singular_n_q=6, intree_regular_quad_degree=7,
            bnd_label=bem_bnd_label, log_fn=progress)
    elif args.wp_bem_backend == "hacapk":
        bem = ScalarBIESIBCSolver(
            bem_input_mesh, order=basis_order, assemble_dense=True,
            use_intree_bem=True, intree_geom_order=vol_curve_order,
            intree_singular_n_q=6, intree_regular_quad_degree=7,
            use_intree_hacapk=True,
            hacapk_aca_eps=args.wp_aca_eps,
            hacapk_leaf=64, hacapk_eta=2.0,
            bnd_label=bem_bnd_label, log_fn=progress)
    else:
        raise ValueError(
            f"unknown --wp-bem-backend {args.wp_bem_backend!r}; "
            f"supported: 'hacapk' (production), 'intree-dense' (debug).")
    t_asm = time.perf_counter() - t0
    # wp DoF context: at P1 the in-tree H1 hat assembler uses one DoF per
    # vertex of the EXTRACTED wp_mesh, at P2 the in-tree Lagrange-P2 path
    # consumes the parent vol_mesh's curved Tri6 nodes (vertex + edge DoFs
    # on the wp_label sideset).  Printing the triangle count alongside
    # ndof makes the basis-vs-mesh density relationship visible in the
    # log (no need to grep both upstream "wp nv=..." and "ndof=...").
    wp_ne_bnd = int(wp_mesh.GetNE(BND))
    progress("BEM",
        f"ndof={bem.ndof} (P{basis_order} on {wp_ne_bnd} BND triangles, "
        f"{wp_mesh.nv} wp_mesh vertices) ({t_asm:.1f}s)")

    # Consolidated DoF summary -- one line that shows the full problem
    # size so users grepping logs after the fact can see coil + wp +
    # total in one place (instead of stitching together upstream
    # "BEMA" lines from _solve_coil_bem_a + this "BEM" line).
    coil_src = coil_data.get("source_type", "?")
    if coil_src == "filament":
        coil_dof_str = (f"PEEC filament n={coil_data.get('n_filaments', 0)}")
        total_dof = int(coil_data.get("n_filaments", 0)) + int(bem.ndof)
    elif coil_src == "surface":
        c_saddle = int(coil_data.get("coil_saddle_ndof", 0))
        c_tris = int(coil_data.get("coil_mesh_n_tris", 0))
        coil_dof_str = f"BEM-A saddle_ndof={c_saddle} on {c_tris} tri"
        total_dof = c_saddle + int(bem.ndof)
    else:
        coil_dof_str = f"{coil_src}"
        total_dof = int(bem.ndof)
    progress("BEM",
        f"SUMMARY: coil[{coil_dof_str}] + wp[ndof={bem.ndof} "
        f"P{basis_order} {wp_ne_bnd}tri] => total_DoF={total_dof}")

    # 4. phi_inc on workpiece — branch on coil source type
    if getattr(bem, "_intree_lagrange_p2", False):
        obs = np.ascontiguousarray(bem.dof_coords)
    else:
        obs = np.array([[wp_mesh.vertices[i].point[j] for j in range(3)]
                         for i in range(wp_mesh.nv)])
    # Incident-potential reconstruction is basis-determined, NOT a knob:
    # P1 -> surface-Poisson psi from the exact vertex H_inc (wall-free,
    # fails loud if H_t,inc is not a surface gradient); Lagrange-P2 ->
    # path integration (the sole implementation for edge-node DOFs; it
    # carries the branch-cut-wall limitation documented on
    # compute_phi_inc_from_filaments).  The legacy P1 path-integration
    # route was REMOVED 2026-07-17 once the surface-Poisson
    # reconstruction was covered by analytic goldens and integration
    # contracts: a superseded twin with a known wall bias must not stay
    # selectable.
    phi_inc_mode = "poisson" if basis_order == 1 else "path"
    progress("BEM", f"phi_inc from coil ({coil_data['source_type']}, "
                    f"mode={phi_inc_mode})")
    t0 = time.perf_counter()
    phi_inc_grad_residual = None
    if phi_inc_mode == "poisson":
        # Surface-Poisson reconstruction from the exact vertex H_inc
        # (wall-free; fails loud when H_t,inc is not a surface
        # gradient).  Early guards in run_inductance restrict this to
        # the weak-coupled P1 nodal path.
        phi_inc, phi_inc_grad_residual = _phi_inc_poisson(
            args, coil_data, wp_mesh, obs)
        progress("BEM",
            f"surface-Poisson psi: grad-consistency residual "
            f"{phi_inc_grad_residual:.1%} "
            f"(gate {_PHI_INC_POISSON_MAX_RESIDUAL:.0%})")
    elif coil_data["source_type"] == "filament":
        # Reverted 2026-05-21: surface-path phi_inc gave P_wp = 1.92 W
        # which is below the weak-coupling lower bound (5.26 W) -- the
        # different gauge convention (phi(ref_vertex on wp) = 0 vs
        # axis-ray's phi(infty) = 0 with dipole tail) altered the BIE
        # solve in a non-physical way.  The actual spatial-peak fix
        # was in the per-vertex |H_t|^2 reconstruction (triangle-wise
        # gradient instead of Galerkin phi*K@phi localization).
        phi_inc = compute_phi_inc_from_filaments(
            obs, coil_data["paths"], coil_data["I_fil"])
    elif coil_data["source_type"] == "surface":
        # Multiply the COMPLEX per-triangle surface J (impedance-EFIE
        # phasor at unit current) by the terminal current; Re/Im are
        # then bridged separately (compute_phi_inc_from_surface_J is
        # real-only and raises on complex input).
        coil_J_complex = (complex(args.current)
                           * coil_data["coil_J_per_tri"]).astype(complex)
        phi_re = compute_phi_inc_from_surface_J(
            obs, coil_data["coil_centroids"], coil_data["coil_areas"],
            np.real(coil_J_complex))
        phi_im = compute_phi_inc_from_surface_J(
            obs, coil_data["coil_centroids"], coil_data["coil_areas"],
            np.imag(coil_J_complex))
        phi_inc = phi_re + 1j * phi_im
    else:
        raise RuntimeError(
            f"unknown coil source_type {coil_data['source_type']!r}")
    t_phi = time.perf_counter() - t0
    progress("BEM", f"phi_inc ({t_phi:.1f}s)")

    # 5. Workpiece SIBC: BIE solve (optionally Karl-iterated for ESIM)
    bh_curve = None
    if args.bh_file:
        from em_material import load_bh_file
        bh_curve = load_bh_file(args.bh_file)
    mat_wp = EMMaterial(name="wp", sigma=args.sigma, mu_r=args.mu_r,
                         bh_curve=bh_curve)
    delta_wp = mat_wp.skin_depth(args.frequency)

    if args.impedance_model == "esim":
        # ESIM cell solver: 1D nonlinear B-H(H) Karl iteration on a
        # cylindrical workpiece slice of radius `half_thickness`.
        # The seed solve at a small H0=5 A/m only needs a rough Z_s
        # to start the outer Karl loop; cap inner Picard at 5 iter
        # to avoid wasting time on a 50-iter solve whose result the
        # next outer iter overrides anyway.
        esim_solver = mat_wp.create_esim_solver(
            args.frequency, args.half_thickness, geometry='cylinder')
        Z_s_seed = complex(esim_solver.solve(5.0, max_iter=5)['Z'])
        if args.esim_per_panel:
            # Per-DOF Z_s works on BOTH intree-dense and HACApK backends
            # (v4.47.2+).  Seed every DOF with the same scalar; subsequent
            # iters refresh Z_s_wp[i] from the per-DOF H_t.
            Z_s_wp = np.full(bem.ndof, Z_s_seed, dtype=complex)
        else:
            Z_s_wp = Z_s_seed
        max_iter = max(int(args.esim_max_iter), 1)
        # Optional Anderson acceleration on the outer Karl loop.
        # m = 0 falls back to plain damped Picard.
        from esim_anderson import AndersonAccelerator
        anderson = AndersonAccelerator(m=int(args.esim_anderson_m),
                                        alpha=float(args.esim_relax))
    else:
        # Linear SIBC (standard Leontovich): Z_s = (1+j) * rho/delta.
        # delta_wp = mat_wp.skin_depth ALREADY includes mu_r
        # (= sqrt(2*rho/(omega*mu0*mu_r))), so rho/delta_wp already equals
        # sqrt(omega*mu0*mu_r/(2*sigma)), i.e. the correct Z_s proportional to
        # sqrt(mu_r) -- matching analytical_formulas.planar_surface_impedance
        # and COMSOL.  A prior `* math.sqrt(args.mu_r)` here double-counted
        # mu_r (Z_s came out proportional to mu_r, sqrt(mu_r)x too large for
        # magnetic workpieces); removed 2026-05-31 after keiko's COMSOL
        # cross-check.  Cu/Al (mu_r=1) were unaffected.
        rho_wp = 1.0 / args.sigma
        Z_s_wp = (1.0 + 1j) * rho_wp / delta_wp
        esim_solver = None
        max_iter = 1

    if isinstance(Z_s_wp, np.ndarray):
        progress("BEM",
            f"BIE solve Z_s=ndarray[{Z_s_wp.shape[0]}] "
            f"(mean|Z_s|={abs(np.mean(Z_s_wp)):.3e}, model=esim per-panel)")
    else:
        progress("BEM",
            f"BIE solve Z_s={Z_s_wp:.3e} (model={args.impedance_model})")
    res_bem = None
    loop_meta = None
    esim_history = []
    esim_converged = (esim_solver is None)
    H_t_rms_iter = 0.0
    t_bie = 0.0
    n_iter_done = 0
    dZ = float("inf")
    for iteration in range(max_iter):
        n_iter_done = iteration + 1
        t0 = time.perf_counter()
        if args.wp_bem_backend == "intree-dense":
            res_bem = bem.solve(phi_inc, Z_s=Z_s_wp, omega=omega)
        else:
            res_bem = bem.solve_hacapk(
                phi_inc, Z_s=Z_s_wp, omega=omega,
                tol=args.wp_gmres_tol, maxiter=500, restart=80)
        t_iter = time.perf_counter() - t0
        t_bie += t_iter
        H_t_rms_iter = float(res_bem["H_t_rms"])

        if esim_solver is None:
            progress("BEM", f"BIE ({t_iter:.1f}s)")
            if wp_loop_apply:
                res_bem, loop_meta = _apply_wp_loop_dof(
                    args, bem, phi_inc, Z_s_wp, omega, coil_data,
                    res_bem, wp_genus, wp_chi)
            break

        # Karl update.
        Z_s_old = Z_s_wp
        if args.esim_per_panel:
            # Extract per-DOF |H_t| via the manual triangle-wise P1
            # gradient of phi (v4.67.0+ formula, matching the q_surf
            # spatial output).  The legacy Galerkin localization
            # ``phi_i * (K @ phi)_i`` is a Laplacian sample, not a
            # gradient-norm sample, and feeds wrong local |H_t| values
            # into the cell solver -- inflating the per-element vs
            # scalar gap by mis-placing the saturation hot-spot.
            from radia.bem_sibc_solver import extract_H_t_per_dof_grad
            phi_vec = np.asarray(res_bem["phi_vec"])
            H_t_per = extract_H_t_per_dof_grad(phi_vec, wp_mesh)
            # Per-DOF ESIM
            Z_s_new = np.empty(bem.ndof, dtype=complex)
            for i in range(bem.ndof):
                Z_s_new[i] = complex(
                    esim_solver.solve(max(float(H_t_per[i]), 1e-3),
                                       max_iter=20)['Z'])
            Z_s_wp = anderson.step(Z_s_old, Z_s_new)
            dZ_per_dof = (np.abs(Z_s_wp - Z_s_old)
                          / np.maximum(np.abs(Z_s_old), 1e-30))
            dZ = float(np.max(dZ_per_dof))
            Zabs_mean = float(np.mean(np.abs(Z_s_wp)))
            Ht_mean = float(np.mean(H_t_per))
            Ht_max = float(np.max(H_t_per))
            esim_history.append({
                "iteration": iteration,
                "Z_s_abs_mean": Zabs_mean,
                "Z_s_abs_max": float(np.max(np.abs(Z_s_wp))),
                "Z_s_abs_min": float(np.min(np.abs(Z_s_wp))),
                "H_t_per_dof_mean": Ht_mean,
                "H_t_per_dof_max": Ht_max,
                "dZ_max": dZ,
                "t_solve": float(t_iter),
            })
            progress("BEM",
                f"ESIM:ITER {iteration} per-panel "
                f"<|Z_s|>={Zabs_mean:.4e} "
                f"<H_t>={Ht_mean:.2f} max(H_t)={Ht_max:.2f} "
                f"max(dZ)={dZ:.4e} t={t_iter:.1f}s")
        else:
            sol_new = esim_solver.solve(max(H_t_rms_iter, 1e-3))
            Z_s_new = complex(sol_new['Z'])
            Z_s_wp = anderson.step(Z_s_old, Z_s_new)
            dZ = abs(Z_s_wp - Z_s_old) / max(abs(Z_s_old), 1e-30)
            esim_history.append({
                "iteration": iteration,
                "Z_s_abs": float(abs(Z_s_wp)),
                "H_t_rms": H_t_rms_iter,
                "dZ": float(dZ),
                "t_solve": float(t_iter),
            })
            progress("BEM",
                f"ESIM:ITER {iteration} |Z_s|={abs(Z_s_wp):.4e} "
                f"H_t={H_t_rms_iter:.2f} dZ={dZ:.4e} t={t_iter:.1f}s")
        # Require iteration > 0 to avoid spurious convergence on the
        # seed Z_s (= esim.solve(5.0)), EXCEPT when max_iter <= 1:
        # the user explicitly asked for one iteration so the dZ we
        # just computed is the only result available.
        if dZ < args.esim_tol and (iteration > 0 or max_iter <= 1):
            progress("BEM", f"ESIM:CONVERGED iter={iteration}")
            esim_converged = True
            break
    else:
        # Loop exhausted max_iter without convergence (ESIM only).
        if esim_solver is not None:
            esim_converged = False
            progress("BEM",
                f"ESIM:NOT-CONVERGED after {max_iter} iter "
                f"(dZ={dZ:.4e} > tol={args.esim_tol:.1e})")
    if esim_solver is not None:
        # Re-solve once at the final Z_s so res_bem reflects the
        # converged impedance (the last iter computed H_t_rms with the
        # PREVIOUS Z_s, not the updated one).
        t0 = time.perf_counter()
        if args.wp_bem_backend == "intree-dense":
            res_bem = bem.solve(phi_inc, Z_s=Z_s_wp, omega=omega)
        else:
            res_bem = bem.solve_hacapk(
                phi_inc, Z_s=Z_s_wp, omega=omega,
                tol=args.wp_gmres_tol, maxiter=500, restart=80)
        t_bie += time.perf_counter() - t0
        # Final per-DOF |H_t| via triangle-gradient (matches Karl loop
        # interior; see extract_H_t_per_dof_grad in bem_sibc_solver.py).
        # Used for the Z_s vs H_t scatter / spatial map in publications.
        if isinstance(Z_s_wp, np.ndarray):
            from radia.bem_sibc_solver import extract_H_t_per_dof_grad
            phi_vec_final = np.asarray(res_bem["phi_vec"])
            H_t_per_final = extract_H_t_per_dof_grad(phi_vec_final, wp_mesh)
        else:
            H_t_per_final = None
    info = res_bem.get("gmres_info", 0)
    progress("BEM",
        f"BIE total ({t_bie:.1f}s over {n_iter_done} iter, gmres info={info})")

    # 5b. Per-DOF |H_t|^2 and q_surf for the spatial qsurf.sol output.
    # This is the data calc_heat.py needs as a Neumann BC for thermal
    # analysis.  Linear AND ESIM paths both need it; H_t_per_final
    # above was ESIM-only.  Same recipe (K @ phi / M_lump) applied
    # unconditionally.  Only meaningful when bem.fes is not None
    # (i.e. basis_order=1 -- H1 P1 on wp_mesh BND); for basis_order>=2
    # phi_vec lives on the parent vol_mesh and the projection back to
    # a wp-surface H1 GridFunction would need explicit coord matching,
    # which the Telegen post-proc below already exercises.  Out of
    # scope for now; warn instead.
    gf_q_wp = None
    q_per_dof = None
    # Spatial q_surf via BIE-calibrated direct Biot-Savart.
    #
    # Rationale (kubota 2026-05-21 follow-up):
    #   - phi_inc via the legacy axis-ray + horizontal-ray path
    #     integral has numerically-induced local-gradient errors
    #     for real coils (leads + multi-loop topology); the BIE
    #     solution phi_vec inherits those errors so the spatial
    #     q distribution comes out scrambled (peak on coil-FAR
    #     face instead of coil-NEAR).
    #   - phi_inc via surface-edge LSQ has correct local gradient
    #     but a different gauge from the BIE's training data, so
    #     the BIE produces an artificially low P_wp (1.9 W vs the
    #     ~6 W axis-ray result on the same geometry).
    #   - Gauge correction (constant shift on phi_inc) does NOT
    #     change BIE energy because the Lagrange multiplier
    #     constraint absorbs constants.
    #
    # Resolution: separate concerns.  Use the BIE's global P_wp
    # (= integrated power, well-trusted) as the magnitude and the
    # direct Biot-Savart |H_t_inc|^2 as the spatial pattern.  This
    # is the weak-coupling approximation rescaled to match the
    # BIE integral.  For typical IH (workpiece inside coil bore,
    # near-uniform shielding factor), |H_t|^2 spatial pattern is
    # very close to |H_t_inc|^2 times a constant -- which is
    # exactly the rescaling we apply.
    #
    # P1 vs P2 (--h1-order): this path is INDEPENDENT of bem.fes
    # because phi_vec is not used.  We only need wp_mesh (for
    # vertex positions + element topology) and res_bem["P_density"]
    # (scalar from the BIE solve).  Both are available regardless
    # of basis_order.
    if True:
        from ngsolve import BND as _BND_imp, H1 as _H1_imp
        n_v = wp_mesh.nv

        # Collect triangle vertex indices + coords.
        tri_v = []
        for el in wp_mesh.Elements(_BND_imp):
            tri_v.append([v.nr for v in el.vertices])
        tri_v = np.asarray(tri_v, dtype=np.int64)
        vert_xyz = np.asarray(
            [list(wp_mesh.vertices[i].point) for i in range(n_v)],
            dtype=float)

        # Triangle areas + outward normals (area-weighted avg to verts).
        p0 = vert_xyz[tri_v[:, 0]]
        p1 = vert_xyz[tri_v[:, 1]]
        p2 = vert_xyz[tri_v[:, 2]]
        nrm = np.cross(p1 - p0, p2 - p0)
        area2 = np.linalg.norm(nrm, axis=1)
        tri_area = 0.5 * area2
        n_hat_tri = nrm / np.maximum(area2[:, None], 1e-30)

        vert_n = np.zeros_like(vert_xyz)
        vert_a = np.zeros(n_v)
        for j in range(3):
            np.add.at(vert_n, tri_v[:, j], tri_area[:, None] * n_hat_tri)
            np.add.at(vert_a, tri_v[:, j], tri_area)
        vert_n /= np.maximum(vert_a[:, None], 1e-30)
        vert_n /= np.maximum(np.linalg.norm(vert_n, axis=1,
                                              keepdims=True), 1e-30)

        # Direct Biot-Savart H at each wp vertex from the coil
        # (curl-free in air, no topology assumption).  Branch on the
        # coil representation, mirroring the phi_inc bridge above:
        # PEEC returns line filaments, BEM-A returns a surface current.
        if coil_data["source_type"] == "filament":
            from radia.bem_sibc_solver import _h_segments_complex as _hbs
            flat_segs, flat_I = [], []
            for fil_segs, Ik in zip(coil_data["paths"], coil_data["I_fil"]):
                Ik_c = complex(Ik)
                for (q1, q2) in fil_segs:
                    flat_segs.append((q1, q2))
                    flat_I.append(Ik_c)
            bs_segs = np.asarray(flat_segs, dtype=float)
            bs_I = np.asarray(flat_I, dtype=complex)
            H_BS = _hbs(bs_segs, vert_xyz, bs_I)        # (n_v, 3) complex
        else:  # "surface" -- BEM-A coil
            coil_J_complex = (complex(args.current)
                              * coil_data["coil_J_per_tri"]).astype(complex)
            H_BS = _H_from_surface_J_complex(
                vert_xyz, coil_data["coil_centroids"],
                coil_data["coil_areas"], coil_J_complex)   # (n_v, 3) complex

        # Tangential H at each vertex: H - (H·n)n.
        Hn = np.sum(H_BS * vert_n, axis=1)
        H_t_vec = H_BS - Hn[:, None] * vert_n
        H_t_sq = (np.abs(H_t_vec[:, 0])**2
                   + np.abs(H_t_vec[:, 1])**2
                   + np.abs(H_t_vec[:, 2])**2)

        # Integrate |H_t|^2 over the wp surface (triangle-mean * area).
        # ∫_S |H_t|^2 dS = sum_tri area_tri * mean_3v(|H_t|^2).
        Hsq_at_verts = H_t_sq[tri_v]                # (n_tri, 3)
        I_norm = float(np.sum(tri_area * Hsq_at_verts.mean(axis=1)))
        # P_wp is the BIE-trusted total -- compute it inline here since
        # the canonical P_wp = res_bem["P_density"] * A_wp is set
        # AFTER this block in the function flow.
        A_wp_local = float(np.sum(tri_area))
        P_wp_local = float(res_bem["P_density"]) * A_wp_local
        if I_norm <= 0 or P_wp_local <= 0:
            progress("BEM",
                f"qsurf.sol skipped: I_norm={I_norm:.3e} "
                f"P_wp={P_wp_local:.3e}")
        else:
            # q_per_dof = (P_wp / ∫|H_t|^2 dS) * |H_t|^2.  Guarantees
            # ∫ q dS = P_wp exactly (BIE energy preserved) AND the
            # spatial pattern follows the topologically-correct
            # Biot-Savart distribution.
            scale = P_wp_local / I_norm
            q_per_dof = scale * H_t_sq
            # Build the wp_mesh-side GridFunction.  For basis_order=1
            # this matches bem.fes; for basis_order>=2 (where bem.fes
            # is None and the BIE uses Lagrange P2 on parent vol_mesh)
            # we build a fresh P1 FES on wp_mesh -- the downstream
            # vol_mesh transfer (via _wp_new_to_old / coord match)
            # is by vertex anyway, independent of basis_order.
            fes_q_local = _H1_imp(wp_mesh, order=1)
            gf_q_wp = _GF(fes_q_local)
            gf_q_wp.vec.FV().NumPy()[:] = q_per_dof
            progress("BEM",
                f"qsurf: BIE-calibrated BS  scale={scale:.3e} "
                f"(P_wp/∫|H_t|²)  ∫|H_t|²dS={I_norm:.3e}")

    # 6. Workpiece scalar outputs
    A_wp = float(Integrate(CF(1), wp_mesh, VOL_or_BND=BND).real)
    P_wp = float(res_bem["P_density"] * A_wp)
    H_t_rms = float(res_bem["H_t_rms"])
    wp_dissipation_R_mOhm = 2.0 * P_wp / (args.current ** 2) * 1e3

    # 7. Telegen ΔL via φ·B form -- aggregate per-triangle wp data first.
    # Two post-proc paths:
    #   * P1 (basis_order=1): phi_vec is indexed by wp_mesh vertices;
    #     bem.fes is the wp_mesh H1 FES; use NGSolve element_wise integrals.
    #   * P2 (basis_order>=2): phi_vec is indexed by Lagrange P2 DOFs on
    #     the parent vol_mesh BND; bem.fes is None.  Build phi-at-wp-
    #     vertices by matching wp_mesh.vertices coords to solver vertex
    #     DOFs (the first n_v entries of dof_coords), then average over
    #     each flat wp_mesh triangle (corner-mean = area-weighted P1 mean).
    n_cf = _scf.normal(3)
    elem_A_int = _Int(_CF(1), wp_mesh, VOL_or_BND=_BND, element_wise=True)
    n_int = [_Int(n_cf[i], wp_mesh, VOL_or_BND=_BND, element_wise=True)
              for i in range(3)]

    if basis_order >= 2:
        # Build phi-at-wp-vertices via coordinate match.
        n_v_solver = len(bem._intree_v_global)
        solver_vert_coords = np.asarray(bem.dof_coords[:n_v_solver])
        coord_to_dof = {}
        for i in range(n_v_solver):
            key = (round(solver_vert_coords[i, 0], 12),
                   round(solver_vert_coords[i, 1], 12),
                   round(solver_vert_coords[i, 2], 12))
            coord_to_dof[key] = i
        phi_at_wp = np.zeros(wp_mesh.nv, dtype=complex)
        n_missed = 0
        for i in range(wp_mesh.nv):
            p = wp_mesh.vertices[i].point
            key = (round(p[0], 12), round(p[1], 12), round(p[2], 12))
            j = coord_to_dof.get(key)
            if j is None:
                n_missed += 1
                continue
            phi_at_wp[i] = res_bem["phi_vec"][j]
        if n_missed > 0:
            raise RuntimeError(
                f"P2 post-proc: {n_missed}/{wp_mesh.nv} wp vertices did "
                f"not match any solver vertex DOF coord.  This indicates "
                f"the extracted wp_mesh vertices and the parent vol_mesh "
                f"BND={args.wp_label!r} vertices are not coincident.")
        wp_phi_re_at_v = phi_at_wp.real
        wp_phi_im_at_v = phi_at_wp.imag
        phi_re_int = None  # not used in this branch
        phi_im_int = None
    else:
        gf_phi_re = _GF(bem.fes)
        gf_phi_re.vec.FV().NumPy()[:] = res_bem["phi_vec"].real
        gf_phi_im = _GF(bem.fes)
        gf_phi_im.vec.FV().NumPy()[:] = res_bem["phi_vec"].imag
        phi_re_int = _Int(gf_phi_re, wp_mesh, VOL_or_BND=_BND,
                           element_wise=True)
        phi_im_int = _Int(gf_phi_im, wp_mesh, VOL_or_BND=_BND,
                           element_wise=True)
        wp_phi_re_at_v = wp_phi_im_at_v = None

    wp_centroids_l, wp_areas_l, wp_norms_l, wp_phi_avg_l = [], [], [], []
    for el in wp_mesh.Elements(_BND):
        a = abs(elem_A_int[el.nr])
        if a < 1e-30:
            continue
        pts = [wp_mesh.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(p[0], p[1], p[2]) for p in pts], axis=0)
        n_at_el = np.array([n_int[i][el.nr] / a for i in range(3)])
        if basis_order >= 2:
            corner_phi_re = np.array(
                [wp_phi_re_at_v[v.nr] for v in el.vertices])
            corner_phi_im = np.array(
                [wp_phi_im_at_v[v.nr] for v in el.vertices])
            phi_avg = corner_phi_re.mean() + 1j * corner_phi_im.mean()
        else:
            phi_avg = (phi_re_int[el.nr] + 1j * phi_im_int[el.nr]) / a
        wp_centroids_l.append(c); wp_areas_l.append(a)
        wp_norms_l.append(n_at_el); wp_phi_avg_l.append(phi_avg)
    wp_centroids_arr = np.asarray(wp_centroids_l, dtype=float)
    wp_areas_arr = np.asarray(wp_areas_l, dtype=float)
    wp_norms_arr = np.asarray(wp_norms_l, dtype=float)
    wp_phi_avg_arr = np.asarray(wp_phi_avg_l, dtype=complex)

    # 8. Telegen ΔL — branch again on coil source type
    if coil_data["source_type"] == "filament":
        from radia.workpiece_surface import delta_L_telegen_phiB
        delta_L_complex = delta_L_telegen_phiB(
            coil_data["paths"], coil_data["I_fil"],
            wp_centroids_arr, wp_areas_arr, wp_norms_arr, wp_phi_avg_arr,
            I_port=args.current)
    else:  # surface
        coil_J_complex = (complex(args.current)
                           * coil_data["coil_J_per_tri"]).astype(complex)
        delta_L_complex = _delta_L_telegen_phiB_from_surface_J(
            coil_data["coil_centroids"], coil_data["coil_areas"],
            coil_J_complex,
            wp_centroids_arr, wp_areas_arr, wp_norms_arr, wp_phi_avg_arr,
            I_port=args.current)
    delta_L_nH = float(delta_L_complex.real) * 1e9
    delta_R_mOhm = -float(delta_L_complex.imag) * omega * 1e3
    progress("BEM",
        f"P_wp={P_wp:.4e} W ΔL={delta_L_nH:+.3f} nH ΔR={delta_R_mOhm:+.4f} mΩ")

    # 9. Optional .msh output
    msh_file = ""
    if args.msh_output:
        try:
            from radia.gmsh_post_export import GmshPostExport
            post = GmshPostExport(wp_mesh, boundary=True)
            post.write(args.msh_output)
            msh_file = args.msh_output
            progress("BEM", f"wp .msh written: {args.msh_output}")
        except Exception as exc:
            progress("BEM", f"msh export failed: {exc}")

    # Z_s_wp may be scalar (default) or ndarray (--esim-per-panel).
    # Report mean Re/Im for the array case to keep the JSON schema
    # backward-compatible; full array goes under esim_per_panel_Zs.
    if isinstance(Z_s_wp, np.ndarray):
        Z_s_mean = complex(np.mean(Z_s_wp))
        Z_s_real_out = float(Z_s_mean.real)
        Z_s_imag_out = float(Z_s_mean.imag)
        per_panel_block = {
            "esim_per_panel": True,
            "esim_per_panel_Z_s_real": Z_s_wp.real.tolist(),
            "esim_per_panel_Z_s_imag": Z_s_wp.imag.tolist(),
        }
        # Save final per-DOF |H_t| at convergence (only when per-panel
        # mode was active and the final re-solve completed).
        try:
            if H_t_per_final is not None:
                per_panel_block["esim_per_panel_H_t"] = (
                    H_t_per_final.tolist())
        except NameError:
            pass
    else:
        Z_s_real_out = float(Z_s_wp.real)
        Z_s_imag_out = float(Z_s_wp.imag)
        per_panel_block = {"esim_per_panel": False}

    # ---- Spatial q_surf .sol output (consumed by calc_heat.py) -------
    # Save q_surf as an H1 GridFunction ON THE PARENT vol_mesh (not on
    # the extracted wp surface mesh).  Why this matters for the
    # downstream cross-mesh transfer:
    #
    #   * wp_mesh is a 2D SURFACE mesh.  If we naively save the GF on
    #     bem.fes (= H1(wp_mesh, order=1)) and ship the wp_mesh as
    #     em_vol, calc_heat's ``em_mesh(xw, yw, zb)`` runs against a
    #     2D mesh embedded in 3D.  NGSolve's mesh-point search on a
    #     thin surface mesh is not robust: a wp_thermal_mesh vertex
    #     that lies on the COIL-NEAR side of the workpiece can match
    #     a surface element on the COIL-FAR side (closest-distance
    #     ambiguity for a wrap-around surface), so the q_surf
    #     distribution applied to the thermal solve is mirror-flipped
    #     vs. the actual EM solution.  Symptom: temperature high on
    #     the coil-far face, low under the coil (2026-05-21 kubota
    #     bug report).
    #
    #   * Saving on H1(vol_mesh, order=1) -- the FULL 3D parent mesh
    #     used by the EM solve -- matches calc_fem_kelvin.py /
    #     calc_fem_coilmesh.py.  The em_vol the thermal panel
    #     consumes is ``args.vol`` (the same input .vol), and
    #     calc_heat's 3D mesh-point search is robust because every
    #     wp_thermal vertex maps to a real volumetric MIP in the EM
    #     mesh.
    #
    # Mapping bem.fes DOFs (wp_mesh vertex indices) to vol_mesh
    # vertex indices uses ``_wp_new_to_old`` returned from
    # ``_extract_surface_mesh_filtered`` for the material-extracted
    # case.  For ``_extract_bnd_only`` (boundary-sideset case) the
    # function does not return a vertex map, so we build one inline
    # via coordinate matching against the parent vol_mesh.
    qsurf_sol_path = ""
    qsurf_vol_path = ""
    if gf_q_wp is not None:
        try:
            from ngsolve import H1 as _H1
            # Build full-3D H1 GridFunction for the qsurf .sol.
            fes_q_vol = _H1(vol_mesh, order=1)
            gf_q_vol = _GF(fes_q_vol)
            gf_q_vol.vec[:] = 0
            vol_vec_np = gf_q_vol.vec.FV().NumPy()

            # Resolve wp_mesh vertex -> vol_mesh vertex mapping.
            if _wp_new_to_old is not None:
                # Material-extracted path: map provided.
                for wp_v_idx in range(wp_mesh.nv):
                    vol_v_idx = _wp_new_to_old.get(wp_v_idx)
                    if vol_v_idx is not None:
                        vol_vec_np[vol_v_idx] = float(q_per_dof[wp_v_idx])
                n_matched = len(_wp_new_to_old)
            else:
                # Boundary-extracted path: build mapping by coord
                # matching against the parent vol_mesh.  Round to
                # 12 sig figs (NGSolve's Save/Load tolerance for
                # vertex coordinates).
                vol_coord_to_idx = {}
                for vol_v_idx in range(vol_mesh.nv):
                    p = vol_mesh.vertices[vol_v_idx].point
                    key = (round(p[0], 12), round(p[1], 12),
                           round(p[2], 12))
                    vol_coord_to_idx[key] = vol_v_idx
                n_matched = 0
                for wp_v_idx in range(wp_mesh.nv):
                    p = wp_mesh.vertices[wp_v_idx].point
                    key = (round(p[0], 12), round(p[1], 12),
                           round(p[2], 12))
                    vol_v_idx = vol_coord_to_idx.get(key)
                    if vol_v_idx is not None:
                        vol_vec_np[vol_v_idx] = float(q_per_dof[wp_v_idx])
                        n_matched += 1

            # Output paths: prefer args.msh_output's basename when
            # the user requested an .msh, else args.vol's basename.
            if args.msh_output:
                base_dir = os.path.dirname(
                    os.path.abspath(args.msh_output)) or "."
                stem = os.path.splitext(
                    os.path.basename(args.msh_output))[0]
            else:
                base_dir = os.path.dirname(
                    os.path.abspath(args.vol)) or "."
                stem = os.path.splitext(
                    os.path.basename(args.vol))[0]
            os.makedirs(base_dir, exist_ok=True)
            sol_Q = os.path.join(base_dir,
                                  f"{stem}_qsurf.sol").replace("\\", "/")
            gf_q_vol.Save(sol_Q)
            qsurf_sol_path = sol_Q
            # em_vol pair is the input .vol itself -- the same mesh
            # the EM solve was on.  No separate _bem.vol needed.
            qsurf_vol_path = os.path.abspath(args.vol).replace("\\", "/")
            progress("BEM",
                f"wrote qsurf.sol on parent vol_mesh: "
                f"{os.path.basename(sol_Q)} "
                f"({n_matched}/{wp_mesh.nv} wp vertices mapped, "
                f"em_vol={os.path.basename(qsurf_vol_path)})")
        except Exception as e:
            import traceback
            progress("BEM",
                f"qsurf.sol save failed: {type(e).__name__}: {e}")
            for line in traceback.format_exc().splitlines()[-4:]:
                progress("BEM", f"  {line}")

    return {
        "P_wp": P_wp, "H_t_rms": H_t_rms, "A_wp": A_wp,
        "delta_L_nH": delta_L_nH, "delta_R_mOhm": delta_R_mOhm,
        "wp_dissipation_R_mOhm": wp_dissipation_R_mOhm,
        "wp_mesh_nv": int(wp_mesh.nv),
        "wp_euler_chi": int(wp_chi),
        "wp_genus": int(wp_genus),
        "wp_loop_dof_mode": wp_loop_req,
        "wp_loop_dof_skip_reason": wp_loop_skip,
        **(loop_meta or {}),
        "wp_mesh_n_tris": int(wp_mesh.GetNE(BND)),
        "wp_ndof": int(bem.ndof),
        "wp_basis_order": int(basis_order),
        "wp_vol_curve_order": int(vol_curve_order),
        "wp_gmres_info": int(info),
        "Z_s_wp_real": Z_s_real_out,
        "Z_s_wp_imag": Z_s_imag_out,
        "skin_depth_wp_mm": float(delta_wp * 1e3),
        "impedance_model": args.impedance_model,
        "esim_iterations": int(n_iter_done),
        "esim_converged": bool(esim_converged),
        "esim_anderson_m": int(getattr(args, "esim_anderson_m", 0)),
        "esim_anderson_restarts": int(anderson.n_restarts) if esim_solver is not None else 0,
        "esim_anderson_clips": int(anderson.n_clips) if esim_solver is not None else 0,
        "esim_history": esim_history,
        **per_panel_block,
        "msh_file": msh_file,
        "qsurf_sol": qsurf_sol_path,
        "qsurf_em_vol": qsurf_vol_path,
        "wp_phi_inc": phi_inc_mode,
        "wp_phi_inc_grad_residual": (
            float(phi_inc_grad_residual)
            if phi_inc_grad_residual is not None else None),
        "t_wp_mesh_s": t_wp_mesh,
        "t_bem_assembly_s": t_asm,
        "t_phi_inc_s": t_phi,
        "t_bem_solve_s": t_bie,
    }


def _extract_bnd_only_inline(vol_mesh, bnd_label):
    """Hole-extractor for wp-as-hole geometry (closed-torus convention).

    Walks the .vol mesh, collects BND elements whose ``el.index``
    matches a FaceDescriptor with name ``bnd_label``, and emits a
    closed surface mesh with outward-oriented triangles (centroid-
    aligned outward).  Multiple FDs may share the same name after
    webcut (e.g. air_top + air_bot halves) and the per-FD winding can
    flip relative to its own domin -- a blanket flip would be wrong.

    Originally lived in ``calc_peec_bem.py``; inlined here when the
    three legacy CLIs were unified into ``calc_inductance.py``.
    """
    from ngsolve import Mesh, BND
    import netgen.meshing as ngm

    bnd_labels = list(vol_mesh.GetBoundaries())
    target_indices = {i for i, n in enumerate(bnd_labels) if n == bnd_label}
    if not target_indices:
        raise ValueError(
            f"Boundary label {bnd_label!r} not found in .vol.  "
            f"Available: {sorted(set(bnd_labels))}")

    ngmesh_new = ngm.Mesh(dim=3)
    ngmesh_new.Add(ngm.FaceDescriptor(bc=1, domin=1))
    ngmesh_new.SetBCName(0, bnd_label)

    used_vtx = set()
    for el in vol_mesh.Elements(BND):
        if el.index in target_indices:
            for v in el.vertices:
                used_vtx.add(v.nr)
    if not used_vtx:
        raise ValueError(
            f"No BND elements found with label {bnd_label!r}")

    old_to_new = {}
    for old in sorted(used_vtx):
        p = vol_mesh.vertices[old].point
        old_to_new[old] = ngmesh_new.Add(
            ngm.MeshPoint(ngm.Pnt(p[0], p[1], p[2])))

    # Collect triangles, then make the winding GLOBALLY consistent +
    # outward via face-BFS + per-component signed volume.  The old
    # per-triangle centroid test breaks genus-1 workpieces: the
    # bore-wall outward normal points TOWARD the centroid, so the whole
    # inner wall is flipped and the double-layer operator is corrupted.
    from surface_mesh_extract import orient_surface_triangles
    order_vtx = sorted(used_vtx)
    compact = {v: i for i, v in enumerate(order_vtx)}
    tri_rows = []
    for el in vol_mesh.Elements(BND):
        if el.index not in target_indices:
            continue
        tri_rows.append([compact[v.nr] for v in el.vertices])
    pts_used = np.array([vol_mesh.vertices[v].point for v in order_vtx])
    tri_oriented, stats = orient_surface_triangles(pts_used, tri_rows)
    for t in tri_oriented:
        ngmesh_new.Add(ngm.Element2D(
            1, [old_to_new[order_vtx[int(c)]] for c in t]))
    progress("BEM",
        f"winding: {stats['n_flips']} flips, "
        f"{stats['components_flipped']}/{stats['n_components']} components "
        f"reversed, conflicts {stats['conflicts_before']} -> "
        f"{stats['conflicts_after']}")
    return Mesh(ngmesh_new)


# ======================================================================
# Output assembly
# ======================================================================
def _assemble_vacuum_output(args, coil_data):
    """L_coil + R_coil only (no workpiece)."""
    out = {
        "status": "ok",
        "method": f"{args.coil_solver}-vacuum",
        "coil_solver_used": args.coil_solver,
        "frequency_hz": float(args.frequency),
        "current_A": float(args.current),
        "L_coil_nH": float(coil_data["L_coil"] * 1e9),
        "R_coil_mOhm": float(coil_data["R_coil"] * 1e3),
        "t_coil_topology_s": float(coil_data["t_coil_topology_s"]),
        "t_coil_solve_s": float(coil_data["t_coil_solve_s"]),
    }
    if args.coil_solver == "peec":
        out["n_filaments"] = int(coil_data["n_filaments"])
    else:
        out["coil_mesh_nv"] = int(coil_data["coil_mesh_nv"])
        out["coil_mesh_n_tris"] = int(coil_data["coil_mesh_n_tris"])
        out["coil_mesh_n_src_tris"] = int(coil_data["coil_mesh_n_src_tris"])
        out["coil_mesh_n_snk_tris"] = int(coil_data["coil_mesh_n_snk_tris"])
        # BEM-A saddle DoF breakdown (matched to the "BEMA coil DoF" log
        # line so users grepping the JSON have the same numbers).
        out["coil_n_J"] = int(coil_data.get("coil_n_J", 0))
        out["coil_n_f"] = int(coil_data.get("coil_n_f", 0))
        out["coil_saddle_ndof"] = int(coil_data.get("coil_saddle_ndof", 0))
        out["bem_a_residual"] = float(coil_data["bem_a_residual"])
    return out


def _assemble_full_output(args, coil_data, wp_data):
    """L_coil + R_coil + workpiece weak-coupled outputs."""
    out = _assemble_vacuum_output(args, coil_data)
    out["method"] = f"{args.coil_solver}-bem-weak"
    out["coupling_mode"] = "weak"
    out["telegen_form"] = "phi-B"
    out["L_total_nH"] = (float(coil_data["L_coil"] * 1e9)
                          + wp_data["delta_L_nH"])
    out["R_total_mOhm"] = (float(coil_data["R_coil"] * 1e3)
                            + wp_data["delta_R_mOhm"])
    out["delta_L_nH"] = wp_data["delta_L_nH"]
    out["delta_R_mOhm"] = wp_data["delta_R_mOhm"]
    out["P_wp_W"] = wp_data["P_wp"]
    out["H_t_rms_A_per_m"] = wp_data["H_t_rms"]
    out["wp_area_m2"] = wp_data["A_wp"]
    # Spatial q_surf .sol pair for downstream thermal analysis.  Empty
    # strings when basis_order>=2 (qsurf save skipped: phi_vec on parent
    # vol_mesh P2 needs explicit coord matching; out of scope for now).
    out["qsurf_sol"] = wp_data.get("qsurf_sol", "")
    out["qsurf_em_vol"] = wp_data.get("qsurf_em_vol", "")
    out["wp_dissipation_R_mOhm"] = wp_data["wp_dissipation_R_mOhm"]
    out["wp_mesh_nv"] = wp_data["wp_mesh_nv"]
    out["wp_mesh_n_tris"] = wp_data["wp_mesh_n_tris"]
    # Topology context + genus caveat (see _wp_genus_check): on a genus>=1
    # workpiece the scalar BIE drops the shorted-turn eddy current, so the
    # absolute P_wp / H_t carry a known over-estimate -- UNLESS the loop
    # DOF is active (--wp-loop-dof), which solves that current explicitly.
    out["wp_euler_chi"] = int(wp_data.get("wp_euler_chi", 2))
    out["wp_genus"] = int(wp_data.get("wp_genus", 0))
    out["wp_loop_dof_mode"] = wp_data.get("wp_loop_dof_mode", "auto")
    if wp_data.get("wp_loop_dof_skip_reason"):
        out["wp_loop_dof_skip_reason"] = wp_data["wp_loop_dof_skip_reason"]
    if wp_data.get("wp_loop_dof"):
        out["wp_loop_dof"] = True
        out["wp_loop_alpha_A"] = float(wp_data["wp_loop_alpha_A"])
        out["wp_loop_alpha_deg"] = float(wp_data["wp_loop_alpha_deg"])
        out["wp_loop_theta_jump"] = float(wp_data["wp_loop_theta_jump"])
        out["wp_loop_cut_n_vertices"] = int(
            wp_data["wp_loop_cut_n_vertices"])
        out["t_loop_dof_s"] = float(wp_data["t_loop_dof_s"])
        out["P_wp_note"] = (
            "genus-1 loop DOF active: the net shorted-turn eddy current "
            "(wp_loop_alpha_A) is solved and its Lenz screening is included "
            "in P_wp / H_t.  delta_L keeps the plain-solve phi convention; "
            "the loop extension is locked by the analytic shorted-ring "
            "golden.")
    elif out["wp_genus"] != 0:
        out["P_wp_caveat"] = _GENUS_P_WP_CAVEAT
    # workpiece BIE DoF context -- mirror of the "BEM ndof=..." log line
    # so the JSON retains the same numbers users see in the console log.
    out["wp_ndof"] = int(wp_data.get("wp_ndof", 0))
    out["wp_basis_order"] = int(wp_data.get("wp_basis_order", 1))
    out["wp_vol_curve_order"] = int(wp_data.get("wp_vol_curve_order", 1))
    out["wp_gmres_info"] = wp_data["wp_gmres_info"]
    out["Z_s_wp_real"] = wp_data["Z_s_wp_real"]
    out["Z_s_wp_imag"] = wp_data["Z_s_wp_imag"]
    out["skin_depth_wp_mm"] = wp_data["skin_depth_wp_mm"]
    out["impedance_model"] = wp_data["impedance_model"]
    out["esim_iterations"] = wp_data["esim_iterations"]
    out["esim_converged"] = wp_data["esim_converged"]
    out["esim_anderson_m"] = wp_data.get("esim_anderson_m", 0)
    out["esim_anderson_restarts"] = wp_data.get("esim_anderson_restarts", 0)
    out["esim_anderson_clips"] = wp_data.get("esim_anderson_clips", 0)
    out["esim_history"] = wp_data["esim_history"]
    out["msh_file"] = wp_data["msh_file"]
    out["t_wp_mesh_s"] = wp_data["t_wp_mesh_s"]
    out["t_bem_assembly_s"] = wp_data["t_bem_assembly_s"]
    out["t_phi_inc_s"] = wp_data["t_phi_inc_s"]
    out["wp_phi_inc"] = wp_data.get("wp_phi_inc", "path")
    if wp_data.get("wp_phi_inc_grad_residual") is not None:
        out["wp_phi_inc_grad_residual"] = float(
            wp_data["wp_phi_inc_grad_residual"])
    out["t_bem_solve_s"] = wp_data["t_bem_solve_s"]
    # Propagate per-panel block (esim_per_panel, esim_per_panel_Z_s_real/imag,
    # esim_per_panel_H_t) — emitted by _solve_workpiece_weak_coupled when
    # --esim-per-panel is set.
    for key in ("esim_per_panel",
                "esim_per_panel_Z_s_real",
                "esim_per_panel_Z_s_imag",
                "esim_per_panel_H_t"):
        if key in wp_data:
            out[key] = wp_data[key]
    return out


# ======================================================================
# Strong-coupled workpiece block (iterative CoupledBEMSolver)
#
# Re-exposes the validated per-DOF back-reaction solver
# (radia.bem_coupled_solver.CoupledBEMSolver) as the --coupling-mode
# strong path.  Unlike the weak Telegen path, the coil surface current is
# recomputed each Picard iteration in response to the workpiece reaction
# field (A_wp projected per-DOF onto the coil), so ΔL includes the
# workpiece magnetic-energy term (the term the weak Telegen form drops)
# and P_wp is self-consistent.  Coil L/R still come from the BEM-A
# impedance-EFIE coil solve (coil_data); the coupled solver contributes
# ΔL, ΔR (from P_wp) and H_t.  Cross-checked vs FEM-Kelvin SIBC in
# validation_test/induction_heating (copper +0.3%, steel μr=100 +1.7% on L).
# ======================================================================
def _solve_workpiece_strong_coupled(args):
    """Iterative coil<->workpiece BEM coupling via CoupledBEMSolver.

    Requires ``--coil-solver bem-a`` (the coupled solver needs a coil
    SURFACE mesh with source/sink labels to drive the coil EFIE; PEEC
    filaments cannot).  Returns a dict with the coupled quantities
    (L_air, L_total, Delta_L, P_total, H_t_rms, iterations, ...) plus mesh
    and timing context, ready for ``_assemble_strong_output``.
    """
    from ngsolve import Mesh, BND
    from surface_mesh_extract import _extract_surface_mesh_filtered
    from radia.bem_coupled_solver import CoupledBEMSolver
    from em_material import EMMaterial

    omega = 2.0 * math.pi * args.frequency

    # --- coil surface mesh (CoupledBEMSolver re-assembles the coil EFIE
    #     saddle internally; we only hand it the mesh + source/sink names) ---
    t0 = time.perf_counter()
    progress("COUPLED", f"load coil surface from {os.path.basename(args.coil_vol)}")
    # Use the SAME loader as the vacuum BEM-A path: it validates the
    # source/sink labels AND converts a volume .vol to the pure surface
    # mesh HDivSurface needs.  Loading Mesh(args.coil_vol) directly here
    # reproduces a singular coupled saddle on volume-tet inputs because
    # HDivSurface also sees the volume modes.
    coil_mesh = _build_bema_coil_mesh(args)
    # Same fail-loud source/sink check as the vacuum BEM-A path: the coupled
    # solver drives the coil EFIE from these labels, so a shortcut placement
    # would silently poison L_air / Delta_L / P_wp too.
    from radia.bem.coil_inductance_ngsolve import check_source_sink_current_path
    _cv, _ct, _sm, _km = _arrays_from_bema_coil_mesh(
        coil_mesh, source_name=args.coil_source_name,
        sink_name=args.coil_sink_name)
    check_source_sink_current_path(
        _cv, _ct, _sm, _km, source_label=args.coil_source_name,
        sink_label=args.coil_sink_name)

    # --- workpiece surface mesh (same extraction convention as the weak
    #     path: --wp-label may name a material or a boundary sideset) ---
    progress("COUPLED", f"load wp surface from {os.path.basename(args.vol)}")
    vol_mesh = Mesh(args.vol)
    mats = set(vol_mesh.GetMaterials())
    bnds = set(vol_mesh.GetBoundaries())
    if args.wp_label in mats:
        wp_mesh, _ = _extract_surface_mesh_filtered(
            vol_mesh, keep_label=args.wp_label, return_vertex_map=True)
    elif args.wp_label in bnds:
        wp_mesh = _extract_bnd_only_inline(vol_mesh, args.wp_label)
    else:
        raise ValueError(
            f"wp_label {args.wp_label!r} not found in materials "
            f"({sorted(mats)}) or boundaries ({sorted(bnds)})")
    t_wp_mesh = time.perf_counter() - t0
    wp_chi, wp_genus = _wp_genus_check(wp_mesh, tag="COUPLED")

    # --- global Leontovich Z_s (delta includes mu_r), same convention as
    #     the linear weak path: Z_s = (1+j) * rho / delta ---
    mat_wp = EMMaterial(name="wp", sigma=args.sigma, mu_r=args.mu_r)
    delta_wp = mat_wp.skin_depth(args.frequency)
    Z_s = (1.0 + 1j) * (1.0 / args.sigma) / delta_wp

    # Workpiece BIE backend: HACApK (O(N log N), scales past the ~12k-tri
    # dense wall) by default; --wp-bem-backend intree-dense forces the dense
    # path (small wp / debug).  Uses the existing weak-path flag.
    wp_hacapk = (args.wp_bem_backend == "hacapk")
    # Coil EFIE backend: dense-LU by default (fastest below ~12k n_J);
    # --coil-saddle-solver hacapk_cocr compresses the coil SL to an O(N log N)
    # H-matrix + loop-COCR saddle solve (O(N r) storage, for large coils).
    coil_hacapk = (args.coil_saddle_solver == "hacapk_cocr")
    progress("COUPLED",
        f"strong coupling: coil {coil_mesh.nv}v / wp {wp_mesh.nv}v, "
        f"f={args.frequency:g} Hz, Z_s={Z_s:.3e}, "
        f"coil_backend={'hacapk' if coil_hacapk else 'dense-lu'}, "
        f"wp_backend={'hacapk' if wp_hacapk else 'dense'}, "
        f"max_iter={args.coupling_max_iter} tol={args.coupling_tol:g} "
        f"relax={args.coupling_relax:g}")
    t0 = time.perf_counter()
    solver = CoupledBEMSolver(
        coil_mesh, wp_mesh,
        source_label=args.coil_source_name,
        sink_label=args.coil_sink_name,
        fes_order=0,
        wp_hacapk=wp_hacapk, wp_aca_eps=args.wp_aca_eps,
        wp_gmres_tol=args.wp_gmres_tol,
        coil_hacapk=coil_hacapk, coil_aca_eps=args.coil_aca_eps)
    sol = solver.solve(
        Z_s=Z_s, omega=omega,
        max_iter=int(args.coupling_max_iter),
        tol=float(args.coupling_tol),
        relax=float(args.coupling_relax))
    t_solve = time.perf_counter() - t0

    # CoupledBEMSolver drives the coil EFIE at UNIT terminal current (its
    # g_red constraint is normalised; solve() takes no current).  The
    # problem is linear, so scale the current-dependent outputs to
    # --current here: fields ~ I, dissipation ~ I^2.  Inductances
    # (L_air / L_total / Delta_L) are per-unit-current energies and stay.
    # Missing this scale was the Kubota 2026-07-16 incident: at 6700 A the
    # panel reported P_wp = 8.3e-4 W / H_t = 9.84 A/m -- 1 A-drive values,
    # I^2 = 4.5e7x low -- while Delta_L came out sane, masking the bug.
    # (The PEEC strong driver passes I_port to its solver, so both drivers
    # now return REAL-current-scaled P_total / H_t_rms / J arrays to the
    # shared _assemble_strong_output.)
    I_term = float(args.current)
    sol["P_total"] = float(sol["P_total"]) * I_term * I_term
    sol["H_t_rms"] = float(sol["H_t_rms"]) * abs(I_term)
    for key in ("wp_J_re", "wp_J_im", "J_coil_re", "J_coil_im"):
        if sol.get(key) is not None:
            sol[key] = sol[key] * I_term

    progress("COUPLED",
        f"converged in {int(sol['iterations'])} iter: "
        f"L_total={sol['L_total'] * 1e9:.2f} nH "
        f"dL={sol['Delta_L'] * 1e9:+.2f} nH "
        f"P_wp={sol['P_total']:.4e} W H_t={sol['H_t_rms']:.2f} A/m "
        f"@ I={I_term:g} A ({t_solve:.1f}s)")

    sol["wp_mesh_nv"] = int(wp_mesh.nv)
    sol["wp_mesh_n_tris"] = int(wp_mesh.GetNE(BND))
    sol["wp_euler_chi"] = int(wp_chi)
    sol["wp_genus"] = int(wp_genus)
    sol["delta_wp_mm"] = float(delta_wp * 1e3)
    sol["t_wp_mesh_s"] = float(t_wp_mesh)
    sol["t_coupled_solve_s"] = float(t_solve)
    sol["wp_hacapk"] = bool(wp_hacapk)
    sol["coil_hacapk"] = bool(coil_hacapk)
    return sol


def _solve_workpiece_strong_coupled_peec(args, coil_data):
    """Iterative PEEC filament coil <-> workpiece BEM coupling.

    The PEEC analogue of ``_solve_workpiece_strong_coupled``: rebuilds the
    filament loop-bundle from ``--coil-step`` and drives
    ``CoupledPEECBEMSolver`` (back-EMF from the workpiece reaction field
    redistributes the filament currents self-consistently).  Returns the
    same key set as the BEM-A strong path so ``_assemble_strong_output``
    renders both uniformly; R is taken from ``2 P_wp / I^2`` (energy),
    ``Delta_L`` from the coupled solve.

    EXPERIMENTAL: the strong-loading response does not yet have a durable
    independent reference case (the committed demo is weakly coupled).  See
    ``radia.peec_coupled_bem_solver`` module docstring.
    """
    from ngsolve import Mesh, BND
    from surface_mesh_extract import _extract_surface_mesh_filtered
    from radia.coil_from_cad import filaments_from_step
    from radia.peec_bundle import build_loop_bundle_impedance
    from radia.peec_coupled_bem_solver import CoupledPEECBEMSolver
    from em_material import EMMaterial

    omega = 2.0 * math.pi * args.frequency

    # --- PEEC coil loop-bundle (same construction as _solve_coil_peec, but
    #     proximity-iterative is OFF for the coupled solve: the workpiece
    #     back-reaction is the coupling of interest here) ---
    progress("COUPLED", f"PEEC bundle from {os.path.basename(args.coil_step)} "
                        f"(n_peri={args.peec_n_peri})")
    t0 = time.perf_counter()
    topo = filaments_from_step(args.coil_step, sigma=args.coil_sigma,
                               n_peri=args.peec_n_peri, use_coil_builder=True)
    paths = topo["filament_paths"]
    R_f, L_f = build_loop_bundle_impedance(topo["solver"],
                                           topo["seg_of_filament"])
    a_eq = _equivalent_wire_radius_m(topo, args.peec_n_peri)
    Zs_fil = _peec_skin_impedance_per_filament(
        paths, a_eq, args.coil_sigma, omega, args.peec_n_peri)

    # --- workpiece surface mesh (same convention as the BEM-A strong path) ---
    vol_mesh = Mesh(args.vol)
    mats = set(vol_mesh.GetMaterials())
    bnds = set(vol_mesh.GetBoundaries())
    if args.wp_label in mats:
        wp_mesh, _ = _extract_surface_mesh_filtered(
            vol_mesh, keep_label=args.wp_label, return_vertex_map=True)
    elif args.wp_label in bnds:
        wp_mesh = _extract_bnd_only_inline(vol_mesh, args.wp_label)
    else:
        raise ValueError(
            f"wp_label {args.wp_label!r} not found in materials "
            f"({sorted(mats)}) or boundaries ({sorted(bnds)})")
    t_wp_mesh = time.perf_counter() - t0
    wp_chi, wp_genus = _wp_genus_check(wp_mesh, tag="COUPLED")

    mat_wp = EMMaterial(name="wp", sigma=args.sigma, mu_r=args.mu_r)
    delta_wp = mat_wp.skin_depth(args.frequency)
    Z_s = (1.0 + 1j) * (1.0 / args.sigma) / delta_wp
    wp_hacapk = (args.wp_bem_backend == "hacapk")

    progress("COUPLED",
        f"strong PEEC coupling [EXPERIMENTAL]: {len(paths)} filaments / "
        f"wp {wp_mesh.nv}v, f={args.frequency:g} Hz, Z_s={Z_s:.3e}, "
        f"wp_backend={'hacapk' if wp_hacapk else 'dense'}, "
        f"max_iter={args.coupling_max_iter} tol={args.coupling_tol:g} "
        f"relax={args.coupling_relax:g}")
    progress("COUPLED",
        "WARNING: strong PEEC loading response has no durable independent "
        "reference case (demo is weakly coupled) -- treat absolute P_wp/dL as "
        "unverified.")
    t0 = time.perf_counter()
    solver = CoupledPEECBEMSolver(
        paths, R_f, L_f, wp_mesh, Zs_fil=Zs_fil, wp_order=1,
        wp_hacapk=wp_hacapk, wp_aca_eps=args.wp_aca_eps,
        wp_gmres_tol=args.wp_gmres_tol)
    sol = solver.solve(
        Z_s=Z_s, omega=omega, I_port=float(args.current),
        max_iter=int(args.coupling_max_iter),
        tol=float(args.coupling_tol),
        relax=float(args.coupling_relax))
    t_solve = time.perf_counter() - t0
    progress("COUPLED",
        f"converged in {int(sol['iterations'])} iter: "
        f"L_total={sol['L_total'] * 1e9:.2f} nH "
        f"dL={sol['Delta_L'] * 1e9:+.2f} nH "
        f"P_wp={sol['P_total']:.4e} W H_t={sol['H_t_rms']:.2f} A/m "
        f"({t_solve:.1f}s)")

    # Context keys expected by _assemble_strong_output (mirror the BEM-A
    # driver).  n_J_coil -> filament count; the coil backend is the PEEC
    # loop-bundle (no H-matrix coil), so coil_hacapk is always False.
    sol["n_J_coil"] = int(sol["n_filaments"])
    sol["wp_mesh_nv"] = int(wp_mesh.nv)
    sol["wp_mesh_n_tris"] = int(wp_mesh.GetNE(BND))
    sol["wp_euler_chi"] = int(wp_chi)
    sol["wp_genus"] = int(wp_genus)
    sol["delta_wp_mm"] = float(delta_wp * 1e3)
    sol["t_wp_mesh_s"] = float(t_wp_mesh)
    sol["t_coupled_solve_s"] = float(t_solve)
    sol["wp_hacapk"] = bool(wp_hacapk)
    sol["coil_hacapk"] = False
    return sol


def _assemble_strong_output(args, coil_data, strong):
    """L_coil + R_coil (BEM-A) + strong-coupled ΔL / ΔR / P_wp / H_t.

    Mirrors the key names of ``_assemble_full_output`` (L_total_nH,
    R_total_mOhm, delta_L_nH, delta_R_mOhm, P_wp_W, H_t_rms_A_per_m) so the
    panel summary renders strong-coupling runs the same way as weak ones,
    with ``coupling_mode = "strong"``.
    """
    out = _assemble_vacuum_output(args, coil_data)
    out["method"] = f"{args.coil_solver}-bem-strong"
    out["coupling_mode"] = "strong"
    I = float(args.current)
    P_wp = float(strong["P_total"])
    dL_nH = float(strong["Delta_L"] * 1e9)
    # ΔR from the self-consistent workpiece dissipation: P_wp = 1/2 I^2 ΔR.
    dR_mOhm = (2.0 * P_wp / (I * I) * 1e3) if I else 0.0
    out["L_total_nH"] = float(coil_data["L_coil"] * 1e9) + dL_nH
    out["R_total_mOhm"] = float(coil_data["R_coil"] * 1e3) + dR_mOhm
    out["delta_L_nH"] = dL_nH
    out["delta_R_mOhm"] = dR_mOhm
    out["P_wp_W"] = P_wp
    out["H_t_rms_A_per_m"] = float(strong["H_t_rms"])
    # Coupled-solver-native quantities (PEC-coil vacuum L + self-consistent
    # total L) for cross-checking against the impedance-EFIE coil L_coil.
    out["coupled_L_air_nH"] = float(strong["L_air"] * 1e9)
    out["coupled_L_total_nH"] = float(strong["L_total"] * 1e9)
    out["coupling_iterations"] = int(strong["iterations"])
    out["wp_ndof"] = int(strong["n_phi_wp"])
    out["coupled_n_J_coil"] = int(strong["n_J_coil"])
    out["wp_mesh_nv"] = int(strong["wp_mesh_nv"])
    out["wp_mesh_n_tris"] = int(strong["wp_mesh_n_tris"])
    out["Z_s_wp_real"] = float(strong["Z_s"].real)
    out["Z_s_wp_imag"] = float(strong["Z_s"].imag)
    out["skin_depth_wp_mm"] = float(strong["delta_wp_mm"])
    out["impedance_model"] = "sibc"
    out["wp_bem_backend"] = "hacapk" if strong.get("wp_hacapk") else "intree-dense"
    out["coil_bem_backend"] = "hacapk_cocr" if strong.get("coil_hacapk") else "dense-lu"
    out["t_wp_mesh_s"] = float(strong["t_wp_mesh_s"])
    out["t_coupled_solve_s"] = float(strong["t_coupled_solve_s"])
    # Topology context + genus-conditional P_wp caveat: the scalar BIE is
    # locked by the analytic genus-0 sphere, but on a genus>=1 workpiece it
    # cannot carry the net shorted-turn eddy current on the handle, so its
    # Lenz screening is lost and P_wp / H_t are over-estimated when the
    # coil flux links the handle.  Strong coupling does not add that
    # workpiece-side loop DOF.
    out["wp_euler_chi"] = int(strong.get("wp_euler_chi", 2))
    out["wp_genus"] = int(strong.get("wp_genus", 0))
    if out["wp_genus"] != 0:
        out["P_wp_caveat"] = _GENUS_P_WP_CAVEAT
    # PEEC filament coil: relabel the coil backend and flag the path as
    # experimental (the strong-loading response has no durable independent
    # reference case; the committed demo is weakly coupled).  Delta_R
    # from the coupled Z_port is surfaced only as a diagnostic (R itself is
    # taken from 2 P_wp / I^2 above, like the BEM-A path).
    if coil_data.get("source_type") == "filament":
        out["coil_bem_backend"] = "peec-loop-bundle"
        out["coupled_n_filaments"] = int(strong.get("n_filaments", 0))
        out["coupling_converged"] = bool(strong.get("converged", False))
        out["coupling_residual"] = float(strong.get(
            "coupling_residual", float("inf")))
        if "Delta_R" in strong:
            out["coupled_delta_R_reaction_mOhm"] = float(strong["Delta_R"] * 1e3)
        out["experimental"] = True
        out["experimental_note"] = (
            "strong PEEC coupling: the coil-current redistribution under "
            "strong loading has no durable independent reference case; the "
            "committed demo is weakly coupled. "
            "Treat absolute P_wp / delta_L as unverified.")
    return out


# ======================================================================
# Coil viz: filament .msh + B-on-bbox / B-on-air-volume .msh
#
# Restored 2026-05-08 -- v4.25.0 unification refactor (commit 16e7e6e7)
# accidentally dropped the visualization branch from the deleted
# calc_peec_inductance.py.  The panel passes --coil-msh-output but the
# unified script silently ignored it.  open_gmsh.py auto-merges any
# sibling <stem>_filaments.msh so we just need to emit the files.
# ======================================================================
def _build_coil_bbox_mesh(filament_paths, expand=2.5, n_cells=18):
    """Generate an OCC Box mesh enclosing the coil bbox * ``expand``.

    Used for the auto-field-domain visualisation when no air-volume
    .vol is available (vacuum case) OR when the user wants a clean
    Box-shaped overlay on top of the workpiece air domain.
    """
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh

    pts_list = []
    for path in filament_paths:
        arr = np.asarray(path, dtype=float)
        if arr.ndim != 3 or arr.shape[-1] != 3:
            continue
        pts_list.append(arr.reshape(-1, 3))
    if not pts_list:
        raise ValueError("_build_coil_bbox_mesh: no valid filament paths.")
    pts = np.concatenate(pts_list, axis=0)
    pmin = pts.min(axis=0)
    pmax = pts.max(axis=0)
    centre = 0.5 * (pmin + pmax)
    half = 0.5 * (pmax - pmin)
    half_exp = np.maximum(half * expand, half.max() * 0.3)
    p_lo = centre - half_exp
    p_hi = centre + half_exp
    box = Box(Pnt(float(p_lo[0]), float(p_lo[1]), float(p_lo[2])),
              Pnt(float(p_hi[0]), float(p_hi[1]), float(p_hi[2])))
    box.faces.name = "outer"
    box.solids.name = "field_domain"
    geo = OCCGeometry(box)
    maxh = float(np.max(p_hi - p_lo)) / max(n_cells, 4)
    with TaskManager():
        return Mesh(geo.GenerateMesh(maxh=maxh))


def _biot_savart_B_on_vertices(mesh, filament_paths, fil_currents):
    """Vectorised Biot-Savart B at every vertex of ``mesh``.

    Returns (n_verts, 3) numpy array of B [T] using |I| (peak amplitude
    of the complex phasor) so the magnitude pattern matches the physical
    AC peak field.
    """
    from biot_savart import h_segments_batch, MU0
    verts = np.asarray([list(v.point) for v in mesh.vertices], dtype=float)
    if verts.shape[1] != 3:
        raise ValueError("_biot_savart_B_on_vertices: 3D mesh required")
    H_total = np.zeros_like(verts)
    for path, I in zip(filament_paths, fil_currents):
        segments = np.asarray(path, dtype=float)
        if segments.ndim != 3 or segments.shape[1:] != (2, 3):
            continue
        Im = float(abs(complex(I)))
        if Im < 1e-30:
            continue
        H_total += h_segments_batch(segments, verts, current=Im)
    return MU0 * H_total


def _resolve_coil_msh_path(args):
    """Decide where to write the coil viz .msh.

    Priority:
      1. ``--coil-msh-output`` (explicit)
      2. ``--msh-output`` if vacuum / coil-only (no workpiece — the wp
         surface output is NOT competing for that path)
      3. ``<msh_output_stem>_coil.msh`` sibling (with-wp case — keeps
         the wp surface at ``args.msh_output``, panel button opens
         the coil viz via gmsh_file priority and open_gmsh.py auto-
         merges the ``_filaments.msh`` sibling)

    Returns "" if no viz path can be derived (no msh-output flags).
    """
    if args.coil_msh_output:
        return args.coil_msh_output
    if not args.msh_output:
        return ""
    if args.coil_only or not args.vol:
        return args.msh_output
    base, ext = os.path.splitext(args.msh_output)
    return base + "_coil" + (ext or ".msh")


def _export_coil_msh_viz(args, coil_data):
    """Write the PEEC/coil visualization .msh files.

    Writes (when a coil viz path can be derived from ``--coil-msh-output``
    or ``--msh-output``, see ``_resolve_coil_msh_path``):
      * ``<coil_msh_path>``         : auto-bbox volume mesh + B vector
      * ``<stem>_filaments.msh``    : 1D filament polylines + |I|
        (auto-merged by open_gmsh.py as a sibling)

    Returns ``{"coil_field_msh": ..., "filament_msh": ..., "field_n_vertices": ...}``
    or empty dict on failure.  Errors are logged via ``progress`` and do
    NOT abort the inductance solve.

    PEEC-only.  The BEM-A coil path returns surface current density on
    a triangle mesh; that has its own viz route (TODO, not in scope).
    """
    if coil_data.get("source_type") != "filament":
        return {}

    out_msh = _resolve_coil_msh_path(args)
    if not out_msh:
        return {}

    paths = coil_data["paths"]
    I_fil = coil_data["I_fil"]
    base_dir = os.path.dirname(os.path.abspath(out_msh)) or "."
    stem = os.path.splitext(os.path.basename(out_msh))[0]
    os.makedirs(base_dir, exist_ok=True)

    # ---- B field on auto-bbox mesh -----------------------------------
    coil_field_msh = ""
    field_n_vertices = 0
    try:
        from ngsolve import H1, GridFunction
        from ngsolve import TaskManager
        from gmsh_post_export import save_vol_sol_pair, vol2msh

        t0 = time.perf_counter()
        field_mesh = _build_coil_bbox_mesh(paths)
        field_n_vertices = field_mesh.nv
        B_arr = _biot_savart_B_on_vertices(field_mesh, paths, I_fil)
        fes_B = H1(field_mesh, order=1, dim=3)
        gf_B = GridFunction(fes_B)
        gf_B.vec.FV().NumPy()[:] = B_arr.reshape(-1)

        sol_B = os.path.join(base_dir, f"{stem}_B.sol").replace("\\", "/")
        vol_B = os.path.join(base_dir, f"{stem}_field.vol").replace("\\", "/")
        save_vol_sol_pair(vol_B, sol_B, field_mesh.ngmesh, gf_B)
        vol2msh(out_msh, vol_B, [
            {"sol": sol_B, "fes": "H1", "fes_order": 1, "fes_dim": 3,
             "name": "B", "ncomp": 3},
        ])
        coil_field_msh = out_msh
        progress("VIZ", f"wrote {os.path.basename(out_msh)} "
                         f"(B field on {field_n_vertices} verts, "
                         f"{time.perf_counter()-t0:.1f}s)")
    except Exception as e:
        progress("VIZ", f"BBOX_ERROR:{type(e).__name__}: {e}")

    # ---- 1D filament polylines + |I| ---------------------------------
    filament_msh = ""
    try:
        from radia.gmsh_post_export import export_filaments_msh
        fil_msh_path = os.path.join(base_dir,
                                     f"{stem}_filaments.msh").replace("\\", "/")
        export_filaments_msh(paths, fil_msh_path, currents=I_fil,
                              label="filament")
        filament_msh = fil_msh_path
        progress("VIZ", f"wrote {os.path.basename(fil_msh_path)} "
                         f"({len(paths)} filaments + |I|)")
    except Exception as e:
        progress("VIZ", f"FIL_ERROR:{type(e).__name__}: {e}")

    return {
        "coil_field_msh": coil_field_msh,
        "filament_msh": filament_msh,
        "field_n_vertices": int(field_n_vertices),
    }


# ======================================================================
# Top-level dispatch
# ======================================================================
def run_inductance(args):
    """Top-level orchestrator.  Returns dict; calc_main wraps to JSON."""
    import radia  # noqa: F401  (DLL path setup for subprocess context)
    from ngsolve import TaskManager

    # A sign typo in --frequency must not silently select the DC
    # vacuum-L branch (R=0, Zs ignored) -- fail fast with the flag name.
    if args.frequency < 0.0 or args.frequency != args.frequency:
        return {"status": "error",
                "error": f"--frequency must be >= 0 Hz (0 selects the "
                          f"DC vacuum-L solve); got {args.frequency!r}"}

    # Validate workpiece args if --vol given
    if args.vol:
        if args.sigma <= 0.0:
            return {"status": "error",
                    "error": "--sigma must be > 0 when --vol is provided"}
    if args.coil_only and args.vol:
        progress("INDUCT", "--coil-only: skipping workpiece despite --vol")

    # Coil-input format must match coil-solver:
    #   peec  -> --coil-step (CAD geometry, filaments derived)
    #   bem-a -> --coil-vol  (pre-meshed surface, no on-the-fly OCC re-mesh)
    if args.coil_solver == "peec" and not args.coil_step:
        return {"status": "error",
                "error": "--coil-step is required for --coil-solver peec"}
    if args.coil_solver == "bem-a" and not args.coil_vol:
        return {"status": "error",
                "error": "--coil-vol is required for --coil-solver bem-a"}

    # Strong coupling (iterative CoupledBEMSolver) needs a coil SURFACE mesh
    # to drive the coil EFIE and a workpiece to couple to -- fail fast on a
    # bad contract before the expensive coil solve.
    if args.coupling_mode == "strong":
        if args.coil_solver not in ("bem-a", "peec"):
            return {"status": "error",
                    "error": f"--coupling-mode strong supports --coil-solver "
                             f"bem-a (CoupledBEMSolver, EFIE coil) or peec "
                             f"(CoupledPEECBEMSolver, filament coil); got "
                             f"{args.coil_solver!r}."}
        if args.coil_only or not args.vol:
            return {"status": "error",
                    "error": "--coupling-mode strong requires a workpiece "
                             "--vol (there is nothing to couple to)."}
        if args.coil_solver == "peec" and args.peec_proximity:
            return {
                "status": "error",
                "error": "--coupling-mode strong with --coil-solver peec "
                         "currently uses the isolated-wire Bessel bundle and "
                         "cannot be combined consistently with "
                         "--peec-proximity.  Pass --no-peec-proximity; "
                         "proximity-aware strong coupling remains unsupported.",
            }

    # Loop-DOF extension: an EXPLICIT "on" fails fast on unsupported
    # combinations BEFORE the expensive coil solve (the genus check itself
    # needs the wp mesh and lives in the weak workpiece path).  The
    # default "auto" never errors here -- it resolves inside the weak
    # driver (apply when genus-1 + prerequisites hold, else skip with
    # wp_loop_dof_skip_reason).
    if args.wp_loop_dof == "on":
        if args.coupling_mode != "weak":
            return {"status": "error",
                    "error": "--wp-loop-dof is a weak-coupling workpiece "
                             "extension; strong coupling integration is not "
                             "implemented yet."}
        if args.coil_only or not args.vol:
            return {"status": "error",
                    "error": "--wp-loop-dof requires a workpiece --vol."}
        if args.impedance_model != "sibc":
            return {"status": "error",
                    "error": "--wp-loop-dof supports the linear SIBC only "
                             "(--impedance-model sibc); the ESIM Karl loop "
                             "is not integrated with the loop DOF."}
        if args.wp_bem_backend != "intree-dense":
            return {"status": "error",
                    "error": "--wp-loop-dof needs the dense operators: pass "
                             "--wp-bem-backend intree-dense (the HACApK "
                             "backend does not expose SL/DL for the loop "
                             "column)."}
        if int(args.h1_order) != 1:
            return {"status": "error",
                    "error": "--wp-loop-dof supports the P1 nodal path only "
                             "(--h1-order 1)."}

    # Coil layer.  Wrap the NGSolve work (BEM-A LaplaceSL dense assembly
    # / PEEC) in TaskManager so it runs in parallel: the helpers were
    # de-wrapped under the "caller wraps, helper does NOT" policy
    # (2026-05-27) and this orchestrator is the caller.  Missing this
    # wrap left BEM-A's LaplaceSL single-threaded (kubota 2026-05-29:
    # "BEM-A BEM is slower than before").
    with TaskManager():
        if args.coil_solver == "peec":
            coil_data = _solve_coil_peec(args)
        elif args.coil_solver == "bem-a":
            coil_data = _solve_coil_bem_a(args)
        else:
            return {"status": "error",
                    "error": f"unknown --coil-solver {args.coil_solver!r}"}

    # Coil viz (filament .msh + bbox B .msh).  PEEC only; BEM-A surface
    # current viz is TODO.  Failures here are non-fatal.
    viz = _export_coil_msh_viz(args, coil_data)

    # Vacuum-only branch
    if args.coil_only or not args.vol:
        out = _assemble_vacuum_output(args, coil_data)
        if viz.get("coil_field_msh"):
            out["gmsh_file"] = viz["coil_field_msh"]
            out["coil_field_msh"] = viz["coil_field_msh"]
            out["field_n_vertices"] = viz["field_n_vertices"]
        if viz.get("filament_msh"):
            out["filament_msh"] = viz["filament_msh"]
        return out

    # Coupled workpiece branch (workpiece BEM-SIBC; parallel under
    # TaskManager).  --coupling-mode strong routes to the iterative
    # CoupledBEMSolver (per-DOF back-reaction: the coil current is
    # recomputed against the workpiece reaction field); weak keeps the
    # one-way Telegen ΔL path.
    with TaskManager():
        if args.coupling_mode == "strong":
            if args.coil_solver == "peec":
                strong = _solve_workpiece_strong_coupled_peec(args, coil_data)
            else:
                strong = _solve_workpiece_strong_coupled(args)
            out = _assemble_strong_output(args, coil_data, strong)
        else:
            wp_data = _solve_workpiece_weak_coupled(args, coil_data)
            out = _assemble_full_output(args, coil_data, wp_data)
    if viz.get("coil_field_msh"):
        # Promote bbox B msh to gmsh_file -- it has volumetric B vectors
        # (Air-volume B viz the user can merge interactively).  The wp
        # surface .msh remains in `msh_file` for compatibility.
        out["gmsh_file"] = viz["coil_field_msh"]
        out["coil_field_msh"] = viz["coil_field_msh"]
        out["field_n_vertices"] = viz["field_n_vertices"]
    if viz.get("filament_msh"):
        out["filament_msh"] = viz["filament_msh"]
    return out


def build_argparser():
    """argparse factory shared by main() and notebook DesignSpec callers."""
    parser = argparse.ArgumentParser(
        description="Coil inductance (vacuum + weak-coupled BEM-SIBC). "
                    "Replaces calc_peec_inductance / calc_peec_bem / "
                    "calc_coil_bem_a_workpiece.")

    # ----- Coil geometry & solver switch -----
    # PEEC consumes a CAD .step (filaments derived from centerline);
    # BEM-A consumes a pre-meshed surface .vol (no on-the-fly OCC re-mesh).
    # Each solver requires its own input format -- no cross-feeding.
    parser.add_argument("--coil-step", default="",
                        help="Coil CAD STEP file (required for --coil-solver peec)")
    parser.add_argument("--coil-vol", default="",
                        help="Coil pre-meshed surface .vol (required for "
                             "--coil-solver bem-a; must contain "
                             "'source'/'sink' boundary labels)")
    parser.add_argument("--coil-solver", required=True,
                        choices=["peec", "bem-a"],
                        help="peec: PEEC perimeter filaments from STEP (fast). "
                             "bem-a: Weggler EFIE on RWG, reads pre-meshed .vol.")

    # ----- PEEC-specific args -----
    parser.add_argument("--peec-n-peri", type=int, default=16,
                        help="PEEC perimeter filament count (peec only)")
    parser.add_argument("--peec-proximity", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Use proximity-aware iterative PEEC (default ON; "
                             "augments the Bessel self-skin Zs_fil with the "
                             "Leontovich surface dissipation evaluated from "
                             "the actual Biot-Savart H field at each "
                             "filament's wire-surface position).  Pass "
                             "--no-peec-proximity to fall back to the "
                             "isolated-wire Bessel self-skin only.")
    parser.add_argument("--peec-proximity-max-iter", type=int, default=60,
                        help="Maximum proximity outer iterations.  On the "
                             "3-turn pancake at 150 kHz the iteration "
                             "exhibits a transient hump near iter 30 then "
                             "settles by iter 50 (relax=0.3).  Default 60 "
                             "gives ~10 iter margin.")
    parser.add_argument("--peec-proximity-tol", type=float, default=1.0e-3,
                        help="Proximity outer-iteration convergence tolerance "
                             "on ||dZ|| / ||Z||.")
    parser.add_argument("--peec-proximity-relax", type=float, default=0.3,
                        help="Proximity outer-iteration under-relaxation "
                             "factor (0 < relax <= 1; 0.3 is conservative).")

    # ----- BEM-A-specific args -----
    parser.add_argument("--coil-maxh", type=float, default=0.012,
                        help="DEPRECATED (kept for back-compat); BEM-A "
                             "now reads pre-meshed --coil-vol instead of "
                             "re-meshing CAD on-the-fly.")
    parser.add_argument("--coil-source-name", default="source",
                        help="Coil source-cap face name (bem-a; falls back "
                             "to area-based detection if absent)")
    parser.add_argument("--coil-sink-name", default="sink",
                        help="Coil sink-cap face name (bem-a; same fallback)")
    parser.add_argument("--coil-rwg-quad-degree", type=int, default=5,
                        help="RWG regular-pair quadrature degree (bem-a)")
    parser.add_argument("--coil-rwg-singular-nq", type=int, default=6,
                        help="RWG Sauter-Schwab Duffy n_q (bem-a)")
    parser.add_argument("--coil-bem-solver", default="auto",
                        choices=["dense-lu", "hacapk-gmres", "auto"],
                        help="dense-lu (current) / hacapk-gmres (Phase C.7 future)")
    parser.add_argument("--coil-aca-eps", type=float, default=1e-10,
                        help="HACApK ACA eps for coil (bem-a, future)")
    parser.add_argument("--coil-gmres-tol", type=float, default=1e-10,
                        help="HACApK GMRES tol for coil (bem-a, future)")
    parser.add_argument("--coil-saddle-solver", default="auto",
                        choices=["auto", "lu", "minres", "gmres",
                                 "cocr", "hacapk_cocr"],
                        help="Saddle-point solver for the BEM-A "
                             "impedance-EFIE: auto (LU below 5000 tris for "
                             "AC / 10000 for frequency=0; above that cocr), "
                             "lu (dense LAPACK, overwrite_a=True), cocr "
                             "(SCALABLE: div-free loop reduction + complex-"
                             "symmetric COCR, ~24 mesh-independent iters, "
                             "exact vs LU, dense SL matvec -- the recommended "
                             "large-N solver), hacapk_cocr (the same COCR "
                             "with the HACApKBEMManager-compressed O(N log N) "
                             "SL matvec, accuracy ~3e-7, RT0 only; identical "
                             "R/L), gmres (complex-capable Krylov -- "
                             "unpreconditioned, stalls on large saddles), "
                             "minres (REAL symmetric Krylov -- VALID ONLY for "
                             "the frequency=0 vacuum-L solve; scipy MINRES "
                             "assumes a Hermitian operator and silently "
                             "returns a wrong solution on the complex-"
                             "symmetric AC saddle, so it is rejected there). "
                             "Defaults to auto.")
    # NOTE (2026-07-02): the --bema-impedance-efie flag was REMOVED --
    # BEM-A is unified on the impedance-EFIE (Zs inside the saddle,
    # J = finite-impedance current).  The legacy post-hoc PEC integral
    # R=Re(Zs) integral|J|^2 dS over-estimated R ~3x on tightly-wound /
    # near-contact coils (kubota 3-turn: 15.1 mOhm vs the physical 4.63)
    # and was deleted per the Discard-the-PoC / No-Fallbacks policies.
    # Recover from git history if ever needed for a comparison study.

    # ----- Coil material -----
    parser.add_argument("--coil-sigma", type=float, default=5.8e7,
                        help="Coil conductivity [S/m] (Cu default; mu_r=1)")

    # ----- Workpiece (optional) -----
    parser.add_argument("--vol", default="",
                        help="Workpiece .vol file (omit for vacuum L only)")
    parser.add_argument("--wp-label", default="sibc",
                        help="Workpiece label (volume material OR sideset)")
    parser.add_argument("--sigma", type=float, default=0.0,
                        help="Workpiece conductivity [S/m] (REQUIRED if --vol)")
    parser.add_argument("--mu-r", type=float, default=1.0,
                        help="Workpiece relative permeability")
    parser.add_argument("--half-thickness", type=float, default=0.01,
                        help="ESIM half-thickness [m] (unused in linear)")
    parser.add_argument("--bh-file", default="", help="ESIM BH table (WIP)")

    # ----- Workpiece BEM-SIBC solver -----
    parser.add_argument("--impedance-model", default="sibc",
                        choices=["sibc", "esim"],
                        help="sibc: linear Dowell tanh formula. "
                             "esim: nonlinear Karl iteration; requires "
                             "--bh-file.")
    parser.add_argument("--esim-max-iter", type=int, default=15)
    parser.add_argument("--esim-tol", type=float, default=1e-3)
    parser.add_argument("--esim-relax", type=float, default=0.5,
                        help="Karl under-relaxation (0,1]; 1.0 = full "
                             "step, 0.5 = half-step (default, matches "
                             "calc_fem_kelvin).  Lower if Karl "
                             "oscillates near saturation.")
    parser.add_argument("--esim-per-panel", action="store_true",
                        help="Per-DOF ESIM Karl iteration (Phase B): the "
                             "ESIM cell solver is called once per BEM "
                             "DOF using a per-node H_t extracted from "
                             "the phi_vec, building a per-node Z_s "
                             "ndarray instead of a single scalar.  "
                             "Resolves the spatial saturation pattern; "
                             "costs N_DOF extra ESIM calls per Karl "
                             "iter (~5 ms each).  Requires "
                             "--wp-bem-backend intree-dense.")
    parser.add_argument("--esim-anderson-m", type=int, default=0,
                        help="Anderson-acceleration memory depth for "
                             "the outer Karl iteration.  0 (default) = "
                             "plain damped Picard.  3-5 closes the "
                             "per-DOF dZ_max noise floor when --esim-"
                             "per-panel is used.  See "
                             "src/radia/esim_anderson.py.")
    parser.add_argument("--h1-order", type=int, default=1,
                        choices=[1, 2],
                        help="Workpiece BEM Lagrange basis order "
                             "(NOT the geometry curving order, which is "
                             "auto-detected from the .vol's companion "
                             ".vol.json -- the .vol's curvedelements is "
                             "fixed at Cubit-export time and a post-load "
                             "mesh.Curve(p) silently falls back to flat). "
                             "1 = P1 hat (flat-friendly default).  "
                             "2 = Lagrange-P2 (uses curved Tri6 geometry "
                             "if the .vol has curve order >= 2, otherwise "
                             "P2 basis on flat Tri3 geometry).")
    parser.add_argument("--wp-bem-backend", default="hacapk",
                        choices=["hacapk", "intree-dense"])
    parser.add_argument("--wp-aca-eps", type=float, default=1e-10)
    parser.add_argument("--wp-gmres-tol", type=float, default=1e-10)
    parser.add_argument("--wp-loop-dof", nargs="?", const="on",
                        default="auto", choices=["auto", "on"],
                        help="Genus-1 loop DOF (net shorted-turn eddy "
                             "current) for the workpiece scalar BIE "
                             "(radia.bem_loop_extension).  Recovers the "
                             "Lenz screening a single-valued potential "
                             "cannot represent on a ring/tube workpiece "
                             "whose bore links the coil flux.  "
                             "auto (default): apply on a genus-1 workpiece "
                             "when the prerequisites hold (weak coupling, "
                             "linear SIBC, --wp-bem-backend intree-dense, "
                             "--h1-order 1); otherwise skip and report "
                             "wp_loop_dof_skip_reason (genus-1 skips also "
                             "keep the P_wp_caveat).  "
                             "on (= bare --wp-loop-dof): require it -- "
                             "unmet prerequisites or genus != 1 fail "
                             "loud.  There is deliberately NO 'off': the "
                             "legacy un-extended genus-1 solve is a known "
                             "+25-30% over-estimate and must not be "
                             "selectable.  delta_L (Telegen) keeps the "
                             "plain-solve phi convention.")

    # ----- Excitation -----
    parser.add_argument("--frequency", type=float, required=True,
                        help="Frequency [Hz]")
    parser.add_argument("--current", type=float, default=1.0,
                        help="Port current [A]")

    # ----- Coupling -----
    parser.add_argument("--coupling-mode", default="weak",
                        choices=["weak", "strong"],
                        help="weak: one-way coil->wp with back-reacted Telegen "
                             "ΔL (coil current NOT recomputed).  strong: "
                             "iterative CoupledBEMSolver (per-DOF back-reaction; "
                             "coil current recomputed against the workpiece "
                             "reaction field, so ΔL includes the wp magnetic-"
                             "energy term and P_wp is self-consistent.  Supports "
                             "bem-a, or experimental peec with "
                             "--no-peec-proximity).")
    parser.add_argument("--coupling-max-iter", type=int, default=10,
                        help="strong: max Picard iterations for the coupled "
                             "coil<->wp back-reaction solve.")
    parser.add_argument("--coupling-tol", type=float, default=1e-3,
                        help="strong: relative coupled-state convergence "
                             "tolerance for the Picard loop (L_total for BEM-A, "
                             "terminal impedance for PEEC).")
    parser.add_argument("--coupling-relax", type=float, default=0.5,
                        help="strong: under-relaxation factor (0<relax<=1) on "
                             "the coupled back-reaction RHS.")
    parser.add_argument("--telegen-form", default="phi-B",
                        choices=["phi-B"],
                        help="Telegen form (only phi.B currently)")

    # ----- Output -----
    parser.add_argument("--msh-output", default="",
                        help="Workpiece .msh viz output path")
    parser.add_argument("--coil-msh-output", default="",
                        help="Coil .msh viz output path. PEEC: writes "
                             "<this>=B-on-bbox + <stem>_filaments.msh "
                             "(1D polylines + |I|).  open_gmsh.py auto-"
                             "merges the filament sibling.")
    parser.add_argument("--write-summary", action="store_true",
                        help="Add extra diagnostics in JSON")
    parser.add_argument("--coil-only", action="store_true",
                        help="Skip workpiece even if --vol given")

    # ----- Performance -----
    parser.add_argument("--n-threads", type=int, default=0,
                        help="0 = auto (TaskManager)")
    return parser


def main():
    parser = build_argparser()
    return calc_main(run_inductance, parser)


if __name__ == "__main__":
    sys.exit(main())
