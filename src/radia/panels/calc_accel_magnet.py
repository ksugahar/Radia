"""
Headless accelerator-magnet solver for the Radia Electromagnet application.

Two formulations:
  omega  -- Omega-reduced scalar potential (H1, fastest)
  a      -- A-formulation vector potential (HCurl, eddy-current ready)

User workflow:
  1. Write coil Python script (defines build_coil() -> CoilBuilder)
  2. Write .jou to create iron yoke + air + mesh + named blocks
  3. Cubit's embedded toolbar exports and validates the .vol
  4. The Simulink block or MCP runner invokes this headless solver

Mesh interface:
  This script takes a Netgen .vol file via --vol (not .cub5).  The .vol
  is the sole interface between Cubit and NGSolve per the Cubit/NGSolve
  Complete Separation Policy in CLAUDE.md — this script does NOT
  import cubit.  Mesh order + element types are embedded in the .vol
  by the Cubit plugin at export time.

Mesh materials / boundaries (named in the .jou before export):
  yoke            - volume material, nonlinear iron (BH curve)
  air             - volume material
  sym_bn=0_<x|y|z>  - symmetry plane where B . n = 0 (field parallel)
  sym_ht=0_<x|y|z>  - symmetry plane where H x n = 0 (field perpendicular)
  sym_tangential   - LEGACY alias for sym_bn=0 (pre-2026-04-25)
  sym_normal       - LEGACY alias for sym_ht=0 (pre-2026-04-25)
  kelvin          - volume material (optional, Periodic Kelvin)
  kelvin_far      - 1/8-reduction "infinity plane" (always Dirichlet)

Symmetry BC (auto-swapped by formulation - the label names the physics,
the solver picks Dirichlet vs Natural based on whether it is A or Omega):
  sym_bn=0_*:  B . n = 0 -> Omega: natural,  A: Dirichlet (A x n = 0)
  sym_ht=0_*:  H x n = 0 -> Omega: Dirichlet (Omega=const), A: natural
  kelvin_far:  always Dirichlet (Omega=0 / A=0) -- represents the same
               infinity plane as the GND vertex, extended to a face.

Created by the Cubit mesh-export workflow when Kelvin geometry is requested:
  kelvin_int  - surface, periodic BC (interior hemisphere)
  kelvin_ext  - surface, periodic BC (exterior hemisphere)

Coil is NOT meshed -- Biot-Savart source from CoilBuilder wire path.
"""

import argparse
import importlib.util
import json
import math
import os
import sys
import time

# Pre-import scipy BEFORE Cubit (Cubit bundles broken scipy)
try:
    import scipy  # noqa: F401
    import scipy.interpolate  # noqa: F401
except ImportError:
    pass

import numpy as np

# Shared utilities
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (MU_0, NU_0, setup_paths,
                          add_periodic_kelvin, detect_kelvin_offset,
                          get_bh_curve, progress, calc_main,
                          EMMaterial, add_material_args, apply_fe_inverse)
from radia.em_design import FEM_SOLVERS


def _log(msg):
    """Write progress to stderr (panel reads these)."""
    progress("ACCEL", msg)


def _select_accel_solver(solver, formulation, order, periodic, ndof):
    if solver not in FEM_SOLVERS:
        raise ValueError(f"Unsupported FE solver: {solver}")
    if solver == "auto":
        if formulation == "a" and not periodic and order in (1, 2, 3):
            return "bddc_ams"
        return "sparsecholesky" if ndof < 200000 else "bddc"
    if solver == "bddc_ams":
        if formulation != "a" or periodic or order not in (1, 2, 3):
            raise ValueError("AMS requires nonperiodic HCurl order 1, 2 or 3")
        # Preserve the existing full HCurl space, including gradient modes.
        # Its lowest-order wirebasket uses AMS; BDDC treats the other modes.
        return "bddc_ams"
    return solver


def _accel_inverse(a, fes, solver):
    """Construct the requested inverse without changing the FE operator."""
    from ngsolve import Preconditioner, TaskManager
    if solver == "sparsecholesky":
        with TaskManager():
            return a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    if solver == "iccg":
        from radia.sparsesolv_ngsolve import SparseSolvSolver
        return SparseSolvSolver(a.mat, method="ICCG", freedofs=fes.FreeDofs(),
                                tol=1e-12, maxiter=2000, shift=1.05,
                                save_best_result=True, printrates=False)
    if solver not in ("bddc", "bddc_ams"):
        raise ValueError(f"Unsupported resolved FE solver: {solver}")
    flags = {"inverse": "sparsecholesky"}
    if solver == "bddc_ams":
        import radia.sparsesolv_ngsolve  # registers the coarse solver
        flags["coarsetype"] = "sparsesolv_ams"
    pre = Preconditioner(a, "bddc", **flags)
    with TaskManager():
        a.Assemble()
    from ngsolve.krylovspace import CGSolver
    return CGSolver(a.mat, pre.mat, maxiter=2000, tol=1e-10, printrates=False)


def _elem_mat(mesh, materials):
    """Build element-nr -> material-name mapping.

    mesh.GetMaterials() returns domain-indexed names (short tuple).
    Each element has el.index -> domain index.
    L2 order=0 DOFs are indexed by el.nr.

    Returns:
        dict: {el.nr: material_name} for all VOL elements
    """
    from ngsolve import VOL
    return {el.nr: materials[el.index]
            for el in mesh.Elements(VOL)}


def _load_coil_script(script_path):
    """Load coil script and call build_coil() -> CoilBuilder.

    The script must define build_coil() returning a CoilBuilder instance.
    """
    if not script_path or not os.path.exists(script_path):
        raise FileNotFoundError(f"Coil script not found: {script_path}")

    setup_paths()

    # EM panel contract (2026-04-25): analytical coils are authored
    # as a Python module exposing `build_coil() -> CoilBuilder`.
    # STEP input is intentionally NOT supported here; it belongs to
    # the PEEC side (IH panel -> coil_from_cad.filaments_from_step).
    # The separation is explicit: .py for Biot-Savart analytical,
    # .step for PEEC filament circuits.
    ext = os.path.splitext(script_path)[1].lower()
    if ext in ('.step', '.stp'):
        raise ValueError(
            f"EM panel does not accept STEP coil inputs.  Use a Python "
            f"module that defines `build_coil() -> CoilBuilder`.  "
            f"PEEC workflows (IH panel) accept STEP via peec_step.  "
            f"Got: {script_path}")

    script_dir = os.path.dirname(os.path.abspath(script_path))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    spec = importlib.util.spec_from_file_location("coil_module", script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    if not hasattr(mod, 'build_coil'):
        raise AttributeError(
            f"Coil script must define build_coil() function: {script_path}")

    coil = mod.build_coil()
    return coil


def _create_bh_interpolators(bh_data):
    """Create chord and differential interpolators from BH curve data.

    The B(H) law is the production one shared by every Radia route
    (:func:`radia.bh_law.monotone_bh_pchip`) inside the table, continued at
    the vacuum slope ``B_max + mu0 (H - H_max)`` beyond it, as the static
    electromagnet solver does.

    Args:
        bh_data: list of [H, B] pairs (H in A/m, B in Tesla), starting at
            [0, 0] with strictly increasing H and non-decreasing B.

    Returns:
        dict with keys:
            mu_chord(H) -> B/H  (chord permeability)
            nu_chord(B) -> H/B  (chord reluctivity)
            mu_diff(H)  -> dB/dH (differential permeability, for Newton)
            nu_diff(B)  -> dH/dB (differential reluctivity, for Newton)
            mu_r_init   -> initial mu_r
    """
    from scipy.optimize import brentq
    from radia.bh_law import monotone_bh_pchip

    bh = np.asarray(bh_data, dtype=float)
    if bh.ndim != 2 or bh.shape[1] < 2 or bh.shape[0] < 2:
        raise ValueError("bh_data must contain at least two [H, B] rows")
    H_tab = bh[:, 0]
    B_tab = bh[:, 1]
    if not (np.all(np.isfinite(H_tab)) and np.all(np.isfinite(B_tab))):
        raise ValueError("bh_data must contain finite H and B values")
    if np.any(np.diff(H_tab) <= 0.0):
        raise ValueError("bh_data H values must be strictly increasing")
    if np.any(np.diff(B_tab) < 0.0):
        raise ValueError("bh_data B values must be non-decreasing")
    if H_tab[0] != 0.0 or B_tab[0] != 0.0:
        raise ValueError("bh_data must start at [0, 0]")

    law = monotone_bh_pchip(H_tab, B_tab, extrapolate=False)
    slope = law.derivative()
    H_max = float(H_tab[-1])
    B_max = float(B_tab[-1])
    # The chord reluctivity needs H(B), so the law must be strictly
    # increasing: dB/dH > 0 on every piece.  dB/dH is a quadratic per piece,
    # so its minimum is exact (endpoints or the vertex); nothing is floored.
    widths = np.diff(H_tab)
    for piece in range(len(widths)):
        a, b, c = slope.c[:, piece]
        candidates = [0.0, float(widths[piece])]
        if a > 0.0:
            vertex = -b / (2.0 * a)
            if 0.0 < vertex < widths[piece]:
                candidates.append(float(vertex))
        minimum = min(a * t * t + b * t + c for t in candidates)
        if not minimum > 0.0:
            raise ValueError(
                "bh_data gives a B(H) law that is not strictly increasing on "
                f"[{H_tab[piece]:g}, {H_tab[piece + 1]:g}] A/m (min dB/dH = "
                f"{minimum:.3e}); H(B) and the chord reluctivity are undefined")

    def b_of_h(H):
        H = abs(float(H))
        if H >= H_max:
            return B_max + MU_0 * (H - H_max)
        return float(law(H))

    def db_dh(H):
        H = abs(float(H))
        if H >= H_max:
            return MU_0
        return float(slope(H))

    def h_of_b(B):
        B = abs(float(B))
        if B >= B_max:
            return H_max + (B - B_max) / MU_0
        if B <= 0.0:
            return 0.0
        # B(H) is strictly increasing on the table, so the root is bracketed
        # and unique.  Stop on relative accuracy only, so a small nonzero B
        # keeps its own nonzero H instead of an absolute-tolerance zero.
        return float(brentq(lambda H: float(law(H)) - B, 0.0, H_max,
                            xtol=1e-300, rtol=4.0 * np.finfo(float).eps,
                            maxiter=500))

    mu_initial = db_dh(0.0)

    def mu_chord(H):
        H = abs(float(H))
        return mu_initial if H == 0.0 else b_of_h(H) / H

    def nu_chord(B):
        B = abs(float(B))
        return 1.0 / mu_initial if B == 0.0 else h_of_b(B) / B

    def nu_diff(B):
        return 1.0 / db_dh(h_of_b(B))

    return {
        'mu_chord': np.vectorize(mu_chord, otypes=[float]),
        'nu_chord': np.vectorize(nu_chord, otypes=[float]),
        'mu_diff': np.vectorize(db_dh, otypes=[float]),
        'nu_diff': np.vectorize(nu_diff, otypes=[float]),
        'mu_r_init': float(mu_initial / MU_0),
    }


def solve_accel(coil_script="", vol_file="", formulation="omega",
                fes_order=1,
                mat=None, n_steps=1,
                max_iter=30, tol=1e-3, relax=0.3,
                newton=False, solver="auto", msh_output=""):
    """Accelerator magnet solver (Omega-reduced or A-formulation).

    Args:
        coil_script: Path to Python script with build_coil()
        vol_file: Netgen .vol file (exported by Cubit `export netgen`
            or produced directly by Netgen; the .vol is the sole interface
            between Cubit and NGSolve per the Cubit/NGSolve separation
            policy — this script no longer imports cubit).
        formulation: "omega" (scalar H1) or "a" (vector HCurl)
        fes_order: FE polynomial order
        mat: EMMaterial instance (iron yoke properties)
        n_steps: Number of quasi-static steps (hysteresis: 0->I->0)
        max_iter: Max Picard/Hantila iterations
        tol: Convergence tolerance
        relax: Fraction of the new material step taken per iteration,
            in (0, 1]: 1 is the full (undamped) step, smaller damps more.
        msh_output: Optional GMSH .msh output path

    Returns:
        dict with B_center, L, W_mag, ndof, iterations, etc.
    """
    if solver not in FEM_SOLVERS:
        raise ValueError(f"Unsupported FE solver: {solver}")
    if max_iter < 1:
        raise ValueError("max_iter must be positive")
    if not 0.0 < float(relax) <= 1.0:
        # relax = 0 would never update the material and then report the
        # initial-permeability solution as converged.
        raise ValueError("relax must lie in (0, 1]; 1 is the full material step")
    if mat is None:
        mat = EMMaterial.from_name("steel")
    material = mat.name
    mu_r = mat.mu_r
    bh_file = ""  # BH data now lives in mat.bh_curve
    hys_file = mat.hys_file

    import ngsolve  # noqa: F401
    from ngsolve import (H1, HCurl, L2, Periodic, BilinearForm, LinearForm,
                         GridFunction, Integrate, InnerProduct,
                         curl, grad, dx, CF, VOL,
                         TaskManager, Mesh)
    from ngsolve import x, y, z, sqrt

    setup_paths()
    import radia as rad

    t_total_start = time.perf_counter()

    # ============================================================
    # Step 1: Load coil -> Radia ObjRecCur/ObjArcCur (exact analytical)
    # ============================================================
    _log("COIL:loading script")
    coil = _load_coil_script(coil_script)
    current = coil.current
    segments, _ = coil.to_wire_segments(n_arc=20)  # for segment count info

    # Create Radia coil objects (ObjRecCur, ObjArcCur) with exact formulas
    radia_coil_objs = coil.to_radia()
    coil_container = rad.ObjCnt(radia_coil_objs)
    _log(f"COIL:{len(radia_coil_objs)} Radia objects, I={current}A")

    # RadiaField CoefficientFunction (exact analytical: BufVect 8-corner formula)
    # Returns B in Tesla directly (not H)
    B_s = rad.RadiaField(coil_container, 'b')
    H_s = rad.RadiaField(coil_container, 'h')
    _log("COIL:RadiaField CF created (exact analytical)")

    # ============================================================
    # Step 2: Load .vol mesh (Cubit-free; .vol is the sole interface)
    # ============================================================
    _log("MESH:loading")
    t0 = time.perf_counter()

    if not vol_file or not os.path.exists(vol_file):
        return {"error": f".vol file required: {vol_file}"}

    mesh = Mesh(vol_file)

    # Kelvin: detect offset, estimate radius, add periodic identification
    a_kelvin = 0.060
    kelvin_center = detect_kelvin_offset(mesh)
    try:
        air_verts = set()
        mats = mesh.GetMaterials()
        for el in mesh.Elements(VOL):
            if "air" in str(mats[el.index]):
                for v in el.vertices:
                    air_verts.add(v.nr)
        if air_verts:
            coords = np.array([mesh.vertices[v].point for v in air_verts])
            a_kelvin = float(np.max(np.linalg.norm(coords, axis=1)))
    except Exception:
        pass
    if add_periodic_kelvin(mesh, kelvin_center):
        mesh = ngsolve.Mesh(mesh.ngmesh)
        _log(f"PERIODIC:added (a={a_kelvin:.4f}, offset={kelvin_center})")

    # The Kelvin transformation needs a curved Kelvin sphere: the
    # (R/rho')^2 reluctivity scaling assumes a smooth sphere, and a
    # polyhedral one produces large errors outside the shell (validated
    # 2026-04-17: 1.15% vs 2D axisym at geometry order 2).  Loading the
    # .vol already applies its stored curvedelements; Mesh.Curve on a
    # loaded mesh would flatten a CAD-less Cubit export or borrow another
    # geometry, so only curve through the guarded helper.
    if fes_order >= 2:
        from radia.mesh_curve import CurveOrderError, ensure_curve_order
        try:
            curve = ensure_curve_order(mesh, fes_order, what=str(vol_file))
        except CurveOrderError as exc:
            return {"error": str(exc)}
        _log(f"MESH:geometry order {curve['curve_order']} ({curve['action']})")

    t_mesh = time.perf_counter() - t0
    materials = mesh.GetMaterials()  # domain-indexed (short tuple)
    # BND (codim-1, surface) + BBBND (codim-3, vertex).  "GND" for the
    # Kelvin-centre Dirichlet is typically a BBBND vertex label set by
    # `add_kelvin_cubit` via the Cubit nodeset -> Netgen CD3 propagation
    # path in export netgen.  Join both namespaces so the existing
    # Dirichlet-pick logic that searches a single `boundaries` iterable
    # still works.
    bnd_list = list(mesh.GetBoundaries())
    try:
        bnd_list += list(mesh.GetBBBoundaries())
    except Exception:
        pass
    boundaries = bnd_list
    ne = mesh.GetNE(VOL)
    elem_mat = _elem_mat(mesh, materials)  # el.nr -> material name
    _log(f"MESH:done ({ne} elems, {t_mesh:.1f}s)")

    has_kelvin = "kelvin" in materials
    has_yoke = any("yoke" in m.lower() for m in materials)

    # ============================================================
    # Step 3: BH curve / material setup
    # ============================================================
    is_nonlinear = False
    is_hysteresis = False
    bh_interp = None  # dict of interpolators
    mu_r_init = mu_r
    hys_mat_handles = None
    nu_rev = None

    if material == "hysteresis" and hys_file and has_yoke:
        # Hysteresis mode: B-input Play model via Hantila decomposition
        # A-formulation: H = nu_rev * B + H_irr(B), B directly available
        # Omega-reduced:  M = alpha * H + R(H),     H directly available
        # Both: constant LHS factored once, only RHS updated per iteration

        from hysteresis_io import hys_to_play_radia

        K_hys, eta, f_k_tables = hys_to_play_radia(hys_file)
        _log(f"HYS:Play model K={K_hys}, file={os.path.basename(hys_file)}")

        # Create per-element material handles (each has independent state)
        yoke_elem_indices = [nr for nr, m in elem_mat.items()
                             if "yoke" in m.lower()]
        n_yoke = len(yoke_elem_indices)
        hys_mat_handles = []
        for _ in range(n_yoke):
            h = rad.MatPlayHysteresis(K_hys, eta, f_k_tables)
            hys_mat_handles.append(h)
        _log(f"HYS:{n_yoke} per-element Play materials created")

        # Get nu_rev from first handle (same for all - same material)
        nu_rev = rad.MatHysGetNuRev(hys_mat_handles[0])
        _log(f"HYS:nu_rev={nu_rev:.4e}")

        is_hysteresis = True
        is_nonlinear = True

    elif has_yoke and material not in ("linear", "hysteresis"):
        bh_data, _ = mat.get_bh_curve()
        if bh_data:
            bh_interp = _create_bh_interpolators(bh_data)
            mu_r_init = bh_interp['mu_r_init']
            is_nonlinear = True
            _log(f"BH:nonlinear, initial mu_r={mu_r_init:.0f}, "
                 f"method={'Newton' if newton else 'Picard'}")
        else:
            _log(f"BH:linear mu_r={mu_r}")
    else:
        _log(f"MATERIAL:linear mu_r={mu_r}")

    # ============================================================
    # Step 4: Kelvin weight CF
    # ============================================================
    kx, ky, kz = kelvin_center
    dx_k = x - kx
    dy_k = y - ky
    dz_k = z - kz
    rp_sq = dx_k * dx_k + dy_k * dy_k + dz_k * dz_k + 1e-20
    kelvin_fac = a_kelvin**2 / rp_sq  # (a/r')^2

    # Build kelvin_weight: 1 everywhere, (a/r')^2 in kelvin domain
    kw_dict = {}
    for m in materials:
        if "kelvin" in m.lower():
            kw_dict[m] = kelvin_fac
        else:
            kw_dict[m] = CF(1.0)
    kelvin_weight = mesh.MaterialCF(kw_dict, default=CF(1.0))

    # ============================================================
    # Step 5: FE space
    # ============================================================
    # Symmetry BC: formulation-dependent Dirichlet assignment.
    # The mesh carries physics-named labels; the solver picks
    # Dirichlet vs Natural based on whether it is A or Omega.
    #   sym_bn=0_*  (B.n=0, field parallel)  -> Omega: natural,  A: Dirichlet
    #   sym_ht=0_*  (Hxn=0, field perp)      -> Omega: Dirichlet, A: natural
    # Legacy labels (pre-2026-04-25) "sym_tangential"/"sym_normal" are
    # also recognised and behave as sym_bn=0/sym_ht=0 respectively.
    # `kelvin_far` (1/8-reduction only) is the "infinity plane"
    # through the Kelvin centre and is ALWAYS Dirichlet (Omega=0 / A=0)
    # regardless of formulation -- it represents the same boundary as
    # the GND vertex, just extended to a flat face.
    dir_parts = []
    if "GND" in boundaries:
        dir_parts.append("GND")
    if "outer" in boundaries:
        dir_parts.append("outer")
    if "kelvin_far" in boundaries:
        dir_parts.append("kelvin_far")

    # Pick Dirichlet per formulation.
    if formulation == "omega":
        dir_bc = "ht=0"       # Omega = const where H x n = 0
        legacy_dir = "sym_normal"
    else:
        dir_bc = "bn=0"       # A x n = 0 where B . n = 0
        legacy_dir = "sym_tangential"

    # Canonical labels: sym_<bc>_<axis> for axis in x/y/z.
    for axis in ("x", "y", "z"):
        label = f"sym_{dir_bc}_{axis}"
        if label in boundaries:
            dir_parts.append(label)

    # Legacy label (axis-less) back-compat.
    if legacy_dir in boundaries:
        dir_parts.append(legacy_dir)

    dirichlet_bnd = "|".join(dir_parts) if dir_parts else ""
    _log(f"BC:dirichlet={dirichlet_bnd or '(none)'}")

    # Omega-reduced (H1) on Kelvin requires a Dirichlet pin at the Kelvin
    # center (= physical infinity in the transformed coordinates) for
    # uniqueness. The .jou is contractually required to define this as a
    # nodeset named "GND" on a vertex (or sideset on a tiny boundary).
    # We do NOT auto-pick the nearest vertex — that result depends on the
    # mesh and changes silently between regenerations.
    if (has_kelvin and "GND" not in boundaries
            and formulation == "omega"):
        return {"error":
                "Omega-reduced formulation with Kelvin transformation "
                "requires a Dirichlet 'GND' label at the Kelvin center "
                "(physical infinity). The .jou must define a nodeset "
                "or sideset named 'GND' at the kelvin_center vertex. "
                "Auto-picking the nearest vertex is disabled because "
                "the choice depends on mesh discretization and would "
                "give silently different results across re-meshes."}

    if formulation == "omega":
        base_fes = H1(mesh, order=fes_order, dirichlet=dirichlet_bnd)
        if has_kelvin:
            fes = Periodic(base_fes)
            _log("FES:Periodic H1 (Omega-reduced)")
        else:
            fes = base_fes
            _log("FES:H1 (Omega-reduced)")
    else:  # "a"
        base_fes = HCurl(mesh, order=fes_order, dirichlet=dirichlet_bnd)
        if has_kelvin:
            fes = Periodic(base_fes)
            _log("FES:Periodic HCurl (A-formulation)")
        else:
            fes = base_fes
            _log("FES:HCurl (A-formulation)")

    u, v = fes.TnT()
    ndof = fes.ndof
    _log(f"FES:ndof={ndof}")

    linear_solver = _select_accel_solver(solver, formulation, fes_order, has_kelvin, ndof)
    _log(f"SOLVER:{linear_solver}")

    # L2 space for per-element material property
    fes_l2 = L2(mesh, order=0)

    # ============================================================
    # Step 6: Nonlinear iteration
    # ============================================================
    gfu = GridFunction(fes)
    history = []
    n_iter_actual = 0
    converged = False
    step_results = []  # for hysteresis quasi-static steps

    if is_hysteresis:
        # ==========================================================
        # Hantila iteration (Play model, both formulations)
        #
        # A-formulation:
        #   H = nu_rev * B + H_irr(B)
        #   LHS: int nu_rev * curl(A) . curl(v) dx  (factored ONCE)
        #   RHS: -int (nu_rev * B_s + H_irr) . curl(v) dx
        #   Update: B = curl(A_r) + B_s → H_irr = MatHysIrreversible(mat, B)
        #
        # Omega-reduced:
        #   M = alpha * H + R(H),  alpha = 1/(mu_0*nu_rev) - 1
        #   LHS: int mu_eff * grad(Omega) . grad(v) dx  (factored ONCE)
        #        where mu_eff = mu_0*(1+alpha) = 1/nu_rev
        #   RHS: int mu_eff * H_s . grad(v) + mu_0 * R . grad(v) dx
        #   Update: H = H_s - grad(Omega) → M = MatMvsH(mat, H), R = M - alpha*H
        # ==========================================================
        yoke_elem_indices = [nr for nr, m in elem_mat.items()
                             if "yoke" in m.lower()]
        alpha_hys = 1.0 / (MU_0 * nu_rev) - 1.0  # reversible susceptibility
        mu_eff = 1.0 / nu_rev  # = mu_0 * (1 + alpha)

        # Constant material CF (Kelvin-weighted)
        gf_const = GridFunction(fes_l2)
        for nr, m in elem_mat.items():
            if "yoke" in m.lower():
                gf_const.vec[nr] = nu_rev if formulation == "a" else mu_eff
            else:
                gf_const.vec[nr] = NU_0 if formulation == "a" else MU_0
        const_cf = gf_const * kelvin_weight

        # Assemble and factor LHS ONCE
        a_bf = BilinearForm(fes)
        if formulation == "a":
            a_bf += const_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
            a_bf += 1e-6 * NU_0 * u * v * dx  # gauge
        else:
            a_bf += const_cf * grad(u) * grad(v) * dx(bonus_intorder=4)
        a_bf.Assemble()

        inv_a = _accel_inverse(a_bf, fes, linear_solver)

        _log(f"HANTILA:LHS factored (nu_rev={nu_rev:.4e}, "
             f"alpha={alpha_hys:.1f}, form={formulation})")

        # Nonlinear residual GridFunction (per-element vector, L2 order=0 dim=3)
        fes_l2_vec = L2(mesh, order=0, dim=3)
        gf_resid = GridFunction(fes_l2_vec)  # H_irr (A-form) or R (Omega)

        # Quasi-static current ramp: 0 -> I_max -> 0
        if n_steps <= 1:
            current_levels = [current]
        else:
            half = n_steps // 2
            ramp_up = [current * (i + 1) / half for i in range(half)]
            ramp_down = [current * (half - 1 - i) / half for i in range(half)]
            current_levels = ramp_up + ramp_down
            if len(current_levels) == 0:
                current_levels = [current]

        for step_idx, I_step in enumerate(current_levels):
            _log(f"STEP:{step_idx+1}/{len(current_levels)} I={I_step:.2f}A")

            # Scale Biot-Savart for this current level
            scale = I_step / current if current != 0 else 0
            B_s_step = scale * B_s
            H_s_step = scale * H_s

            # Save hysteresis state at step start
            saved_states = [rad.MatHysSaveState(h) for h in hys_mat_handles]

            for iteration in range(max_iter):
                t0_iter = time.perf_counter()

                # Build RHS (constant source + nonlinear residual)
                f_lf = LinearForm(fes)
                if formulation == "a":
                    # RHS = -nu_rev * B_s . curl(v) - H_irr . curl(v)
                    f_lf += -const_cf * B_s_step * curl(v) * dx(bonus_intorder=4)
                    f_lf += -gf_resid * curl(v) * dx(bonus_intorder=4)
                else:
                    # RHS = mu_eff * H_s . grad(v) + mu_0 * R . grad(v)
                    f_lf += const_cf * H_s_step * grad(v) * dx(bonus_intorder=4)
                    f_lf += MU_0 * gf_resid * grad(v) * dx(bonus_intorder=4)
                f_lf.Assemble()

                # Back-substitution (no re-factorization)
                with TaskManager():
                    linear_residual = apply_fe_inverse(
                        a_bf.mat, inv_a, f_lf.vec, gfu.vec, fes.FreeDofs())

                t_solve_iter = time.perf_counter() - t0_iter

                # Update nonlinear residual per yoke element
                max_change = 0.0

                for j, elem_idx in enumerate(yoke_elem_indices):
                    try:
                        el = mesh.Elements(VOL)[elem_idx]
                        mp = mesh(*el.vertices[0].point)
                    except Exception:
                        continue

                    if formulation == "a":
                        # B = curl(A_r) + B_s → H_irr from Play Forward
                        B_total_cf = curl(gfu) + B_s_step
                        B_vec = [float(c) for c in B_total_cf(mp)]
                        resid_new = np.array(rad.MatHysIrreversible(
                            hys_mat_handles[j], B_vec), dtype=float)
                    else:
                        # H = H_s - grad(Omega) → M from MatMvsH, R = M - alpha*H
                        H_total_cf = H_s_step - grad(gfu)
                        H_vec = [float(c) for c in H_total_cf(mp)]
                        M_vec = np.array(rad.MatMvsH(
                            hys_mat_handles[j], H_vec), dtype=float)
                        resid_new = M_vec - alpha_hys * np.array(H_vec)

                    # Convergence check
                    resid_old = np.array([
                        float(gf_resid.vec[elem_idx * 3 + c])
                        for c in range(3)])
                    diff = np.linalg.norm(resid_new - resid_old)
                    norm_old = np.linalg.norm(resid_old) + 1e-30
                    max_change = max(max_change, diff / norm_old)

                    # Update residual GridFunction
                    for c in range(3):
                        gf_resid.vec[elem_idx * 3 + c] = float(resid_new[c])

                n_iter_actual = iteration + 1
                history.append({
                    'step': step_idx,
                    'iteration': iteration,
                    'max_change': float(max_change),
                    't_solve': t_solve_iter,
                })
                _log(f"HANTILA:{step_idx+1}/{len(current_levels)} "
                     f"iter={iteration} d_resid={max_change:.4e} "
                     f"t={t_solve_iter:.1f}s")

                if max_change < tol and iteration > 0:
                    _log(f"CONVERGED:step={step_idx+1} iter={iteration}")
                    converged = True
                    break

            # Commit hysteresis state for next step
            for h in hys_mat_handles:
                rad.MatHysCommitState(h)

            # Record step result
            try:
                if formulation == "a":
                    B_cf = curl(gfu) + B_s_step
                else:
                    B_cf = gf_const * (H_s_step - grad(gfu))
                B_step = [float(c) for c in B_cf(mesh(0, 0, 0))]
                B_step_mag = math.sqrt(sum(b**2 for b in B_step))
            except Exception:
                B_step_mag = 0
            step_results.append({
                'step': step_idx, 'current': I_step,
                'B_center': B_step_mag, 'iterations': n_iter_actual,
            })

    else:
        # ==========================================================
        # Picard or Newton iteration (BH curve or linear)
        #
        # Newton tangent for A-formulation:
        #   dH/dB = nu_chord*(I - B_hat x B_hat) + nu_diff * B_hat x B_hat
        # Newton tangent for Omega-reduced:
        #   dB/dH = mu_chord*(I - H_hat x H_hat) + mu_diff * H_hat x H_hat
        # ==========================================================
        yoke_elem_indices = [nr for nr, m in elem_mat.items()
                             if "yoke" in m.lower()]
        gf_prop = GridFunction(fes_l2)       # chord: mu or nu
        gf_prop_diff = GridFunction(fes_l2)  # differential: mu_diff or nu_diff
        for nr, m in elem_mat.items():
            if "yoke" in m.lower():
                if formulation == "omega":
                    gf_prop.vec[nr] = MU_0 * mu_r_init
                    gf_prop_diff.vec[nr] = MU_0 * mu_r_init
                else:
                    gf_prop.vec[nr] = 1.0 / (MU_0 * mu_r_init)
                    gf_prop_diff.vec[nr] = 1.0 / (MU_0 * mu_r_init)
            else:
                if formulation == "omega":
                    gf_prop.vec[nr] = MU_0
                    gf_prop_diff.vec[nr] = MU_0
                else:
                    gf_prop.vec[nr] = NU_0
                    gf_prop_diff.vec[nr] = NU_0

        for iteration in range(max_iter):
            t0_iter = time.perf_counter()

            chord_cf = gf_prop * kelvin_weight

            if newton and is_nonlinear and iteration > 0:
                # Newton: anisotropic tangent matrix
                # dH/dB = chord*(I - B_hat*B_hat^T) + diff * B_hat*B_hat^T
                #       = chord*I + (diff - chord) * B_hat*B_hat^T
                diff_cf = gf_prop_diff * kelvin_weight
                if formulation == "a":
                    B_total = curl(gfu) + B_s
                    B_mag_sq = InnerProduct(B_total, B_total) + 1e-30
                    a_bf = BilinearForm(fes)
                    a_bf += chord_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
                    a_bf += ((diff_cf - chord_cf) / B_mag_sq
                             * InnerProduct(B_total, curl(u))
                             * InnerProduct(B_total, curl(v))
                             ) * dx(bonus_intorder=4)
                    a_bf += 1e-6 * NU_0 * u * v * dx
                else:
                    H_total = H_s - grad(gfu)
                    H_mag_sq = InnerProduct(H_total, H_total) + 1e-30
                    a_bf = BilinearForm(fes)
                    a_bf += chord_cf * grad(u) * grad(v) * dx(bonus_intorder=4)
                    a_bf += ((diff_cf - chord_cf) / H_mag_sq
                             * InnerProduct(H_total, grad(u))
                             * InnerProduct(H_total, grad(v))
                             ) * dx(bonus_intorder=4)
            else:
                # Picard: isotropic chord
                a_bf = BilinearForm(fes)
                if formulation == "omega":
                    a_bf += chord_cf * grad(u) * grad(v) * dx(bonus_intorder=4)
                else:
                    a_bf += chord_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
                    a_bf += 1e-6 * NU_0 * u * v * dx

            a_bf.Assemble()

            newton_step = newton and is_nonlinear and iteration > 0
            f_lf = LinearForm(fes)
            if newton_step:
                # Newton: the tangent solves for the correction to the
                # current iterate from the current nonlinear residual,
                # not for a new iterate from the chord source term.
                if formulation == "omega":
                    f_lf += chord_cf * (H_s - grad(gfu)) * grad(v) * dx(bonus_intorder=4)
                else:
                    f_lf += -chord_cf * (curl(gfu) + B_s) * curl(v) * dx(bonus_intorder=4)
                    f_lf += -1e-6 * NU_0 * gfu * v * dx
            elif formulation == "omega":
                f_lf += chord_cf * H_s * grad(v) * dx(bonus_intorder=4)
            else:
                f_lf += -chord_cf * B_s * curl(v) * dx(bonus_intorder=4)
            f_lf.Assemble()

            inv_a = _accel_inverse(a_bf, fes, linear_solver)
            if newton_step:
                # The correction is measured against the system load, the
                # chord source term a full (Picard) solve would use.
                f_load = LinearForm(fes)
                if formulation == "omega":
                    f_load += chord_cf * H_s * grad(v) * dx(bonus_intorder=4)
                else:
                    f_load += -chord_cf * B_s * curl(v) * dx(bonus_intorder=4)
                f_load.Assemble()
                load_free = np.array([bool(flag) for flag in fes.FreeDofs()])
                load_norm = float(np.linalg.norm(
                    np.asarray(f_load.vec.FV().NumPy())[load_free]))
                gf_delta = GridFunction(fes)
                with TaskManager():
                    linear_residual = apply_fe_inverse(
                        a_bf.mat, inv_a, f_lf.vec, gf_delta.vec, fes.FreeDofs(),
                        reference_norm=load_norm)
                gfu.vec.data += gf_delta.vec
            else:
                with TaskManager():
                    linear_residual = apply_fe_inverse(
                        a_bf.mat, inv_a, f_lf.vec, gfu.vec, fes.FreeDofs())

            t_solve_iter = time.perf_counter() - t0_iter
            n_iter_actual = iteration + 1

            if not is_nonlinear:
                _log(f"ITER:{iteration} linear t={t_solve_iter:.1f}s")
                converged = True
                break

            # Update material properties from BH curve
            max_change = 0.0
            if formulation == "omega":
                H_total = H_s - grad(gfu)
                gf_Hmag = GridFunction(fes_l2)
                gf_Hmag.Set(sqrt(InnerProduct(H_total, H_total)))
                for nr in yoke_elem_indices:
                    H_val = max(float(gf_Hmag.vec[nr]), 1e-10)
                    mu_new = float(bh_interp['mu_chord'](H_val))
                    mu_old = float(gf_prop.vec[nr])
                    mu_upd = (1 - relax) * mu_old + relax * mu_new
                    change = abs(mu_upd - mu_old) / max(abs(mu_old), 1e-30)
                    max_change = max(max_change, change)
                    gf_prop.vec[nr] = mu_upd
                    if newton:
                        gf_prop_diff.vec[nr] = float(
                            bh_interp['mu_diff'](H_val))
            else:
                B_total = curl(gfu) + B_s
                gf_Bmag = GridFunction(fes_l2)
                gf_Bmag.Set(sqrt(InnerProduct(B_total, B_total)))
                for nr in yoke_elem_indices:
                    B_val = max(float(gf_Bmag.vec[nr]), 1e-10)
                    nu_new = float(bh_interp['nu_chord'](B_val))
                    nu_old = float(gf_prop.vec[nr])
                    nu_upd = (1 - relax) * nu_old + relax * nu_new
                    change = abs(nu_upd - nu_old) / max(abs(nu_old), 1e-30)
                    max_change = max(max_change, change)
                    gf_prop.vec[nr] = nu_upd
                    if newton:
                        gf_prop_diff.vec[nr] = float(
                            bh_interp['nu_diff'](B_val))

            method = "Newton" if (newton and iteration > 0) else "Picard"
            history.append({
                'iteration': iteration,
                'max_change': float(max_change),
                't_solve': t_solve_iter,
                'method': method,
            })
            _log(f"ITER:{iteration} [{method}] max_change={max_change:.4e} "
                 f"t={t_solve_iter:.1f}s")

            if max_change < tol and iteration > 0:
                _log(f"CONVERGED:iter={iteration}")
                converged = True
                break

    # Residual of the discrete nonlinear equation at the returned iterate.
    # The loop stops on the relaxed material change; this re-evaluates the
    # material from the final field (no relaxation) and measures how well
    # that field satisfies the panel's own equations.
    nonlinear_relative_residual = None
    if is_nonlinear and not is_hysteresis:
        gf_check = GridFunction(fes_l2)
        gf_check.vec.data = gf_prop.vec
        if formulation == "omega":
            check_field = H_s - grad(gfu)
            law = bh_interp['mu_chord']
        else:
            check_field = curl(gfu) + B_s
            law = bh_interp['nu_chord']
        gf_mag = GridFunction(fes_l2)
        gf_mag.Set(sqrt(InnerProduct(check_field, check_field)))
        for nr in yoke_elem_indices:
            gf_check.vec[nr] = float(law(max(float(gf_mag.vec[nr]), 1e-10)))
        check_cf = gf_check * kelvin_weight
        a_check = BilinearForm(fes)
        f_check = LinearForm(fes)
        if formulation == "omega":
            a_check += check_cf * grad(u) * grad(v) * dx(bonus_intorder=4)
            f_check += check_cf * H_s * grad(v) * dx(bonus_intorder=4)
        else:
            a_check += check_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
            a_check += 1e-6 * NU_0 * u * v * dx
            f_check += -check_cf * B_s * curl(v) * dx(bonus_intorder=4)
        a_check.Assemble()
        f_check.Assemble()
        r_vec = f_check.vec.CreateVector()
        r_vec.data = a_check.mat * gfu.vec - f_check.vec
        free = np.array([bool(flag) for flag in fes.FreeDofs()])
        r_np = np.asarray(r_vec.FV().NumPy())[free]
        f_np = np.asarray(f_check.vec.FV().NumPy())[free]
        nonlinear_relative_residual = float(
            np.linalg.norm(r_np) / max(np.linalg.norm(f_np), 1e-300))
        _log(f"NONLINEAR:relative residual of the final iterate "
             f"{nonlinear_relative_residual:.3e}")

    # ============================================================
    # Step 7: Post-process
    # ============================================================
    if formulation == "omega":
        H_total = H_s - grad(gfu)
        B_field = gf_prop * H_total  # B = mu * H
    elif is_hysteresis:
        B_field = curl(gfu) + B_s  # use last B_s (final step)
    else:
        B_field = curl(gfu) + B_s

    # Magnetic energy (physical domain only, exclude kelvin)
    phys_dict = {}
    for m in materials:
        if "kelvin" in m.lower():
            phys_dict[m] = CF(0.0)
        else:
            phys_dict[m] = CF(1.0)
    phys_mask = mesh.MaterialCF(phys_dict, default=CF(1.0))

    if formulation == "omega":
        H_total = H_s - grad(gfu)
        W_mag = 0.5 * Integrate(
            phys_mask * gf_prop * InnerProduct(H_total, H_total),
            mesh, order=10).real
    elif is_hysteresis:
        # W = 0.5 * int nu_rev * |B|^2 dx (approximate, ignoring H_irr energy)
        W_mag = 0.5 * Integrate(
            phys_mask * gf_nu_const * InnerProduct(B_field, B_field),
            mesh, order=10).real
    else:
        W_mag = 0.5 * Integrate(
            phys_mask * gf_prop * InnerProduct(B_field, B_field),
            mesh, order=10).real

    # 0.5*int(B.H) equals the magnetic energy only for linear material.  It is
    # kept as W_half_BH; L = 2*W_half_BH/I^2 is the retained finite-domain B.H
    # inductance estimate (Kelvin region excluded), not a validated flux linkage.
    W_half_BH = W_mag
    W_coenergy = W_mag
    energy_convention = ("linear material: energy = coenergy = 0.5*int(B.H)")
    if is_hysteresis:
        energy_convention = ("hysteresis: 0.5*int(nu_rev*|B|^2) approximation; no unique "
                             "reversible energy or coenergy exists and none is reported")
        W_coenergy = None
    elif is_nonlinear:
        # Reversible nonlinear law: energy int(H dB) and coenergy int(B dH)
        # from the shared production B(H) law and its antiderivative,
        # evaluated on the solved field H.
        from radia.scalar_potential_solver import _build_bh_spline_law

        b_of, coenergy_of, _ = _build_bh_spline_law(bh_data)
        H_cf = H_s - grad(gfu) if formulation == "omega" else gf_prop * B_field
        H_mag = sqrt(InnerProduct(H_cf, H_cf))
        iron = mesh.MaterialCF({m: CF(1.0 if "yoke" in m.lower() else 0.0)
                                for m in materials}, default=CF(0.0))
        air_coenergy = 0.5 * gf_prop * InnerProduct(H_cf, H_cf) if formulation == "omega" \
            else 0.5 * InnerProduct(B_field, H_cf)
        iron_coenergy = coenergy_of(H_mag)
        iron_bh = b_of(H_mag) * H_mag
        W_coenergy = Integrate(phys_mask * (iron * iron_coenergy + (1 - iron) * air_coenergy),
                               mesh, order=10).real
        W_bh = Integrate(phys_mask * (iron * iron_bh + (1 - iron) * 2.0 * air_coenergy),
                         mesh, order=10).real
        W_mag = W_bh - W_coenergy
        energy_convention = ("reversible nonlinear: W_mag = int(H dB) and W_coenergy = "
                             "int(B dH) from the shared radia.bh_law antiderivative; "
                             "W_half_BH = 0.5*int(B.H) is not an energy")

    L = 2 * W_half_BH / current**2
    _log(f"ENERGY:W={W_mag:.6e} J, W'={W_coenergy}, "
         f"L(B.H estimate)={L*1e6:.4f} uH")

    # B at origin
    try:
        B_origin = [float(v) for v in B_field(mesh(0, 0, 0))]
        B_origin_mag = math.sqrt(sum(b**2 for b in B_origin))
    except Exception:
        B_origin = [0, 0, 0]
        B_origin_mag = 0

    t_total = time.perf_counter() - t_total_start
    _log(f"DONE:B0={B_origin_mag:.6f}T L={L*1e6:.4f}uH t={t_total:.1f}s")

    # ============================================================
    # Step 8: GMSH export
    # ============================================================
    # Phase B: save .vol + |B|.sol (NGSolve native), then vol2msh.
    gmsh_file = ""
    if msh_output:
        try:
            from gmsh_post_export import save_vol_sol_pair, vol2msh
            base_dir = os.path.dirname(os.path.abspath(msh_output))
            name_stem = os.path.splitext(os.path.basename(msh_output))[0]

            B_mag_cf = sqrt(InnerProduct(B_field, B_field))
            fes_h1 = H1(mesh, order=1)
            gf_B = GridFunction(fes_h1)
            gf_B.Set(B_mag_cf)

            vol_path = os.path.join(base_dir, f"{name_stem}_mesh.vol").replace("\\", "/")
            sol_path = os.path.join(base_dir, f"{name_stem}_Bmag.sol").replace("\\", "/")
            save_vol_sol_pair(vol_path, sol_path, mesh.ngmesh, gf_B)

            vol2msh(msh_output, vol_path, [
                {"sol": sol_path, "fes": "H1", "fes_order": 1, "fes_dim": 1,
                 "name": "|B|", "ncomp": 1},
            ])
            gmsh_file = msh_output
            _log(f"GMSH:{msh_output} (via vol2msh)")
        except Exception as e:
            _log(f"GMSH_ERROR:{e}")

    # ============================================================
    # Step 9: Result JSON
    # ============================================================
    result = {
        "B_origin": B_origin,
        "B_origin_mag": B_origin_mag,
        "W_mag": float(W_mag),
        "W_coenergy": None if W_coenergy is None else float(W_coenergy),
        "W_half_BH": float(W_half_BH),
        "energy_convention": energy_convention,
        "L": float(L),
        "L_convention": ("finite-domain estimate int(B.H)/I^2 over the physical "
                         "(non-Kelvin) region; not a validated flux linkage"),
        "ndof": ndof,
        "ne": ne,
        "formulation": formulation,
        "linear_solver": linear_solver,
        "linear_relative_residual": linear_residual,
        "iterations": n_iter_actual,
        "converged": converged,
        "nonlinear_relative_residual": nonlinear_relative_residual,
        "nonlinear": is_nonlinear,
        "t_mesh": round(t_mesh, 2),
        "t_total": round(t_total, 2),
        "current": current,
        "n_wire_segments": len(segments),
        "has_kelvin": has_kelvin,
        "hysteresis": is_hysteresis,
        "n_steps": len(step_results) if step_results else 1,
        "step_results": step_results,
        "msh_file": gmsh_file,
    }
    return result


def build_argparser():
    """argparse factory shared by main() and notebook DesignSpec callers."""
    parser = argparse.ArgumentParser(
        description="Accelerator magnet solver (Omega-reduced or A-formulation)")
    parser.add_argument("--coil-script", required=True,
                        help="Python script with build_coil() -> CoilBuilder")
    parser.add_argument("--vol", default="",
                        help="Netgen .vol file (sole interface between "
                             "Cubit and NGSolve per Cubit/NGSolve separation "
                             "policy -- this script no longer accepts .cub5)")
    parser.add_argument("--formulation", default="omega",
                        choices=["omega", "a"],
                        help="omega = scalar H1, a = vector HCurl")
    parser.add_argument("--fes-order", type=int, default=1,
                        help="FE polynomial order")
    add_material_args(parser, default_material="steel",
                      include_custom=False, include_hys=True)
    parser.add_argument("--n-steps", type=int, default=1,
                        help="Quasi-static steps (hysteresis: 0->I->0)")
    parser.add_argument("--max-iter", type=int, default=30,
                        help="Max Picard iterations")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="Convergence tolerance")
    parser.add_argument("--relax", type=float, default=0.3,
                        help="Fraction of the new material step per iteration, "
                             "in (0, 1]; 1 is the full step")
    parser.add_argument("--newton", action="store_true",
                        help="Use Newton iteration (default: Picard)")
    parser.add_argument("--solver", default="auto",
                        choices=FEM_SOLVERS,
                        help="auto (nonperiodic A p=1..3: BDDC+AMS; otherwise "
                             "SparseCholesky/BDDC by size), sparsecholesky, "
                             "bddc, iccg, bddc_ams (BDDC with AMS coarse solver). "
                             "Standalone ams is unsupported on this full HCurl space.")
    parser.add_argument("--msh-output", default="",
                        help="GMSH .msh output path")
    parser.add_argument("--output", default="",
                        help="JSON output file")
    return parser


def main():
    parser = build_argparser()

    def run(args):
        return solve_accel(
            coil_script=args.coil_script,
            vol_file=args.vol,
            formulation=args.formulation,
            fes_order=args.fes_order,
            mat=EMMaterial.from_args(args),
            n_steps=args.n_steps,
            max_iter=args.max_iter,
            tol=args.tol,
            relax=args.relax,
            newton=args.newton,
            solver=args.solver,
            msh_output=args.msh_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
