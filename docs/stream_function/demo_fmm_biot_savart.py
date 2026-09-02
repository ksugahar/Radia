"""FMM Biot-Savart benchmark: ngsolve.bem.BiotSavartCF vs LinearForm matrix path.

Compares two methods for evaluating the B field of a designed stream-function
coil at exterior observation points (r > coil bounding sphere):

  (A) LinearForm assembly (current design-matrix path)
      A_matrix @ psi_f  -- the same kernel as _assemble_biot_savart in
      calc_streamfunction.py.  Integrates K = n x grad(psi) analytically
      via a LinearForm for each eval point.  Valid at ANY distance.

  (B) ngsolve.bem.BiotSavartCF -- multipole-accelerated path.
      The surface current K = n x grad(psi) is represented as short elemental
      wire dipoles (one per surface triangle, centroid +/- K hat * eps),
      loaded via BiotSavartCF.AddCurrent.  Valid for eval_pts OUTSIDE the
      expansion sphere that encloses the coil.

      NOTE on BiotSavartCF(order, kappa, center, rad) in NGSolve 6.2.2604:
        This is a SPHERICAL-HARMONIC MULTIPOLE EXPANSION.  Sources (AddCurrent
        segments) must lie INSIDE the sphere (center, rad).  Evaluation at
        points OUTSIDE rad is the valid regime -- inside rad the series diverges
        and returns nan.  For a coil of radius 0.15 m and half-length 0.25 m,
        the bounding sphere radius is ~0.29 m; with the 1.5x safety margin the
        expansion sphere radius is ~0.44 m.  Exterior eval points at r >= 0.60 m
        are safely outside.

      Both paths also include a DSV INTERIOR check (Path A only) to confirm
      the designed coil actually produces the target Gz field inside the DSV.

Workflow:
  1. Mesh generation in subprocess (validation_test/panels/fixtures/
     make_streamfunction_vol.py).
  2. Design via in-process import of calc_streamfunction._build_problem.
  3. Evaluation:
     (A) LinearForm @ psi at DSV interior pts and exterior shell pts.
     (B) BiotSavartCF elemental K-dipoles at exterior shell pts only.
  4. Accuracy + timing table, JSON + timing bar chart.

Writes:
    demo_fmm_biot_savart.json  -- timing + accuracy results
    demo_fmm_biot_savart.png   -- timing bar chart
"""
import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FIXTURE = REPO / "tests" / "panels" / "fixtures" / "make_streamfunction_vol.py"
CALC = REPO / "src" / "radia" / "panels" / "calc_streamfunction.py"

sys.path.insert(0, str(REPO / "src" / "radia" / "panels"))
sys.path.insert(0, str(REPO / "src" / "radia"))

_MU0_4PI = 1.0e-7


# ---------------------------------------------------------------------------
# Subprocess helper for mesh generation
# ---------------------------------------------------------------------------

def _make_meshes(outdir):
    """Generate coil + eval .vol meshes via the canonical fixture (subprocess)."""
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([sys.executable, str(FIXTURE), str(outdir)],
                       capture_output=True, text=True, env=env, timeout=600)
    if r.returncode != 0:
        raise RuntimeError(f"mesh gen failed (rc={r.returncode}):\n{r.stderr[-600:]}")
    coil = os.path.join(outdir, "coil_cyl_surf.vol")
    eval_v = os.path.join(outdir, "eval_dsv.vol")
    if not (os.path.exists(coil) and os.path.exists(eval_v)):
        raise RuntimeError(f"mesh gen did not produce .vol files")
    return coil, eval_v


# ---------------------------------------------------------------------------
# Path A: LinearForm assembly @ psi_full
# ---------------------------------------------------------------------------

def _eval_matrix_path(coil_vol, psi_full, order, eval_pts):
    """Evaluate Bz at eval_pts via the LinearForm matrix path.

    Assembles one LinearForm row per eval point (the _assemble_biot_savart
    kernel) and computes A_row @ psi_full.  Returns (Bz, elapsed_s).

    This is the EXACT surface-current integral valid at ANY point.
    """
    from ngsolve import (Mesh, H1, specialcf, grad, LinearForm, ds,
                         sqrt as ng_sqrt, Cross, x as ng_x, y as ng_y,
                         z as ng_z, TaskManager)

    with TaskManager():
        coil = Mesh(coil_vol)
        fes = H1(coil, order=order, definedon=coil.Boundaries(".*"))
        n = specialcf.normal(3)
        ndof = fes.ndof
        v = fes.TestFunction()

        t0 = time.perf_counter()
        A_row = np.zeros((len(eval_pts), ndof))
        for i, p in enumerate(eval_pts):
            dxt = p[0] - ng_x
            dyt = p[1] - ng_y
            dzt = p[2] - ng_z
            r2 = dxt * dxt + dyt * dyt + dzt * dzt
            r3 = r2 * ng_sqrt(r2)
            Kv = Cross(n, grad(v).Trace())
            # Bz = (K x r)_z = Kx * dy - Ky * dx (with r pointing FROM source TO obs)
            # r = obs - src, Biot-Savart: dB = mu0/(4pi) K x r / r^3 dS
            cross_z = Kv[0] * dyt - Kv[1] * dxt
            Lf = LinearForm(fes)
            Lf += (_MU0_4PI / r3) * cross_z * ds
            Lf.Assemble()
            A_row[i, :] = Lf.vec.FV().NumPy()
        elapsed = time.perf_counter() - t0

    Bz = A_row @ psi_full
    return Bz, elapsed


# ---------------------------------------------------------------------------
# Path B: BiotSavartCF elemental K-dipoles
# ---------------------------------------------------------------------------

def _build_K_dipoles(coil_vol, psi_full, order):
    """Represent surface current K = n x grad(psi) as elemental wire dipoles.

    Each surface triangle -> one short filament at its centroid with length
    eps ~ sqrt(area) and current I = |K_avg| * area / (2 eps).  This is the
    standard wire-approximation of a surface current element used in the
    IH BEM Biot-Savart post-processor.  Returns list of (p1, p2, I) tuples
    where p1, p2 are numpy arrays.
    """
    from ngsolve import (Mesh, H1, specialcf, grad, GridFunction, Cross,
                         Integrate, CF, BND, HDivSurface, TaskManager)

    with TaskManager():
        coil = Mesh(coil_vol)
        fes = H1(coil, order=order, definedon=coil.Boundaries(".*"))
        n = specialcf.normal(3)

        # Recover psi GridFunction
        gfpsi = GridFunction(fes)
        gfpsi.vec.FV().NumPy()[:] = psi_full

        # Project K = n x grad(psi) onto HDivSurface (RT0 surface current)
        fesJ = HDivSurface(coil, order=0)
        gfJ = GridFunction(fesJ)
        gfJ.Set(Cross(n, grad(gfpsi).Trace()),
                definedon=coil.Boundaries(".*"), dual=True)

        # Per-element K averages and areas
        elem_A = Integrate(CF(1), coil, VOL_or_BND=BND, element_wise=True)
        elem_Jx = Integrate(gfJ[0], coil, VOL_or_BND=BND, element_wise=True)
        elem_Jy = Integrate(gfJ[1], coil, VOL_or_BND=BND, element_wise=True)
        elem_Jz = Integrate(gfJ[2], coil, VOL_or_BND=BND, element_wise=True)

    dipoles = []
    for el in coil.Elements(BND):
        area = abs(float(elem_A[el.nr]))
        if area < 1e-30:
            continue
        Jvec = np.array([float(elem_Jx[el.nr]),
                         float(elem_Jy[el.nr]),
                         float(elem_Jz[el.nr])]) / area
        J_mag = np.linalg.norm(Jvec)
        if J_mag < 1e-30:
            continue

        verts = [coil.vertices[v.nr].point for v in el.vertices]
        centroid = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)

        eps = math.sqrt(area) * 0.5
        J_hat = Jvec / J_mag
        p1 = centroid - eps * J_hat
        p2 = centroid + eps * J_hat
        I_fil = J_mag * area / (2.0 * eps)

        dipoles.append((p1, p2, I_fil))

    return dipoles


def _eval_bem_path(dipoles, eval_pts, order=8):
    """Evaluate Bz at eval_pts using ngsolve.bem.BiotSavartCF.

    Loads elemental K-dipoles as AddCurrent segments and evaluates via
    the spherical-harmonic multipole expansion.

    Returns (Bz, elapsed_s, expansion_rad).

    Valid ONLY for eval_pts with distance from expansion center > expansion_rad.
    Points inside return nan.
    """
    from ngsolve import Mesh, H1, GridFunction, TaskManager
    from ngsolve.bem import BiotSavartCF
    from ngsolve.bla import Vec3D
    from netgen.occ import Pnt, OCCGeometry, Sphere

    MU0 = 4.0e-7 * math.pi

    # Expansion sphere: enclose all dipole source points.
    all_pts = np.array([pt for (p1, p2, _) in dipoles for pt in (p1, p2)])
    center_xyz = all_pts.mean(axis=0)
    dists = np.linalg.norm(all_pts - center_xyz[np.newaxis, :], axis=1)
    # Standard safety margin: 1.5x the source bounding radius.
    rad = float(np.max(dists)) * 1.5 + 1e-6

    eval_arr = np.asarray(eval_pts, float)
    eval_dists = np.linalg.norm(eval_arr - center_xyz[np.newaxis, :], axis=1)
    n_inside = int(np.sum(eval_dists < rad))
    if n_inside > 0:
        print(f"    WARNING: {n_inside}/{len(eval_arr)} eval pts INSIDE "
              f"expansion sphere (rad={rad:.3f} m from center) -- nan expected.")

    t0 = time.perf_counter()
    bs = BiotSavartCF(order=order, kappa=1e-10,
                      center=Vec3D(*center_xyz.tolist()), rad=rad)
    for p1, p2, I in dipoles:
        bs.AddCurrent(Vec3D(*p1.tolist()), Vec3D(*p2.tolist()), complex(I))

    # BiotSavartCF returns H [A/m]; B = mu0 * H.
    # Evaluate via GridFunction.Set on a thin spherical shell mesh that
    # encloses the eval_pts.  Using a shell (not a filled box) reduces
    # the mesh size dramatically and avoids the slow 3D box meshing.
    eval_R = float(np.max(np.linalg.norm(eval_arr, axis=1))) * 1.05
    shell_outer = Sphere(Pnt(0, 0, 0), eval_R + 0.05)
    shell_inner = Sphere(Pnt(0, 0, 0), eval_R - 0.05)
    shell_geom = shell_outer - shell_inner

    with TaskManager():
        ng_shell = OCCGeometry(shell_geom).GenerateMesh(maxh=0.15)

    from ngsolve import Mesh as NGMesh, H1 as NGH1, GridFunction as NGGF
    shell_mesh = NGMesh(ng_shell)
    fes_sh = NGH1(shell_mesh, order=1, complex=True)
    gf_Bz = NGGF(fes_sh)

    with TaskManager():
        gf_Bz.Set((MU0 * bs)[2])

    Bz = np.empty(len(eval_arr))
    for i, p in enumerate(eval_arr):
        mp = shell_mesh(p[0], p[1], p[2])
        val = gf_Bz(mp)
        Bz[i] = float(val.real) if hasattr(val, "real") else float(val)

    elapsed = time.perf_counter() - t0
    return Bz, elapsed, rad


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot(result, png):
    """Timing bar chart (no in-figure title per lab convention)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
    except ImportError:
        print("(matplotlib not available; skipping plot)")
        return

    labels = ["LinearForm\n(A, DSV)", "LinearForm\n(A, ext)", "BiotSavartCF\n(B, ext)"]
    times = [result["t_linear_form_dsv_s"],
             result["t_linear_form_ext_s"],
             result["t_bem_exterior_s"]]

    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    colors = ["steelblue", "cornflowerblue", "tomato"]
    bars = ax.bar(labels, times, color=colors, width=0.5)
    ax.set_ylabel("Wall time (s)")
    ax.set_yscale("log")
    for bar, t in zip(bars, times):
        if t > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, t * 1.15,
                    f"{t:.2f} s", ha="center", va="bottom", fontsize=8)
    ax.text(0.97, 0.95,
            f"N_eval_dsv={result['n_eval_dsv']}\n"
            f"N_eval_ext={result['n_eval_ext']}\n"
            f"ndof={result['coil']['ndof']}\n"
            f"n_dipoles={result['n_dipoles']}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7))
    ax.grid(axis="y", alpha=0.4)
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)
    print(f"Saved plot: {Path(png).name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Benchmark: LinearForm vs BiotSavartCF for SF coil Bz.")
    ap.add_argument("--order", type=int, default=2,
                    help="psi FE order (default 2)")
    ap.add_argument("--n-eval-dsv", type=int, default=60,
                    help="evaluation points inside DSV for Path A (default 60)")
    ap.add_argument("--n-eval-ext", type=int, default=60,
                    help="exterior eval points for both paths (default 60)")
    ap.add_argument("--ext-radius", type=float, default=0.60,
                    help="exterior shell radius [m]. Must be outside the "
                         "BiotSavartCF expansion sphere (~0.44 m for this coil). "
                         "(default 0.60)")
    ap.add_argument("--bem-order", type=int, default=8,
                    help="BiotSavartCF multipole order (default 8)")
    ap.add_argument("--outdir", type=str, default=None,
                    help="output directory for temp mesh files (default: temp, cleaned)")
    args = ap.parse_args()

    # ---- mesh generation (subprocess) -------------------------------------
    own_tmpdir = args.outdir is None
    if own_tmpdir:
        tmpdir_obj = tempfile.TemporaryDirectory(prefix="demo_fmm_")
        tmpdir = tmpdir_obj.name
    else:
        tmpdir = args.outdir
        os.makedirs(tmpdir, exist_ok=True)

    print("Generating coil + eval meshes via fixture (subprocess)...")
    t_mesh0 = time.perf_counter()
    coil_vol, eval_vol = _make_meshes(tmpdir)
    t_mesh = time.perf_counter() - t_mesh0
    print(f"  coil_vol: {Path(coil_vol).name}  eval_vol: {Path(eval_vol).name}")
    print(f"  t_mesh = {t_mesh:.1f} s")

    # ---- in-process design via calc_streamfunction import -----------------
    print(f"\nBuilding design problem (Gz coil, order={args.order})...")
    from ngsolve import Mesh, H1, specialcf, TaskManager, GridFunction
    import calc_streamfunction as CS

    class _A:
        pass

    sf = _A()
    sf.coil_vol = coil_vol
    sf.eval_vol = eval_vol
    sf.target_cf = "z"
    sf.target_harmonic = ""
    sf._harmonic_terms = None
    sf.order = args.order
    sf.regularize = "h1"
    sf.alpha = 0.0
    sf.confine = "off"
    sf.eval_max = 400

    t_build0 = time.perf_counter()
    with TaskManager():
        P = CS._build_problem(sf)
        psi_f, _, homo = CS._solve_and_metrics(P, 0.0)
    psi_full = P["R"] @ psi_f
    t_build = time.perf_counter() - t_build0
    ndof = int(P["fes"].ndof)
    print(f"  ndof={ndof}, homogeneity_rms={homo:.3e}, t_build={t_build:.1f} s")

    # ---- evaluation points ------------------------------------------------
    dsv_mesh = Mesh(eval_vol)
    dsv_pts_all = np.array([list(v.point) for v in dsv_mesh.vertices])
    rng = np.random.default_rng(42)
    idx = rng.choice(len(dsv_pts_all),
                     size=min(args.n_eval_dsv, len(dsv_pts_all)), replace=False)
    dsv_pts = dsv_pts_all[idx]

    R_ext = args.ext_radius
    n_ext = args.n_eval_ext
    golden = (1.0 + math.sqrt(5.0)) / 2.0
    i_arr = np.arange(n_ext)
    theta = np.arccos(1.0 - 2.0 * (i_arr + 0.5) / n_ext)
    phi = 2.0 * math.pi * i_arr / golden
    ext_pts = R_ext * np.column_stack([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta)])

    print(f"\nEvaluation: {len(dsv_pts)} DSV interior pts + "
          f"{len(ext_pts)} exterior shell pts (R={R_ext:.2f} m)")

    # ---- Path A: LinearForm @ psi at DSV interior -------------------------
    print("\n(A) LinearForm at DSV interior pts...")
    Bz_A_dsv, t_A_dsv = _eval_matrix_path(coil_vol, psi_full, args.order, dsv_pts)
    bz_rms_dsv = float(np.sqrt(np.mean(Bz_A_dsv ** 2)))
    print(f"    Bz_rms(DSV)={bz_rms_dsv:.3e} T   t={t_A_dsv:.2f} s")

    # ---- Path A: LinearForm @ psi at exterior shell -----------------------
    print("\n(A) LinearForm at exterior shell pts...")
    Bz_A_ext, t_A_ext = _eval_matrix_path(coil_vol, psi_full, args.order, ext_pts)
    bz_rms_ext_A = float(np.sqrt(np.mean(Bz_A_ext ** 2)))
    print(f"    Bz_rms(ext)={bz_rms_ext_A:.3e} T   t={t_A_ext:.2f} s")

    # ---- Build elemental K-dipoles for Path B -----------------------------
    print("\nBuilding elemental K-dipoles from psi (for BiotSavartCF)...")
    t_dip0 = time.perf_counter()
    dipoles = _build_K_dipoles(coil_vol, psi_full, args.order)
    t_dip = time.perf_counter() - t_dip0
    print(f"  n_dipoles={len(dipoles)}  t_build_dipoles={t_dip:.2f} s")

    # ---- Path B: BiotSavartCF at exterior shell ---------------------------
    print(f"\n(B) BiotSavartCF (ngsolve.bem, order={args.bem_order}) at exterior shell...")
    if len(dipoles) < 2:
        print("    No dipoles -- skipping BiotSavartCF")
        Bz_B_ext = np.full(len(ext_pts), float("nan"))
        t_B_ext = 0.0
        rad_exp = float("nan")
    else:
        Bz_B_ext, t_B_ext, rad_exp = _eval_bem_path(
            dipoles, ext_pts, order=args.bem_order)
        print(f"    expansion_rad={rad_exp:.3f} m  (eval R={R_ext:.2f} > {rad_exp:.3f}? "
              f"{'YES -- exterior OK' if R_ext > rad_exp else 'NO -- inside sphere!'})")

    valid = ~np.isnan(Bz_B_ext)
    if valid.any():
        bz_rms_ext_B = float(np.sqrt(np.mean(Bz_B_ext[valid] ** 2)))
        print(f"    Bz_rms(ext)={bz_rms_ext_B:.3e} T   t={t_B_ext:.2f} s")
    else:
        bz_rms_ext_B = float("nan")
        print(f"    All nan  t={t_B_ext:.2f} s")

    # ---- relative error A vs B at exterior --------------------------------
    if valid.any() and np.linalg.norm(Bz_A_ext[valid]) > 1e-30:
        rel_err_ext = float(
            np.linalg.norm(Bz_A_ext[valid] - Bz_B_ext[valid])
            / np.linalg.norm(Bz_A_ext[valid]))
    else:
        rel_err_ext = float("nan")

    # ---- summary table ----------------------------------------------------
    print()
    print("=== Biot-Savart methods (Gz coil, order=%d, ndof=%d) ===" % (args.order, ndof))
    print("    N_eval_dsv=%d (interior),  N_eval_ext=%d (exterior R=%.2f m)"
          % (len(dsv_pts), len(ext_pts), R_ext))
    print("  (A) LinearForm (DSV interior)     : t=%6.2f s   Bz_rms=%9.2e T"
          % (t_A_dsv, bz_rms_dsv))
    print("  (A) LinearForm (exterior shell)   : t=%6.2f s   Bz_rms=%9.2e T"
          % (t_A_ext, bz_rms_ext_A))
    if not math.isnan(bz_rms_ext_B):
        print("  (B) BiotSavartCF (exterior shell) : t=%6.2f s   Bz_rms=%9.2e T"
              "  rel_err=%.2e" % (t_B_ext, bz_rms_ext_B, rel_err_ext))
    else:
        print("  (B) BiotSavartCF (exterior shell) : t=%6.2f s   Bz_rms=nan"
              % t_B_ext)

    # ---- save results -----------------------------------------------------
    out_json = HERE / "demo_fmm_biot_savart.json"
    out_png = HERE / "demo_fmm_biot_savart.png"

    result = {
        "coil": {"radius_m": 0.15, "length_m": 0.50,
                 "ndof": ndof, "order": args.order},
        "design": {"target_cf": "z", "homogeneity_rms": homo,
                   "t_build_s": round(t_build, 2)},
        "n_eval_dsv": int(len(dsv_pts)),
        "n_eval_ext": int(len(ext_pts)),
        "ext_radius_m": R_ext,
        "n_dipoles": len(dipoles),
        "expansion_rad_m": rad_exp if not math.isnan(rad_exp) else None,
        "t_mesh_s": round(t_mesh, 2),
        "t_build_dipoles_s": round(t_dip, 3),
        "t_linear_form_dsv_s": round(t_A_dsv, 3),
        "t_linear_form_ext_s": round(t_A_ext, 3),
        "t_bem_exterior_s": round(t_B_ext, 3),
        "Bz_rms_dsv_A": bz_rms_dsv,
        "Bz_rms_ext_A": bz_rms_ext_A,
        "Bz_rms_ext_B": bz_rms_ext_B,
        "rel_err_ext_A_vs_B": rel_err_ext,
        "note": (
            "Path A uses continuous surface K = n x grad(psi) via LinearForm. "
            "Path B uses per-triangle elemental K-dipoles (same source, different "
            "discretisation) via BiotSavartCF multipole expansion. "
            "Relative error measures the dipole discretisation error vs the "
            "continuous-K LinearForm reference."
        ),
    }

    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved {out_json.name} + {out_png.name}")

    _plot(result, str(out_png))

    if own_tmpdir:
        try:
            tmpdir_obj.cleanup()
        except Exception:
            pass


if __name__ == "__main__":
    main()
