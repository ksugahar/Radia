#!/usr/bin/env python3
"""Active-shielding demo: primary + shield gradient coil (Turner 1986).

Designs a Gz gradient coil two ways on the SAME inner former and compares the
STRAY field outside the assembly:

  (1) UNSHIELDED -- primary cylinder only (the inner coil).
  (2) SHIELDED   -- primary + an outer SHIELD cylinder designed JOINTLY so
      that primary + shield together still produce the Gz target inside the
      DSV while nulling the field at external sample points (the
      calc_streamfunction --shield-vol path, _build_shielded_problem).

Both designs are evaluated at INDEPENDENT external points (element centroids
that are NOT in the constraint set) so the reported stray is the HONEST
generalising field, not the circular constraint-point residual.

Honest result (primary r=0.15 m, shield r=0.20 m, DSV r=0.05 m, external null
region = a LARGE annular shell r=0.235-0.27 m spanning the full coil length,
~800 constraint points):
  unshielded stray ~ 1.65   (relative to the DSV target field)
  shielded   stray ~ 0.019
  => a genuine ~86x (~39 dB) external-field reduction, DSV homogeneity
     preserved; the stray-vs-radius profile shows 100-1700x suppression in
     the well-covered far region (r > 0.25 m).

CRITICAL LESSON (measured): the external null region must COVER the exterior,
not just a thin mid-plane slice.  With a too-small region (e.g. r=0.21-0.23,
L=0.18) the point-sampled nulling overfits that slice and gives only ~4x
locally while making the field WORSE at larger z.  With the full-length shell
here the shielding is real and strong everywhere outside it.  The remaining
limitation is the un-sampled GAP between the shield (r=0.20) and the null
region (r=0.235): the r=0.21 point sees only ~2-3x (it is not constrained).
An analytic external-multipole-moment constraint would remove the sampling
dependence entirely; that is the documented next step.

The design uses the SAME production code paths as the CLI/panel
(calc_streamfunction._build_shielded_problem / _build_problem); mesh-gen runs
in a subprocess so the orchestrator setup is clean.

Usage::

    python demo_active_shield.py [--order 2] [--shield-weight 1.0]

Outputs (committed next to this script, per the Data Persistence Policy):
    demo_active_shield.json   -- stray + homogeneity for both designs
    demo_active_shield.png    -- stray bar chart + stray-vs-radius profile
"""
import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FIXTURE = REPO / "tests" / "panels" / "fixtures" / "make_shielded_coil_vol.py"

sys.path.insert(0, str(REPO / "src" / "radia" / "panels"))
sys.path.insert(0, str(REPO / "src" / "radia"))


def _make_meshes(outdir):
    """Generate primary + shield + DSV + external .vol via the fixture."""
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([sys.executable, str(FIXTURE), str(outdir)],
                       capture_output=True, text=True, env=env, timeout=600)
    if r.returncode != 0:
        raise RuntimeError(f"mesh gen failed (rc={r.returncode}):\n{r.stderr[-800:]}")
    paths = {k: os.path.join(outdir, v) for k, v in {
        "primary": "primary_cyl.vol", "shield": "shield_cyl.vol",
        "dsv": "dsv.vol", "ext": "external.vol"}.items()}
    if not all(os.path.exists(p) for p in paths.values()):
        raise RuntimeError("fixture did not produce all 4 .vol files")
    return paths


class _Args:
    """Minimal argparse-Namespace stand-in for _build_problem helpers."""
    def __init__(self, **kw):
        self.order = 2
        self.regularize = "h1"
        self.alpha = 0.0
        self.confine = "abe"
        self.eval_max = 800
        self.target_cf = ""
        self.target_harmonic = "Z"
        self._harmonic_terms = None
        self.shield_weight = 1.0
        for k, v in kw.items():
            setattr(self, k, v)


def _field_at(CS, fes, R, n, pts, psi_free):
    """Bz (z-component) at pts from one coil's free-DOF psi."""
    A = CS._assemble_biot_savart(fes, n, pts, [2])
    AR = np.asarray((R.T @ A.T).T)
    return AR @ psi_free


def _shell_points(radius, n_phi=16, n_z=7, z_max=0.22):
    """Points on a cylindrical shell at the given radius, spanning the coil's
    z-extent.  A Gz coil has Bz=0 at the midplane (z=0) by symmetry, so a
    midplane-only ring is a misleading probe -- the shell spans z to capture
    the real stray."""
    t = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    zs = np.linspace(-z_max, z_max, n_z)
    pts = []
    for z in zs:
        pts.append(np.column_stack([radius * np.cos(t), radius * np.sin(t),
                                    np.full(n_phi, z)]))
    return np.vstack(pts)


def main():
    ap = argparse.ArgumentParser(
        description="Active-shielding demo (primary + shield Gz coil).")
    ap.add_argument("--order", type=int, default=2, help="psi FE order")
    ap.add_argument("--shield-weight", type=float, default=1.0,
                    help="weight on the stray-null constraint rows")
    ap.add_argument("--out-dir", default=str(HERE))
    ap.add_argument("--work-dir", default=None,
                    help="mesh scratch (default: fresh temp dir)")
    args = ap.parse_args()

    work = args.work_dir or tempfile.mkdtemp(prefix="sf_shield_")
    os.makedirs(work, exist_ok=True)
    print("Generating primary + shield + DSV + external meshes (subprocess)...")
    vols = _make_meshes(work)

    from ngsolve import Mesh, TaskManager
    import calc_streamfunction as CS

    with TaskManager():
        # --- SHIELDED design (primary + shield, joint solve) ---
        a_sh = _Args(coil_vol=vols["primary"], eval_vol=vols["dsv"],
                     shield_vol=vols["shield"], shield_eval_vol=vols["ext"],
                     order=args.order, shield_weight=args.shield_weight)
        P_sh = CS._build_shielded_problem(a_sh)
        psi_sh, _, homo_sh = CS._solve_and_metrics(P_sh, 0.0)
        n_fp = P_sh["n_primary"]
        psi_sh_p, psi_sh_s = psi_sh[:n_fp], psi_sh[n_fp:]

        # --- UNSHIELDED design (primary only, same target) ---
        a_un = _Args(coil_vol=vols["primary"], eval_vol=vols["dsv"],
                     order=args.order)
        P_un = CS._build_problem(a_un)
        psi_un, _, homo_un = CS._solve_and_metrics(P_un, 0.0)

        # --- independent external measure points (centroids) ---
        ext_mesh = Mesh(vols["ext"])
        e_all = CS._element_centroids(ext_mesh)
        e_mpts = e_all[np.linspace(0, len(e_all) - 1,
                                   min(400, len(e_all))).astype(int)]
        B_dsv_norm = float(np.linalg.norm(P_sh["Bc_dsv"])) + 1e-30

        # unshielded stray = primary field at e_mpts
        Bz_un_ext = _field_at(CS, P_un["fes"], P_un["R"], P_un["n"],
                              e_mpts, psi_un)
        stray_un = float(np.linalg.norm(Bz_un_ext) / B_dsv_norm)

        # shielded stray = primary + shield field at e_mpts
        Bz_sh_ext = (_field_at(CS, P_sh["fes"], P_sh["R_p"], P_sh["n"],
                               e_mpts, psi_sh_p)
                     + _field_at(CS, P_sh["shield_fes"], P_sh["R_s"],
                                 P_sh["shield_n"], e_mpts, psi_sh_s))
        stray_sh = float(np.linalg.norm(Bz_sh_ext) / B_dsv_norm)

        # --- stray-vs-radius profile (cylindrical shells spanning z) ---
        radii = [0.21, 0.25, 0.30, 0.40, 0.55, 0.75]
        prof_un, prof_sh = [], []
        for rr in radii:
            pts = _shell_points(rr)
            bz_u = _field_at(CS, P_un["fes"], P_un["R"], P_un["n"], pts, psi_un)
            bz_s = (_field_at(CS, P_sh["fes"], P_sh["R_p"], P_sh["n"],
                              pts, psi_sh_p)
                    + _field_at(CS, P_sh["shield_fes"], P_sh["R_s"],
                                P_sh["shield_n"], pts, psi_sh_s))
            prof_un.append(float(np.sqrt(np.mean(bz_u ** 2)) / B_dsv_norm))
            prof_sh.append(float(np.sqrt(np.mean(bz_s ** 2)) / B_dsv_norm))

    factor = stray_un / max(stray_sh, 1e-30)
    print("\n=== Active shielding (Gz coil, primary r=0.15, shield r=0.20) ===")
    print(f"  DSV homogeneity:  unshielded {homo_un:.2e}   shielded {homo_sh:.2e}")
    print(f"  stray @ external: unshielded {stray_un:.3e}  shielded {stray_sh:.3e}")
    print(f"  shielding factor  = {factor:.2f}x  ({20*math.log10(max(factor,1e-9)):.1f} dB)")
    print("  stray-vs-radius (relative to DSV target):")
    for rr, u, s in zip(radii, prof_un, prof_sh):
        print(f"    r={rr*1e3:4.0f} mm:  unshielded {u:.3e}   shielded {s:.3e}"
              f"   ({u/max(s,1e-30):.1f}x)")

    result = {
        "order": args.order, "shield_weight": args.shield_weight,
        "geometry": {"primary_r_m": 0.15, "shield_r_m": 0.20,
                     "dsv_r_m": 0.05},
        "homogeneity_unshielded": homo_un, "homogeneity_shielded": homo_sh,
        "stray_unshielded": stray_un, "stray_shielded": stray_sh,
        "shielding_factor": factor,
        "n_external_measure": int(len(e_mpts)),
        "profile_radii_m": radii,
        "profile_stray_unshielded": prof_un,
        "profile_stray_shielded": prof_sh,
        "note": (
            "Genuine ~86x (39 dB) stray reduction with a LARGE external null "
            "region (full-length shell r=0.235-0.27).  Reported stray is at "
            "INDEPENDENT external points (centroids), the honest generalising "
            "field -- NOT the circular constraint-point residual.  The "
            "external region MUST cover the exterior: a thin mid-plane slice "
            "overfits locally (~4x) and worsens the field at larger z.  The "
            "un-sampled gap between shield (r=0.20) and null region (r=0.235) "
            "is only ~2-3x shielded (not constrained)."
        ),
    }
    jpath = os.path.join(args.out_dir, "demo_active_shield.json")
    with open(jpath, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved data: {os.path.basename(jpath)}")

    _plot(result, os.path.join(args.out_dir, "demo_active_shield.png"))


def _plot(result, png):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
    except ImportError:
        print("(matplotlib not installed; skipped plot)")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))

    # (a) stray bar chart: unshielded vs shielded
    vals = [result["stray_unshielded"], result["stray_shielded"]]
    ax1.bar(["unshielded\n(primary only)", "shielded\n(primary + shield)"],
            vals, color=["#d95f0e", "#2c7fb8"], width=0.55)
    ax1.set_ylabel("stray field @ external pts  (rel. to DSV target)")
    for i, v in enumerate(vals):
        ax1.text(i, v * 1.02, f"{v:.3f}", ha="center", va="bottom", fontsize=10)
    ax1.text(0.5, 0.90, f"{result['shielding_factor']:.1f}x reduction",
             transform=ax1.transAxes, ha="center", fontsize=11,
             bbox=dict(boxstyle="round,pad=0.3", fc="#ffffcc", ec="gray"))
    ax1.set_ylim(0, max(vals) * 1.25)
    ax1.grid(axis="y", alpha=0.3)
    ax1.text(0.02, 0.96, "(a)", transform=ax1.transAxes, fontsize=11, va="top")

    # (b) stray-vs-radius profile
    radii_mm = [r * 1e3 for r in result["profile_radii_m"]]
    ax2.semilogy(radii_mm, result["profile_stray_unshielded"], "o-",
                 color="#d95f0e", label="unshielded")
    ax2.semilogy(radii_mm, result["profile_stray_shielded"], "s-",
                 color="#2c7fb8", label="shielded")
    ax2.axvline(150, color="gray", ls=":", lw=0.8)
    ax2.axvline(200, color="gray", ls=":", lw=0.8)
    ax2.set_xlabel("radius (mm)   [primary=150, shield=200]")
    ax2.set_ylabel("shell-RMS stray |Bz|  (rel. to DSV target)")
    ax2.legend()
    ax2.grid(alpha=0.3, which="both")
    ax2.text(0.02, 0.06, "(b)", transform=ax2.transAxes, fontsize=11,
             va="bottom")

    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)
    print(f"Saved plot: {os.path.basename(png)}")


if __name__ == "__main__":
    main()
