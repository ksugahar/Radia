"""Can a surface impedance be repaired where it fails, instead of abandoned?

The 2-D surface-impedance solve of the beak section is converged to 1024
perimeter panels and is still wrong, because the error is in the model rather
than the mesh.  A Leontovich impedance is the leading term of an expansion in
the skin depth over the local radius of curvature; at 150 kHz this section's
tip radius is 0.25 mm against a 0.171 mm skin, so that ratio is 0.68, and its
beak root and far end carry square corners, where the radius is zero and no
term of the expansion applies at all.

This measures the damage and then repairs it the way Proekt, Yuferev,
Tsukerman and Ida (2002) propose: resolve the conductor's interior in a
neighbourhood of each offending feature, leave the surface impedance in charge
only of the flat band between them, and match the two across a chord of the
thick body, where the outer solution is accurate in the sense its own
expansion claims.  Dauge, Dular, Krahenbuhl, Peron, Perrussel and Poignard
(2014) place the corner layer at the scale of the skin depth, which is why a
patch of a few skin depths contains it and a match several hundred skin depths
back is safely outside it.

Everything is scored against a control that resolves the whole section in the
plane, driven by the exact exterior field of the same total current.  That
control is anchored on the Bessel round wire, which it reproduces to five
significant figures, so it is not itself a thing to be trusted on assertion.
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

MU0 = 4e-7 * np.pi
BEAK_X_M = 1.6e-3
TIP_X_M = 3.5e-3


def section_metres(tmp_step):
    """The tracked beak profile as an OCC face, in metres."""
    from build123d import Plane, export_step
    from make_beak_fin_step import _beak_face
    from netgen.occ import OCCGeometry, Pnt

    export_step(_beak_face(Plane.XY), str(tmp_step))
    shape = OCCGeometry(str(tmp_step), dim=2).shape.Scale(Pnt(0, 0, 0), 1e-3)
    tmp_step.unlink(missing_ok=True)
    return shape


def round_wire_anchor(sigma, order, maxh_over_delta):
    """Reproduce the exact Bessel round-wire resistance with the same solver.

    The exterior of a round wire is that of a line current whatever the skin
    depth, so the open-boundary data is known in closed form; the solve is
    affine in the uniform axial field, so two runs fix the one that carries
    exactly 1 A.  Nothing is fitted between the solver and the exact result.
    """
    import ngsolve as ng
    from scipy.special import jv

    from netgen.occ import WorkPlane
    from radia.sibc_corner_patch import solve_cut_patch

    a = 1.0e-3
    rows = []
    for frequency in (50_000.0, 200_000.0, 800_000.0):
        omega = 2.0 * np.pi * frequency
        delta = float(np.sqrt(2.0 / (omega * MU0 * sigma)))
        radius_cf = ng.sqrt(ng.x * ng.x + ng.y * ng.y)
        exterior = (MU0 / (2.0 * np.pi)) * ng.log(a / radius_cf)
        theta = np.linspace(0.0, 2.0 * np.pi, 512, endpoint=False)
        panel_xy = np.stack([a * np.cos(theta), a * np.sin(theta)], axis=1)
        panel_ds = np.full(len(theta), 2.0 * np.pi * a / len(theta))
        panel_current = np.full(len(theta), 1.0 / len(theta), dtype=complex)
        face = WorkPlane().Circle(0.0, 0.0, a).Face()
        common = dict(x_cut=None, panel_xy=panel_xy, panel_ds=panel_ds,
                      panel_current=panel_current,
                      surface_impedance=(1 + 1j) / (sigma * delta),
                      omega=omega, sigma=sigma,
                      maxh_conductor=maxh_over_delta * delta,
                      air_pad=20.0 * a, maxh_air=0.5 * a, order=order,
                      exterior_potential=exterior)
        currents = [solve_cut_patch(face, axial_E_field=drive,
                                    **common).total_current()
                    for drive in (0.0 + 0j, 1.0 + 0j)]
        drive = (1.0 - currents[0]) / (currents[1] - currents[0])
        patch = solve_cut_patch(face, axial_E_field=drive, **common)
        current = patch.total_current()
        resistance = 2.0 * patch.total_loss() / abs(current) ** 2
        k = np.sqrt(-1j * omega * MU0 * sigma)
        exact = float(((k / (2.0 * np.pi * a * sigma))
                       * jv(0, k * a) / jv(1, k * a)).real)
        rows.append({
            "frequency_hz": frequency,
            "radius_over_skin_depth": a / delta,
            "resistance_ohm_per_m": resistance,
            "exact_bessel_ohm_per_m": exact,
            "relative_error": resistance / exact - 1.0,
            "recovered_current_A": abs(current),
        })
        print(json.dumps({"phase": "anchor", **rows[-1]}), flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, default=(
        ROOT / "tests/coil_from_cad/fixtures/beak_fin_straight.step"))
    parser.add_argument("--frequency", type=float, default=150_000.0)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--n-peri", type=int, default=1024)
    parser.add_argument("--maxh-over-delta", type=float, default=0.25)
    parser.add_argument("--order", type=int, default=3)
    parser.add_argument("--control-pad-mm", type=float, nargs="+",
                        default=[8.0, 24.0, 72.0])
    parser.add_argument("--patch-pad-mm", type=float, default=6.0)
    parser.add_argument("--single-cut-mm", type=float, nargs="+",
                        default=[1.0, 0.0, -2.0])
    parser.add_argument("--band-mm", type=float, nargs=2, action="append",
                        default=None,
                        help="outer band as LOW HIGH, repeatable")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    bands = args.band_mm or [(-2.0, 0.0), (-2.5, 0.5), (-1.5, -0.5)]

    from make_beak_fin_step import TIP_RADIUS
    from radia.sibc_corner_patch import solve_cut_patch
    from sibc2d_beak_reference import solve_beak_sibc_2d

    anchor = round_wire_anchor(args.sigma, args.order, args.maxh_over_delta)

    outer = solve_beak_sibc_2d(args.step, frequency=args.frequency,
                               sigma=args.sigma, n_peri=args.n_peri)
    delta = outer["skin_depth_m"]
    xy = np.asarray(outer["panel_xy_m"], dtype=float)
    ds = np.asarray(outer["panel_ds_m"], dtype=float)
    current = np.asarray(outer["panel_current_A"])
    rs = complex(outer["surface_impedance_ohm"]).real
    surface_loss = 0.5 * rs * np.abs(current) ** 2 / ds

    face = section_metres(HERE / "_patch_section_tmp.step")
    common = dict(panel_xy=xy, panel_ds=ds, panel_current=current,
                  axial_E_field=outer["axial_E_field_V_per_m"],
                  surface_impedance=outer["surface_impedance_ohm"],
                  omega=outer["omega_rad_per_s"], sigma=args.sigma,
                  maxh_conductor=args.maxh_over_delta * delta,
                  maxh_air=4.0e-3, order=args.order)

    # --- the control, and the evidence that it is converged ----------------
    controls = []
    for pad_mm in args.control_pad_mm:
        t0 = time.perf_counter()
        solution = solve_cut_patch(face, x_cut=None, air_pad=pad_mm * 1e-3,
                                   **common)
        scale = 1.0 / abs(solution.total_current()) ** 2
        total = solution.total_loss() * scale
        controls.append({
            "air_pad_mm": pad_mm, "ndof": solution.ndof,
            "seconds": time.perf_counter() - t0,
            "recovered_current_A": abs(solution.total_current()),
            "total_loss_W_per_m": total,
            "beak_loss_fraction": solution.loss_beyond(BEAK_X_M) * scale / total,
            "tip_loss_fraction": solution.loss_beyond(TIP_X_M) * scale / total,
        })
        print(json.dumps({"phase": "control", **controls[-1]}), flush=True)
        control, control_scale = solution, scale
    truth_total = controls[-1]["total_loss_W_per_m"]
    truth_beak = control.loss_beyond(BEAK_X_M) * control_scale
    truth_tip = control.loss_beyond(TIP_X_M) * control_scale

    sibc_total = float(np.sum(surface_loss))
    outer_only = {
        "total_loss_W_per_m": sibc_total,
        "total_loss_relative": sibc_total / truth_total - 1.0,
        "beak_loss_fraction": outer["beak_loss_fraction"],
        "beak_relative": outer["beak_loss_fraction"] * truth_total / truth_beak - 1.0,
        "tip_loss_fraction": outer["tip_loss_fraction"],
        "tip_relative": outer["tip_loss_fraction"] * truth_total / truth_tip - 1.0,
    }
    print(json.dumps({"phase": "outer_only", **outer_only}), flush=True)

    # --- one patch, scored on the region it actually owns -------------------
    single = []
    for cut_mm in args.single_cut_mm:
        cut = cut_mm * 1e-3
        t0 = time.perf_counter()
        patch = solve_cut_patch(face, x_cut=cut, side="high",
                                air_pad=args.patch_pad_mm * 1e-3, **common)
        exact_inside = control.loss_beyond(cut) * control_scale
        single.append({
            "x_cut_mm": cut_mm,
            "cut_over_skin_depth_from_root": (BEAK_X_M - cut) / delta,
            "ndof": patch.ndof, "seconds": time.perf_counter() - t0,
            "exact_loss_beyond_cut_W_per_m": exact_inside,
            "patched_relative": patch.loss_beyond(cut) / exact_inside - 1.0,
            "outer_only_relative": float(
                np.sum(surface_loss[xy[:, 0] >= cut])) / exact_inside - 1.0,
            "outer_only_relative_outside_cut": float(
                np.sum(surface_loss[xy[:, 0] < cut]))
            / (truth_total - exact_inside) - 1.0,
        })
        print(json.dumps({"phase": "single", **single[-1]}), flush=True)

    # --- both ends patched, the surface impedance left the flat band --------
    composite = []
    for lo_mm, hi_mm in bands:
        lo, hi = lo_mm * 1e-3, hi_mm * 1e-3
        t0 = time.perf_counter()
        right = solve_cut_patch(face, x_cut=hi, side="high",
                                air_pad=args.patch_pad_mm * 1e-3, **common)
        left = solve_cut_patch(face, x_cut=lo, side="low",
                               air_pad=args.patch_pad_mm * 1e-3, **common)
        band = float(np.sum(surface_loss[(xy[:, 0] >= lo) & (xy[:, 0] < hi)]))
        total = right.loss_beyond(hi) + left.total_loss() + band
        composite.append({
            "outer_band_mm": [lo_mm, hi_mm],
            "outer_band_share": band / total,
            "ndof": right.ndof + left.ndof,
            "seconds": time.perf_counter() - t0,
            "total_loss_W_per_m": total,
            "total_loss_relative": total / truth_total - 1.0,
            "beak_loss_fraction": right.loss_beyond(BEAK_X_M) / total,
            "beak_relative": (right.loss_beyond(BEAK_X_M) / total
                              * truth_total / truth_beak - 1.0),
            "tip_loss_fraction": right.loss_beyond(TIP_X_M) / total,
            "tip_relative": (right.loss_beyond(TIP_X_M) / total
                             * truth_total / truth_tip - 1.0),
        })
        print(json.dumps({"phase": "composite", **composite[-1]}), flush=True)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.sibc_patch_corrected_beak.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": ("the surface impedance fails where the surface curves "
                     "inside a skin depth.  Can it be repaired there instead "
                     "of abandoned, and by how much?"),
        "method": ("overlapping patches, Proekt, Yuferev, Tsukerman and Ida "
                   "(2002); patch scale from the corner-layer result of "
                   "Dauge, Dular, Krahenbuhl, Peron, Perrussel and Poignard "
                   "(2014)"),
        "frequency_hz": args.frequency,
        "sigma_S_per_m": args.sigma,
        "skin_depth_mm": delta * 1e3,
        "tip_radius_mm": float(TIP_RADIUS),
        "skin_depth_over_tip_radius": delta * 1e3 / float(TIP_RADIUS),
        "cuts_mm": {"beak": BEAK_X_M * 1e3, "tip": TIP_X_M * 1e3},
        "control": {
            "what": ("the whole section resolved in the plane, driven by the "
                     "exterior field of the same total current"),
            "anchored_on": ("the exact Bessel round wire, solved by the same "
                            "routine with a closed-form exterior"),
            "round_wire": anchor,
            "levels": controls,
        },
        "surface_impedance_only": outer_only,
        "one_patch": single,
        "two_patches": composite,
        "not_claimed": (
            "this is one section at one frequency, and it is a cost "
            "demonstration for nothing: in two dimensions the patches cover "
            "most of the conductor, so the composite is not cheaper than "
            "resolving all of it.  What is demonstrated is the correction, "
            "not the economy.  The economy is a three-dimensional claim -- a "
            "patch around a feature LINE rather than a whole cross-section -- "
            "and it is not tested here.  End effects are absent by "
            "construction, so these are mid-span quantities"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({"phase": "complete", "output": str(args.output)}))


if __name__ == "__main__":
    main()
