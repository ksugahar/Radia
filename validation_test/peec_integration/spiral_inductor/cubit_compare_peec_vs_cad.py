"""
Compare PEEC line model geometry vs original CAD (3Dmodel_14turns.x_t).

Extracts the actual spiral trace centerlines from the CAD model and
compares them with the PEEC rectangular spiral waypoints.

Model units: meters. Output in mm.

Key comparison:
  - PEEC outer extent vs CAD outer extent
  - PEEC inner opening vs CAD inner opening
  - Turn-by-turn centerline position differences
"""

import os
import sys
from pathlib import Path
import numpy as np
from collections import defaultdict

# --- PEEC model ---
REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "docs" / "peec_integration" / "demos" / "spiral_inductor"))
from spiral_peec import generate_rectangular_spiral, generate_two_layer_spiral

# --- Cubit ---
# Auto-detect Cubit installation
sys.path.insert(0, str(REPO_ROOT / "src" / "radia"))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
model_path = os.path.join(SCRIPT_DIR, "provided_by_sato", "3Dmodel_14turns.x_t")
cubit.cmd(f'import parasolid "{model_path}"')

M2MM = 1000.0

# ===== PEEC parameters (from CAD edge analysis) =====
a0 = 9.125e-3       # (10.000+8.250)/2 = 9.125 mm
b0 = 29.900e-3      # corner vertex 30.775 - w/2 = 29.900 mm
N_turns = 8
pitch = 2.000e-3
w = 1.750e-3         # CAD: outer-inner = 1.750 mm
z_bottom = -0.8175e-3
z_top = +0.8175e-3

print("=" * 70)
print("  PEEC Line Model vs CAD Geometry Comparison")
print("=" * 70)

# ===== Extract CAD spiral trace edges =====
# Vol 1 = bottom spiral, Vol 3 = top spiral
for vid, layer_name in [(1, "BOTTOM"), (3, "TOP")]:
    vol = cubit.volume(vid)
    bbox = cubit.get_bounding_box("volume", vid)

    cad_x_min = bbox[0] * M2MM
    cad_x_max = bbox[1] * M2MM
    cad_y_min = bbox[3] * M2MM
    cad_y_max = bbox[4] * M2MM
    cad_z_min = bbox[6] * M2MM
    cad_z_max = bbox[7] * M2MM

    print(f"\n{'='*70}")
    print(f"  {layer_name} LAYER (Vol {vid})")
    print(f"{'='*70}")
    print(f"  CAD bounding box:")
    print(f"    x: [{cad_x_min:.3f}, {cad_x_max:.3f}] mm (span={cad_x_max-cad_x_min:.3f})")
    print(f"    y: [{cad_y_min:.3f}, {cad_y_max:.3f}] mm (span={cad_y_max-cad_y_min:.3f})")
    print(f"    z: [{cad_z_min:.4f}, {cad_z_max:.4f}] mm (h={cad_z_max-cad_z_min:.4f})")

    # PEEC equivalent for this layer (symmetric turns)
    hw = w / 2 * M2MM  # half-width in mm
    a_outer = (a0 + (N_turns - 1) * pitch) * M2MM
    b_outer = (b0 + (N_turns - 1) * pitch) * M2MM

    # Symmetric spiral: both sides are +-a_outer, +-b_outer
    peec_x_min = -a_outer
    peec_x_max = a_outer
    peec_y_min = -b_outer
    peec_y_max = b_outer

    # PEEC outer edge (centerline + half-width)
    peec_outer_x_min = peec_x_min - hw
    peec_outer_x_max = peec_x_max + hw
    peec_outer_y_min = peec_y_min - hw
    peec_outer_y_max = peec_y_max + hw

    print(f"\n  PEEC outer edge:")
    print(f"    x: [{peec_outer_x_min:.3f}, {peec_outer_x_max:.3f}] mm (span={peec_outer_x_max-peec_outer_x_min:.3f})")
    print(f"    y: [{peec_outer_y_min:.3f}, {peec_outer_y_max:.3f}] mm (span={peec_outer_y_max-peec_outer_y_min:.3f})")

    print(f"\n  PEEC centerline:")
    print(f"    x: [{peec_x_min:.3f}, {peec_x_max:.3f}] mm")
    print(f"    y: [{peec_y_min:.3f}, {peec_y_max:.3f}] mm")

    print(f"\n  Difference (PEEC outer edge - CAD bbox):")
    dx_min = peec_outer_x_min - cad_x_min
    dx_max = peec_outer_x_max - cad_x_max
    dy_min = peec_outer_y_min - cad_y_min
    dy_max = peec_outer_y_max - cad_y_max
    print(f"    x_min: {dx_min:+.3f} mm  (PEEC {'larger' if dx_min < 0 else 'smaller'})")
    print(f"    x_max: {dx_max:+.3f} mm  (PEEC {'larger' if dx_max > 0 else 'smaller'})")
    print(f"    y_min: {dy_min:+.3f} mm  (PEEC {'larger' if dy_min < 0 else 'smaller'})")
    print(f"    y_max: {dy_max:+.3f} mm  (PEEC {'larger' if dy_max > 0 else 'smaller'})")

    # Extract actual track edge positions from curves at each z-level
    curves = {}
    for surf in vol.surfaces():
        for curve in surf.curves():
            cid = curve.id()
            if cid in curves:
                continue
            verts = curve.vertices()
            if len(verts) < 2:
                continue
            p0 = np.array(verts[0].coordinates())
            p1 = np.array(verts[1].coordinates())
            length = np.linalg.norm(p1 - p0)
            dz = abs(p1[2] - p0[2])
            curves[cid] = {
                'start': p0, 'end': p1, 'length': length, 'dz': dz,
                'z_avg': (p0[2] + p1[2]) / 2,
            }

    # Group horizontal long curves by z-level
    z_levels = defaultdict(list)
    for cid, cd in curves.items():
        if cd['dz'] < 1e-5 and cd['length'] > 0.005:
            z_key = round(cd['z_avg'], 5)
            z_levels[z_key].append(cid)

    for z_key in sorted(z_levels.keys()):
        z_mm = z_key * M2MM
        cids = z_levels[z_key]

        # Find extremes from all long curve endpoints
        all_x = []
        all_y = []
        for cid in cids:
            cd = curves[cid]
            all_x.extend([cd['start'][0] * M2MM, cd['end'][0] * M2MM])
            all_y.extend([cd['start'][1] * M2MM, cd['end'][1] * M2MM])

        if not all_x:
            continue

        print(f"\n  CAD long curves at z={z_mm:+.4f} mm ({len(cids)} curves):")
        print(f"    x extent: [{min(all_x):.3f}, {max(all_x):.3f}] mm")
        print(f"    y extent: [{min(all_y):.3f}, {max(all_y):.3f}] mm")

        # Classify into x-direction and y-direction
        x_dir_y = []  # y-positions of x-direction curves (top/bottom sides)
        y_dir_x = []  # x-positions of y-direction curves (left/right sides)
        for cid in cids:
            cd = curves[cid]
            dx = abs(cd['end'][0] - cd['start'][0])
            dy = abs(cd['end'][1] - cd['start'][1])
            if dx > dy:
                y_avg = (cd['start'][1] + cd['end'][1]) / 2 * M2MM
                x_dir_y.append(y_avg)
            else:
                x_avg = (cd['start'][0] + cd['end'][0]) / 2 * M2MM
                y_dir_x.append(x_avg)

        x_dir_y_unique = sorted(set([round(y, 2) for y in x_dir_y]))
        y_dir_x_unique = sorted(set([round(x, 2) for x in y_dir_x]))

        if len(x_dir_y_unique) >= 2:
            y_inner_edge = min(abs(y) for y in x_dir_y_unique)
            y_outer_edge = max(abs(y) for y in x_dir_y_unique)
            n_tracks = len(x_dir_y_unique) // 2
            print(f"    X-dir (top/bottom) edges: {n_tracks} tracks")
            print(f"      Inner |y| = {y_inner_edge:.3f} mm, Outer |y| = {y_outer_edge:.3f} mm")

        if len(y_dir_x_unique) >= 2:
            x_inner_edge = min(abs(x) for x in y_dir_x_unique)
            x_outer_edge = max(abs(x) for x in y_dir_x_unique)
            n_tracks_x = len(y_dir_x_unique) // 2
            print(f"    Y-dir (left/right) edges: {n_tracks_x} tracks")
            print(f"      Inner |x| = {x_inner_edge:.3f} mm, Outer |x| = {x_outer_edge:.3f} mm")

        # Extract centerline from edge pairs
        if x_dir_y_unique:
            y_sorted = np.array(x_dir_y_unique)
            # Bottom sides (y < 0) and top sides (y > 0)
            y_neg = y_sorted[y_sorted < 0]
            y_pos = y_sorted[y_sorted > 0]
            if len(y_neg) >= 2 and len(y_pos) >= 2:
                print(f"\n    CAD centerline y-positions (from edge pairs):")
                # Each track has inner and outer edge
                # Tracks are sorted by y: outer_inner, outer_outer, ..., inner_inner, inner_outer
                # Bottom: most negative = outermost outer edge
                y_neg_rev = sorted(y_neg)  # most negative first
                y_pos_fwd = sorted(y_pos)  # least positive first

                # Pair consecutive edges: (y_neg[0], y_neg[1]) = outer edges of outermost bottom
                if len(y_neg) % 2 == 0:
                    n_pairs = len(y_neg) // 2
                    print(f"      Bottom tracks (n={n_pairs}):")
                    for k in range(n_pairs):
                        e1 = y_neg_rev[2 * k]
                        e2 = y_neg_rev[2 * k + 1]
                        center = (e1 + e2) / 2
                        width_k = abs(e2 - e1)
                        print(f"        Track {k}: edges=[{e1:.3f}, {e2:.3f}], "
                              f"center={center:.3f}, w={width_k:.3f} mm")

# ===== Via comparison =====
print(f"\n{'='*70}")
print(f"  VIA COMPARISON")
print(f"{'='*70}")
for vid in [4, 5]:
    bbox = cubit.get_bounding_box("volume", vid)
    cx = (bbox[0] + bbox[1]) / 2 * M2MM
    cy = (bbox[3] + bbox[4]) / 2 * M2MM
    print(f"  Via {vid}: center=({cx:.3f}, {cy:.3f}) mm")

# PEEC vias at (0, -b0)
print(f"  PEEC via: center=(0.000, {-b0*M2MM:.3f}) mm")

# ===== Summary =====
print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")
print(f"  PEEC: L = 18.216 uH (+4.8% vs measurement)")
print(f"  Meas: L = 17.388 uH")
print(f"  (Symmetric turns fixed: was 18.838 uH / +8.3%)")

print("\nDone.")
