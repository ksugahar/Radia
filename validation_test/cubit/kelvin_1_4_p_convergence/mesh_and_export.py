"""Build Cubit 1/4 sector mesh (x>=0, y>=0, full z) of magnetic sphere
+ air + Kelvin and export Netgen .vol at orders 1, 2, 3.

Reuses `add_kelvin_cubit(reduction={"x":"bn=0","y":"bn=0"}, offset_dir="z")`
which is the verified 1/4 path (matches em_elf_quarter pattern).

Why 1/4 not 1/8: Cubit's copy mesh on a 3-cut octant cap produces
~5% of vertex pairs with O(2 cm) translation error and creates a
hash-table inconsistency that NGSolve cannot rewrap (NgException:
"Ask for unused hash-value").  The blocker is in Cubit's copy-mesh
projection, not in NGSolve.  Documented in memory/feedback_kelvin_1_8_blocker.md.

Geometry (units: meters):
  Magnetic sphere (radius a)                a    = 0.05
  Air outer sphere (kelvin_int boundary)    R    = 0.20
  Kelvin sphere offset (offset_dir="z")    3*R  = 0.60

1/4 sector: x>=0, y>=0, all z
  x = 0 plane:  sym_bn=0_x  (Natural BC for Omega -- dOmega/dx = 0)
  y = 0 plane:  sym_bn=0_y  (Natural BC for Omega -- dOmega/dy = 0)
  z = free      (no reduction along z -- preserves uniform Hz physics)
  GND vertex    Dirichlet (Omega = 0 at kelvin centre = "infinity")

Run:
  python mesh_and_export.py [--out-dir DIR] [--orders 1,2,3]
"""

import argparse
import os
import sys

# --- Boilerplate to attach to in-tree Cubit + plugin ---
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.join(THIS_DIR, "..", "..", "..", "src", "radia")
sys.path.insert(0, SRC_RADIA)
sys.path.insert(0, os.path.join(SRC_RADIA, "panels"))
from install_panels import find_cubit_bin

_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)
_plugin_dir = os.path.join(os.path.dirname(_cubit_path), 'bin', 'plugins') \
    if _cubit_path else None
if _plugin_dir and os.path.isdir(_plugin_dir):
    os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

import cubit
cubit.init(['cubit', '-nojournal', '-batch', '-nographics',
            '-commandplugindir', _plugin_dir or ''])

# Import add_kelvin_cubit AFTER cubit is initialised
from add_kelvin import add_kelvin_cubit


def probe_volumes(label=""):
    print(f"\n--- volumes after: {label} ---")
    for vid in cubit.parse_cubit_list("volume", "all"):
        v = cubit.volume(vid)
        c = v.centroid()
        bb = v.bounding_box()
        print(f"  vol {vid}: c=({c[0]:.3f},{c[1]:.3f},{c[2]:.3f}), "
              f"extent=({bb[3]:.3f},{bb[4]:.3f},{bb[5]:.3f}), "
              f"vol_meas={v.volume():.3e}")


def build_and_export(out_dir, orders=(1, 2, 3),
                     a_mag=0.05, R_air=0.20, mesh_size=0.025):
    """Build geometry, mesh, blocks, kelvin, and export .vol at each order."""
    os.makedirs(out_dir, exist_ok=True)

    cubit.cmd("reset")

    # === Step 1: build magnetic + air spheres, subtract -> mag + air-shell ===
    # Create air sphere FIRST, then mag.  Use `subtract ... keep` then
    # explicitly create a fresh mag sphere that occupies the cavity --
    # because `subtract ... keep` in Cubit leaves the subtrahend
    # (mag) but does NOT actually carve a cavity in the minuend (air).
    # Workaround: subtract WITHOUT keep, then re-create mag at the cavity.
    cubit.cmd(f"create sphere radius {R_air}")
    air_id_initial = cubit.get_last_id("volume")
    cubit.cmd(f"create sphere radius {a_mag}")
    mag_subtrahend = cubit.get_last_id("volume")
    cubit.cmd(f"subtract volume {mag_subtrahend} from volume {air_id_initial}")
    # air now has the mag-shaped cavity; mag was consumed.  Re-create mag.
    cubit.cmd(f"create sphere radius {a_mag}")
    mag_id_initial = cubit.get_last_id("volume")
    probe_volumes("after subtract + recreate mag")

    # === Step 2: webcut to 1/4 sector (x>=0, y>=0, full z); delete neg ===
    cubit.cmd("webcut volume all with plane xplane offset 0")
    cubit.cmd("webcut volume all with plane yplane offset 0")
    drop = []
    for vid in cubit.parse_cubit_list("volume", "all"):
        cx, cy, cz = cubit.volume(vid).centroid()
        if cx < -1e-9 or cy < -1e-9:
            drop.append(vid)
    if drop:
        cubit.cmd("delete volume " + " ".join(str(v) for v in drop))
    probe_volumes("after webcut to 1/4 sector + delete neg")

    # === Step 3: identify magnetic + air octants ===
    survivors = list(cubit.parse_cubit_list("volume", "all"))
    by_meas = sorted(survivors, key=lambda v: cubit.volume(v).volume())
    # subtract...keep on Cubit may leave a duplicate.  Take smallest as
    # magnetic (radius a octant) and largest as air shell octant; drop
    # any duplicates between them.
    mag_vol = by_meas[0]
    air_vol = by_meas[-1]
    extras = [v for v in survivors if v not in (mag_vol, air_vol)]
    if extras:
        print(f"WARNING: deleting extra near-origin volumes {extras} "
              f"(subtract...keep artifact)")
        cubit.cmd("delete volume " + " ".join(str(v) for v in extras))
    print(f"Volume assignment: magnetic={mag_vol}, air={air_vol}")

    # === Step 4: imprint + merge ===
    cubit.cmd("imprint volume all")
    cubit.cmd("merge volume all")

    # === Step 5: mesh interior (magnetic + air shell) ===
    cubit.cmd(f"volume {mag_vol} {air_vol} scheme tetmesh")
    cubit.cmd(f"volume {mag_vol} size {mesh_size * 0.6}")
    cubit.cmd(f"volume {air_vol} size {mesh_size}")
    cubit.cmd(f"mesh volume {mag_vol} {air_vol}")

    # === Step 6: blocks + sphere sideset (mag-air interface) ===
    cubit.cmd("set duplicate block elements on")
    # Sphere sideset: mag-air interface, the curved (non-cut) surface
    # of the magnetic volume.  Probe magnetic volume's surfaces and
    # exclude flat cut faces (zero extent in one coordinate direction).
    mag_surfs_all = list(cubit.get_relatives("volume", mag_vol, "surface"))
    print(f"\nmag_vol={mag_vol} surfaces:")
    sphere_surfs = []
    for sid in mag_surfs_all:
        s = cubit.surface(sid)
        cx, cy, cz = s.center_point()
        bb = s.bounding_box()
        bx, by, bz = bb[3:6]
        is_x_cut = (abs(cx) < 1e-6 and bx < 1e-6)
        is_y_cut = (abs(cy) < 1e-6 and by < 1e-6)
        is_z_cut = (abs(cz) < 1e-6 and bz < 1e-6)
        print(f"  surf {sid}: c=({cx:.3f},{cy:.3f},{cz:.3f}), "
              f"area={s.area():.3e}, "
              f"x_cut={is_x_cut}, y_cut={is_y_cut}, z_cut={is_z_cut}")
        if not (is_x_cut or is_y_cut or is_z_cut):
            sphere_surfs.append(sid)
    print(f"sphere sideset candidates: {sphere_surfs}")
    if sphere_surfs:
        cubit.cmd(f"sideset 100 add surface {' '.join(map(str, sphere_surfs))}")
        cubit.cmd('sideset 100 name "sphere"')
    cubit.cmd(f"block 1 add volume {mag_vol}")
    cubit.cmd('block 1 name "magnetic"')
    cubit.cmd(f"block 2 add volume {air_vol}")
    cubit.cmd('block 2 name "air"')

    # === Step 7: Kelvin sphere via add_kelvin_cubit (1/4 reduction) ===
    add_kelvin_cubit(
        R=R_air, air_block="air",
        reduction={"x": "bn=0", "y": "bn=0"},
        offset_dir="z", offset_dist=3.0 * R_air,
        mesh_size=mesh_size, kelvin_block="kelvin",
        gnd_nodeset=100,
    )

    # === Step 8: export .vol at each order ===
    saved = []
    for order in orders:
        out_path = os.path.join(out_dir, f"sphere_1_4_p{order}.vol").replace(
            "\\", "/")
        cubit.cmd(f'export netgen "{out_path}" order {order} overwrite')
        saved.append(out_path)
        print(f"  Exported order={order}: {out_path}")
    return saved


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", default=THIS_DIR)
    p.add_argument("--orders", default="1,2,3")
    p.add_argument("--mesh-size", type=float, default=0.025)
    args = p.parse_args()
    orders = tuple(int(x) for x in args.orders.split(","))
    files = build_and_export(args.out_dir, orders=orders,
                             mesh_size=args.mesh_size)
    print("\nAll exports:")
    for f in files:
        print(f"  {f}")
