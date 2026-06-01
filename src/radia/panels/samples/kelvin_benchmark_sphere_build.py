"""Build a Cubit-meshed magnetic-sphere + Kelvin model for the EM panel
Kelvin Benchmark mode at any of the three supported symmetry fractions
(1/2, 1/4, 1/8).  Exports Netgen `.vol` at orders 1, 2, 3.

Geometry (units: meters):
  Magnetic sphere (radius a)               a    = 0.05
  Air outer sphere (kelvin_int boundary)   R    = 0.20
  Kelvin sphere offset                     3*R  = 0.60

Symmetry fractions and BC choice (uniform Hz field, field_axis = z):

  --frac 1_2 : reduction = {"y": "bn=0"}                offset along z
               (one symmetry plane: y=0; B is parallel to y=0 plane.)

  --frac 1_4 : reduction = {"x": "bn=0", "y": "bn=0"}   offset along z
               (two symmetry planes: x=0 and y=0; B parallel to both.)

  --frac 1_8 : reduction = {"x": "bn=0",
                            "y": "bn=0",
                            "z": "ht=0"}                offset along x
               (three planes; z=0 is the field-perpendicular plane,
               so it must be ht=0 -- bn=0 there would force B = 0.
               offset_dir cannot be a free axis (none exists), so the
               kelvin_far face becomes the Dirichlet "infinity" plane.)

Kelvin-far singularity (1/8 only):
  The reluctivity Mu = mu0 * (R/r')^2 used in the Kelvin exterior is
  singular at r' = 0 (the Kelvin sphere centre).  In the 1/2 and 1/4
  cases the offset axis is free (perpendicular to all reduction
  planes), so r'=0 sits inside the Kelvin volume, away from the
  Dirichlet boundary.  In the 1/8 case the offset axis IS one of the
  reduction axes, which means the kelvin_far cut plane passes
  through r'=0 -- mesh elements there have a 1/0 reluctivity at the
  Dirichlet boundary, leading to ill-conditioning that the standard
  H1 + pardiso linear solve cannot handle on this mesh density.
  The 1/8 build still produces a clean .vol; the solver issue is a
  formulation problem, not a mesh problem.  See README.

Run:
  python kelvin_benchmark_sphere_build.py --frac 1_4 [--out-dir DIR]
                                          [--orders 1,2,3] [--mesh-size 0.025]
"""

import argparse
import os
import sys


# --- Boilerplate to attach to in-tree Cubit + plugin ---
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.join(THIS_DIR, "..", "..")
sys.path.insert(0, SRC_RADIA)
sys.path.insert(0, os.path.join(SRC_RADIA, "panels"))


def _init_cubit():
    """Initialise Cubit batch + plugin path; idempotent."""
    from install_panels import find_cubit_bin
    cubit_path = find_cubit_bin()
    if cubit_path and cubit_path not in sys.path:
        sys.path.append(cubit_path)
    plugin_dir = os.path.join(os.path.dirname(cubit_path), 'bin', 'plugins') \
        if cubit_path else None
    if plugin_dir and os.path.isdir(plugin_dir):
        os.environ['CUBIT_PLUGIN_DIR'] = plugin_dir

    import cubit
    if not getattr(cubit, "_initialized_for_kelvin_build", False):
        cubit.init(['cubit', '-nojournal', '-batch', '-nographics',
                    '-commandplugindir', plugin_dir or ''])
        cubit._initialized_for_kelvin_build = True
    return cubit


# Per-fraction parameter table.
#
# axes         : list of reduction axes (order matters: matches webcut order)
# drop_axes    : list of axes whose negative side is dropped (always == axes)
# reduction    : reduction dict passed to add_kelvin_cubit
# offset_dir   : Kelvin-sphere offset direction (must be a free axis for
#                1/2 and 1/4; for 1/8 must be one of the reduction axes)
# field_axis   : default uniform-field direction for the benchmark
_FRAC_TABLE = {
    "1_2": {
        "axes":       ["y"],
        "reduction":  {"y": "bn=0"},
        "offset_dir": "z",
        "field_axis": "z",
        "label":      "1/2",
    },
    "1_4": {
        "axes":       ["x", "y"],
        "reduction":  {"x": "bn=0", "y": "bn=0"},
        "offset_dir": "z",
        "field_axis": "z",
        "label":      "1/4",
    },
    "1_8": {
        "axes":       ["x", "y", "z"],
        "reduction":  {"x": "bn=0", "y": "bn=0", "z": "ht=0"},
        "offset_dir": "x",
        "field_axis": "z",
        "label":      "1/8",
    },
}


def _probe_volumes(cubit, label=""):
    print(f"\n--- volumes after: {label} ---")
    for vid in cubit.parse_cubit_list("volume", "all"):
        v = cubit.volume(vid)
        c = v.centroid()
        bb = v.bounding_box()
        print(f"  vol {vid}: c=({c[0]:.3f},{c[1]:.3f},{c[2]:.3f}), "
              f"extent=({bb[3]:.3f},{bb[4]:.3f},{bb[5]:.3f}), "
              f"vol_meas={v.volume():.3e}")


def build_and_export(frac, out_dir, orders=(1, 2, 3),
                     a_mag=0.05, R_air=0.20, mesh_size=0.025,
                     filename_prefix=None):
    """Build geometry, mesh, blocks, kelvin and export .vol at each order.

    Args:
        frac: One of "1_2", "1_4", "1_8".
        out_dir: Directory to write `<prefix>_p{order}.vol`.
        orders: Tuple of polynomial orders to export.
        a_mag: Magnetic sphere radius [m].
        R_air: Air outer sphere radius [m] (= kelvin_int boundary).
        mesh_size: Tet edge length in the air shell [m]. Magnetic
            volume gets 0.6x this (denser inside the sphere).
        filename_prefix: Output basename prefix.  Defaults to
            `sphere_<frac>` (e.g. `sphere_1_4_p2.vol`).
    """
    if frac not in _FRAC_TABLE:
        raise ValueError(
            f"frac must be one of {sorted(_FRAC_TABLE)}, got {frac!r}")
    params = _FRAC_TABLE[frac]
    axes = params["axes"]

    cubit = _init_cubit()
    from add_kelvin import add_kelvin_cubit

    os.makedirs(out_dir, exist_ok=True)
    if filename_prefix is None:
        filename_prefix = f"sphere_{frac}"

    cubit.cmd("reset")

    # === Step 1: build air + magnetic spheres, subtract -> air-shell + mag ===
    # In Cubit 2025.3, `subtract ... keep` does NOT actually carve the
    # minuend.  Workaround: subtract WITHOUT `keep`, then re-create mag
    # at the (now-empty) cavity.
    cubit.cmd(f"create sphere radius {R_air}")
    air_id_initial = cubit.get_last_id("volume")
    cubit.cmd(f"create sphere radius {a_mag}")
    mag_subtrahend = cubit.get_last_id("volume")
    cubit.cmd(f"subtract volume {mag_subtrahend} from volume {air_id_initial}")
    cubit.cmd(f"create sphere radius {a_mag}")
    _probe_volumes(cubit, "after subtract + recreate mag")

    # === Step 2: webcut on each reduction axis; drop negative-side ===
    for axis in axes:
        cubit.cmd(f"webcut volume all with plane {axis}plane offset 0")
    drop = []
    axis_idx = {"x": 0, "y": 1, "z": 2}
    for vid in cubit.parse_cubit_list("volume", "all"):
        c = cubit.volume(vid).centroid()
        for axis in axes:
            if c[axis_idx[axis]] < -1e-9:
                drop.append(vid)
                break
    if drop:
        cubit.cmd("delete volume " + " ".join(str(v) for v in drop))
    _probe_volumes(cubit, f"after webcut to {params['label']} sector + delete neg")

    # === Step 3: identify mag (smallest) + air (largest) survivors ===
    survivors = list(cubit.parse_cubit_list("volume", "all"))
    by_meas = sorted(survivors, key=lambda v: cubit.volume(v).volume())
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

    # === Step 5: mesh interior (mag + air shell) ===
    cubit.cmd(f"volume {mag_vol} {air_vol} scheme tetmesh")
    cubit.cmd(f"volume {mag_vol} size {mesh_size * 0.6}")
    cubit.cmd(f"volume {air_vol} size {mesh_size}")
    cubit.cmd(f"mesh volume {mag_vol} {air_vol}")

    # === Step 6: blocks + sphere sideset (mag-air interface) ===
    # The mag-air interface is the curved surface of the mag volume
    # (the non-cut face).  All other mag-volume surfaces are flat cut
    # faces in one of the reduction planes.
    cubit.cmd("set duplicate block elements on")
    mag_surfs_all = list(cubit.get_relatives("volume", mag_vol, "surface"))
    print(f"\nmag_vol={mag_vol} surfaces:")
    sphere_surfs = []
    for sid in mag_surfs_all:
        s = cubit.surface(sid)
        cx, cy, cz = s.center_point()
        bb = s.bounding_box()
        bx, by, bz = bb[3:6]
        coords = (cx, cy, cz)
        extents = (bx, by, bz)
        cut_axes = []
        for axis in axes:
            ai = axis_idx[axis]
            if abs(coords[ai]) < 1e-6 and extents[ai] < 1e-6:
                cut_axes.append(axis)
        is_cut = bool(cut_axes)
        print(f"  surf {sid}: c=({cx:.3f},{cy:.3f},{cz:.3f}), "
              f"area={s.area():.3e}, cut_on={cut_axes if cut_axes else '-'}")
        if not is_cut:
            sphere_surfs.append(sid)
    print(f"sphere sideset candidates: {sphere_surfs}")
    if sphere_surfs:
        cubit.cmd(f"sideset 100 add surface {' '.join(map(str, sphere_surfs))}")
        cubit.cmd('sideset 100 name "sphere"')
    cubit.cmd(f"block 1 add volume {mag_vol}")
    cubit.cmd('block 1 name "magnetic"')
    cubit.cmd(f"block 2 add volume {air_vol}")
    cubit.cmd('block 2 name "air"')

    # === Step 7: Kelvin sphere via add_kelvin_cubit (reduction mode) ===
    add_kelvin_cubit(
        R=R_air, air_block="air",
        reduction=params["reduction"],
        offset_dir=params["offset_dir"], offset_dist=3.0 * R_air,
        mesh_size=mesh_size, kelvin_block="kelvin",
        gnd_nodeset=100,
    )

    # === Step 8: export .vol at each order ===
    saved = []
    for order in orders:
        out_path = os.path.join(
            out_dir, f"{filename_prefix}_p{order}.vol").replace("\\", "/")
        cubit.cmd(f'export netgen "{out_path}" order {order} overwrite')
        saved.append(out_path)
        print(f"  Exported order={order}: {out_path}")
    return saved


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--frac", required=True, choices=sorted(_FRAC_TABLE),
                   help="Symmetry fraction (1/2, 1/4, 1/8)")
    p.add_argument("--out-dir", default=THIS_DIR,
                   help="Directory to write the .vol files")
    p.add_argument("--orders", default="1,2,3",
                   help="Comma-separated polynomial orders to export")
    p.add_argument("--mesh-size", type=float, default=0.025,
                   help="Tet edge length in the air shell [m]")
    p.add_argument("--prefix", default=None,
                   help="Output filename prefix (defaults to sphere_<frac>)")
    args = p.parse_args(argv)
    orders = tuple(int(x) for x in args.orders.split(","))
    files = build_and_export(args.frac, args.out_dir, orders=orders,
                             mesh_size=args.mesh_size,
                             filename_prefix=args.prefix)
    print("\nAll exports:")
    for f in files:
        print(f"  {f}")


if __name__ == "__main__":
    main()
