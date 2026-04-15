"""
Create induction heating model in Cubit: coil + air + Kelvin exterior domain.

Geometry (Go-Tech toymodel inspired):
  - Coil: single-turn torus (J imposed, not eddy-current solved)
  - Workpiece: cylindrical hole (SIBC Robin BC, not meshed)
  - Air: inner sphere (physical domain)
  - Kelvin: outer spherical shell (truncated exterior, modified nu)

Blocks:
  coil        - volume (HCurl source region)
  air         - volume (may be multiple volumes due to boolean splits)
  kelvin      - volume (modified nu = nu0*(r/a)^4)
  wp_surface  - surface (Robin BC)
  outer       - surface (Dirichlet A=0)

Usage:
  python create_induction_model.py                    # default toymodel
  python create_induction_model.py --mesh-size 0.005  # finer mesh
  # Or in Cubit GUI: play "create_induction_model.py"
"""

import math
import os
import sys

# ============================================================
# Parameters (Go-Tech toymodel inspired)
# ============================================================
R_COIL = 0.030       # Coil center (major) radius [m]
A_COIL = 0.005       # Coil wire (minor) radius [m]
R_WP = 0.015         # Workpiece cylinder radius [m]
H_WP = 0.040         # Workpiece cylinder height [m]
R_AIR = 0.100        # Air sphere / Kelvin boundary radius [m]
R_KELVIN = 0.200     # Kelvin exterior domain outer radius [m] (Dirichlet at physical infinity)
MESH_SIZE = 0.008    # Global mesh size [m]
MESH_SIZE_COIL = 0.004  # Coil mesh size [m]


def _get_vols(cubit):
    return set(cubit.get_entities("volume"))


def _classify_volumes(cubit, R_coil, a_coil, R_wp, H_wp, R_air, R_kelvin):
    """Classify all existing volumes by geometry (centroid + volume size).

    Cubit boolean operations create volumes with unpredictable IDs.
    This function inspects each volume to determine its role.

    Returns:
        dict mapping volume ID -> role ('coil', 'air', 'kelvin', 'workpiece')
    """
    coil_vol_expected = math.pi * a_coil**2 * 2 * math.pi * R_coil
    kelvin_vol_expected = 4/3 * math.pi * (R_kelvin**3 - R_air**3)
    wp_vol_expected = math.pi * R_wp**2 * H_wp

    classified = {}
    for vid in sorted(cubit.get_entities("volume")):
        v = cubit.volume(vid)
        vol = abs(v.volume())
        cx, cy, cz = v.center_point()
        r_cen = math.sqrt(cx**2 + cy**2)
        d_cen = math.sqrt(cx**2 + cy**2 + cz**2)

        # Coil: volume matches torus, centroid at radius ~R_coil
        if (abs(vol - coil_vol_expected) / max(coil_vol_expected, 1e-30) < 0.5
                and abs(r_cen - R_coil) < R_coil * 0.5):
            classified[vid] = 'coil'
        # Kelvin exterior domain: very large volume, centroid near origin
        elif vol > kelvin_vol_expected * 0.3 and d_cen < R_kelvin * 0.3:
            classified[vid] = 'kelvin'
        # Workpiece: small, centroid near origin, matches cylinder
        elif (abs(vol - wp_vol_expected) / max(wp_vol_expected, 1e-30) < 0.5
              and d_cen < R_wp * 1.5):
            classified[vid] = 'workpiece'
        else:
            # Everything else is air
            classified[vid] = 'air'

    return classified


def create_induction_model(cubit,
                           R_coil=R_COIL, a_coil=A_COIL,
                           R_wp=R_WP, H_wp=H_WP,
                           R_air=R_AIR, R_kelvin=R_KELVIN,
                           mesh_size=MESH_SIZE, mesh_size_coil=MESH_SIZE_COIL):
    """Create induction heating model in Cubit.

    After boolean operations, volumes are classified by geometry (centroid
    and volume size) rather than tracking IDs, since Cubit may create/destroy
    volumes unpredictably.

    Returns:
        dict with volume IDs, surface IDs, geometry params
    """
    cubit.cmd('reset')
    cubit.cmd('set journal off')

    print("=" * 60)
    print("Creating induction heating model")
    print("=" * 60)
    print(f"  Coil:      R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm")
    print(f"  Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm (hole)")
    print(f"  Air:       R={R_air*1e3:.0f}mm")
    print(f"  Kelvin:    R={R_kelvin*1e3:.0f}mm")

    # ============================================================
    # Step 1: Kelvin exterior domain (optional: outer sphere - inner sphere)
    # ============================================================
    use_kelvin = R_kelvin > R_air
    if use_kelvin:
        vb = _get_vols(cubit)
        cubit.cmd(f'create sphere radius {R_kelvin}')
        outer_id = max(_get_vols(cubit) - vb)
        vb = _get_vols(cubit)
        cubit.cmd(f'create sphere radius {R_air}')
        inner_id = max(_get_vols(cubit) - vb)
        cubit.cmd(f'subtract volume {inner_id} from volume {outer_id}')
    shell_vols = _get_vols(cubit)
    print(f"  After shell: volumes = {sorted(shell_vols)}"
          + (" (Kelvin)" if use_kelvin else " (no Kelvin)"))

    # ============================================================
    # Step 2: Air sphere with workpiece hole
    # ============================================================
    vb = _get_vols(cubit)
    cubit.cmd(f'create sphere radius {R_air}')
    air_id = max(_get_vols(cubit) - vb)

    vb = _get_vols(cubit)
    cubit.cmd(f'create cylinder height {H_wp} radius {R_wp}')
    wp_id = max(_get_vols(cubit) - vb)

    cubit.cmd(f'subtract volume {wp_id} from volume {air_id}')
    # air may have new ID
    air_after_wp = _get_vols(cubit) - shell_vols
    print(f"  After wp hole: air candidates = {sorted(air_after_wp)}")

    # ============================================================
    # Step 3: Remove coil shape from air
    # ============================================================
    vb = _get_vols(cubit)
    cubit.cmd(f'create torus major radius {R_coil} minor radius {a_coil}')
    coil_tool_id = max(_get_vols(cubit) - vb)

    # Find the air volume (largest non-shell, non-coil-tool volume)
    air_candidates = sorted(_get_vols(cubit) - shell_vols - {coil_tool_id})
    air_id = max(air_candidates) if air_candidates else air_id
    cubit.cmd(f'subtract volume {coil_tool_id} from volume {air_id}')
    print(f"  After coil subtract: volumes = {sorted(_get_vols(cubit))}")

    # ============================================================
    # Step 4: Recreate coil torus (consumed by subtract)
    # ============================================================
    vb = _get_vols(cubit)
    cubit.cmd(f'create torus major radius {R_coil} minor radius {a_coil}')
    coil_id = max(_get_vols(cubit) - vb)
    print(f"  After coil recreate: volumes = {sorted(_get_vols(cubit))}")
    print(f"  Coil vol ID (explicit) = {coil_id}")

    # ============================================================
    # Step 5: Classify all volumes by geometry
    # ============================================================
    classified = _classify_volumes(
        cubit, R_coil, a_coil, R_wp, H_wp, R_air, R_kelvin)
    # Override coil: torus centroid is at origin, classifier may confuse
    # it with workpiece. We know the exact ID from creation.
    classified[coil_id] = 'coil'

    # Collect by role
    coil_vols = [v for v, r in classified.items() if r == 'coil']
    air_vols = [v for v, r in classified.items() if r == 'air']
    kelvin_vols = [v for v, r in classified.items() if r == 'kelvin']
    wp_vols = [v for v, r in classified.items() if r == 'workpiece']

    print(f"\n  Classification:")
    for vid, role in sorted(classified.items()):
        v = cubit.volume(vid)
        vol = abs(v.volume())
        cx, cy, cz = v.center_point()
        print(f"    vol {vid}: {role:10s} (V={vol:.4e}, "
              f"centroid=({cx:.4f},{cy:.4f},{cz:.4f}))")

    if not coil_vols:
        print("  ERROR: No coil volume found!")
        return None
    if not air_vols:
        print("  ERROR: No air volume found!")
        return None

    # Delete workpiece volumes (they should be holes)
    for vid in wp_vols:
        cubit.cmd(f'delete volume {vid}')
        print(f"  Deleted workpiece volume {vid}")

    # ============================================================
    # Step 3: Imprint and merge for conformal interfaces
    # ============================================================
    all_mesh_vols = coil_vols + air_vols + kelvin_vols
    vol_str = " ".join(str(v) for v in all_mesh_vols)
    cubit.cmd(f'imprint volume {vol_str}')
    cubit.cmd(f'merge volume {vol_str}')

    # ============================================================
    # Step 4: Mesh
    # ============================================================
    for vid in air_vols:
        cubit.cmd(f'volume {vid} scheme tetmesh')
        cubit.cmd(f'volume {vid} size {mesh_size}')
    for vid in coil_vols:
        cubit.cmd(f'volume {vid} scheme tetmesh')
        cubit.cmd(f'volume {vid} size {mesh_size_coil}')
    for vid in kelvin_vols:
        cubit.cmd(f'volume {vid} scheme tetmesh')
        cubit.cmd(f'volume {vid} size {mesh_size * 2}')

    cubit.cmd(f'mesh volume {vol_str}')

    ne_total = 0
    for vid in all_mesh_vols:
        ne = len(cubit.get_volume_tets(vid))
        ne_total += ne

    print(f"\n  Mesh: {ne_total} tets total")

    # ============================================================
    # Step 5: Identify boundary surfaces
    # ============================================================
    wp_surfaces = []
    outer_surfaces = []

    # Workpiece surface: free surfaces of air near the cylinder hole
    for vid in air_vols:
        for sid in cubit.get_relatives("volume", vid, "surface"):
            adj_vols = cubit.get_relatives("surface", sid, "volume")
            if len(adj_vols) > 1:
                continue  # shared surface
            s = cubit.surface(sid)
            cx, cy, cz = s.center_point()
            r_cen = math.sqrt(cx**2 + cy**2)
            d_cen = math.sqrt(cx**2 + cy**2 + cz**2)

            if d_cen > R_air * 0.7:
                continue  # near outer boundary, not wp

            # Near workpiece cylinder
            is_lateral = (abs(r_cen - R_wp) < R_wp * 0.4
                          and abs(cz) <= H_wp / 2 + 0.002)
            is_cap = (r_cen < R_wp * 1.2
                      and abs(abs(cz) - H_wp / 2) < H_wp * 0.15)
            if is_lateral or is_cap:
                wp_surfaces.append(sid)

    # Outer boundary: free surfaces at outermost radius
    outer_radius = R_kelvin if use_kelvin else R_air
    outer_search_vols = kelvin_vols if kelvin_vols else air_vols
    for vid in outer_search_vols:
        for sid in cubit.get_relatives("volume", vid, "surface"):
            adj_vols = cubit.get_relatives("surface", sid, "volume")
            if len(adj_vols) > 1:
                continue
            s = cubit.surface(sid)
            cx, cy, cz = s.center_point()
            d_cen = math.sqrt(cx**2 + cy**2 + cz**2)
            if d_cen > outer_radius * 0.7:
                outer_surfaces.append(sid)

    print(f"  Workpiece surfaces: {wp_surfaces}")
    print(f"  Outer surfaces: {outer_surfaces}")

    if not wp_surfaces:
        print("  WARNING: No workpiece surfaces detected!")
    if not outer_surfaces:
        print("  WARNING: No outer surfaces detected!")

    # ============================================================
    # Step 6: Define blocks
    # ============================================================
    cubit.cmd('set duplicate block elements on')

    # Volume blocks
    for vid in coil_vols:
        cubit.cmd(f'block 1 add volume {vid}')
    cubit.cmd('block 1 name "coil"')

    for vid in air_vols:
        cubit.cmd(f'block 2 add volume {vid}')
    cubit.cmd('block 2 name "air"')

    for vid in kelvin_vols:
        cubit.cmd(f'block 3 add volume {vid}')
    cubit.cmd('block 3 name "kelvin"')

    # Surface blocks
    for sid in wp_surfaces:
        tris = cubit.get_surface_tris(sid)
        if tris:
            cubit.cmd(f'block 4 add tri in surface {sid}')
    if wp_surfaces:
        cubit.cmd('block 4 name "wp_surface"')

    for sid in outer_surfaces:
        tris = cubit.get_surface_tris(sid)
        if tris:
            cubit.cmd(f'block 5 add tri in surface {sid}')
    if outer_surfaces:
        cubit.cmd('block 5 name "outer"')

    info = {
        'coil_vols': coil_vols, 'air_vols': air_vols,
        'kelvin_vols': kelvin_vols,
        'wp_surfaces': wp_surfaces, 'outer_surfaces': outer_surfaces,
        'R_coil': R_coil, 'a_coil': a_coil,
        'R_wp': R_wp, 'H_wp': H_wp,
        'R_air': R_air, 'R_kelvin': R_kelvin,
        'ne_total': ne_total,
    }
    print(f"  Total elements: {ne_total}")
    print("  Done.")
    return info


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    import argparse
    try:
        parser = argparse.ArgumentParser(description="Create induction model")
        parser.add_argument("--mesh-size", type=float, default=MESH_SIZE)
        parser.add_argument("--r-air", type=float, default=R_AIR)
        parser.add_argument("--r-kelvin", type=float, default=R_KELVIN)
        args, _ = parser.parse_known_args()
    except:
        class Args:
            mesh_size = MESH_SIZE
            r_air = R_AIR
            r_kelvin = R_KELVIN
        args = Args()

    cubit_path = os.environ.get("CUBIT_PATH")
    if cubit_path:
        sys.path.append(cubit_path)

    import cubit
    try:
        cubit.init(['cubit', '-nojournal', '-batch'])
    except:
        pass

    info = create_induction_model(
        cubit, R_air=args.r_air, R_kelvin=args.r_kelvin,
        mesh_size=args.mesh_size)

    script_dir = os.path.dirname(os.path.abspath(__file__)) \
        if '__file__' in dir() else os.getcwd()
    cub5 = os.path.join(script_dir, 'induction_model.cub5')
    cubit.cmd(f'save cub5 "{cub5}" overwrite')
    print(f"\nSaved: {cub5}")
