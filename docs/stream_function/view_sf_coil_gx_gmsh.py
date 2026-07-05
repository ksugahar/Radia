"""Open the Gx fingerprint coil structure in GMSH for visualisation.

Three modes:

  --mode contours (default) regenerate the SF design's equal-current
                            CLOSED CONTOURS only (each fingerprint loop
                            as its own 1D Physical Group, no connection
                            arcs) + the host cylinder.  This is what the
                            stream function method *designs* -- the
                            single-stroke chain's connection arcs are a
                            manufacturability afterthought that is NOT
                            part of the SF solution.
  --mode chain              regenerate the SINGLE-STROKE chain (all
                            contours threaded into one continuous wire
                            with cylinder-geodesic connection arcs).
                            The "filament-ified" coil used by PEEC.
                            Use this to inspect the chain topology
                            (and to see the connection-arc artefacts
                            that contribute the ~22% DSV RMS error).
  --mode step               merge the multi-piece loft-chain STEP file
                            from `demo_sf_to_peec_gx.py --with-peec`.

Per Radia policy: GMSH is used for VISUALISATION only.  This script writes
display files and optionally opens them with a standalone gmsh executable;
it does not import the pip gmsh Python runtime.

Run:  python view_sf_coil_gx_gmsh.py [--mode contours|chain|step] [--no-gui]
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))


def _gmsh_executable():
    """Prefer the lab standalone gmsh.exe; fall back to PATH if needed."""
    preferred = Path(r"C:\gmsh.exe")
    if preferred.exists():
        return str(preferred)
    return shutil.which("gmsh.exe") or shutil.which("gmsh")


def _open_standalone_gmsh(display_geo):
    """Open a .geo display companion with standalone GMSH."""
    exe = _gmsh_executable()
    if exe is None:
        print(f"standalone gmsh executable not found; open manually: {display_geo}")
        return
    subprocess.run([exe, "-noconfig", str(display_geo)], check=False)


def _display_geo_path(source_path):
    source_path = Path(source_path)
    return source_path.with_name(f"{source_path.stem}_gmsh_display.geo")


def _write_display_geo(display_geo, merge_paths, option_lines):
    display_geo = Path(display_geo)
    lines = [
        f"// Open with: gmsh {display_geo.name}",
        "Mesh.NumSubEdges = 4;",
        "General.GraphicsPositionX = 100;",
        "General.GraphicsPositionY = 100;",
        "General.GraphicsWidth = 1280;",
        "General.GraphicsHeight = 800;",
        "General.MenuPositionX = 100;",
        "General.MenuPositionY = 100;",
    ]
    for path in merge_paths:
        lines.append(f'Merge "{Path(path).name.replace("\\\\", "/")}";')
    lines.extend(option_lines)
    display_geo.write_text("\n".join(lines) + "\n", encoding="ascii")
    print(f"wrote display companion {display_geo}")
    return display_geo


def view_step(step_path, show_gui=True):
    """Write a STEP display companion and optionally open standalone GMSH."""
    step_path = Path(step_path)
    display_geo = _display_geo_path(step_path)
    _write_display_geo(
        display_geo,
        [step_path],
        [
            "Geometry.Surfaces = 1;",
            "Geometry.SurfaceLabels = 0;",
        ],
    )
    print(f"STEP display source: {step_path}")
    if show_gui:
        _open_standalone_gmsh(display_geo)


def _sf_design():
    """Run the SF + contour extraction; return (polylines_phi_z, dI, a, L).

    This is the design step before the single-stroke chain -- the *true*
    fingerprint filament structure.  Each polyline is one closed contour
    of psi(phi, z) at one equal-current level, in (phi, z) coordinates.
    """
    from demo_sf_to_peec_gx import (
        loop_corners, make_dsv, _loop_Hz, contour_polylines_phi_z,
    )
    from radia.stream_function import aca_tsvd, pseudo_inverse_solve

    a, L, Gx, dsv = 0.15, 0.50, 1.0, 0.10
    nphi, nz, nlevels = 24, 40, 12
    phi_grid = np.linspace(0.0, 2 * np.pi, nphi, endpoint=False)
    z_grid = np.linspace(-L / 2, L / 2, nz)
    dphi = 2 * np.pi / nphi
    dz = L / (nz - 1)
    corners = [loop_corners(a, pp, zq, dphi, dz)
               for zq in z_grid for pp in phi_grid]
    N = len(corners)
    obs = make_dsv(dsv, 7)
    M = obs.shape[0]
    B = Gx * obs[:, 0]

    def entry(i, j):
        return _loop_Hz(obs[i], corners[j])

    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1e-8)
    psi = pseudo_inverse_solve(res, B, k_mode=res.modes).reshape(nz, nphi)
    polylines, dI, _ = contour_polylines_phi_z(psi, phi_grid, z_grid, nlevels)
    return polylines, dI, a, L


def _phi_z_to_3d(poly_phi_z, a):
    return np.column_stack([
        a * np.cos(poly_phi_z[:, 0]),
        a * np.sin(poly_phi_z[:, 0]),
        poly_phi_z[:, 1],
    ])


def _cylinder_mesh(a, L, n_phi=60, n_z=20):
    cyl_phi = np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False)
    cyl_z = np.linspace(-L / 2, L / 2, n_z)
    nodes = []
    for zz in cyl_z:
        for pp in cyl_phi:
            nodes.append([a * np.cos(pp), a * np.sin(pp), zz])
    return nodes, n_phi, n_z


def build_contours_msh(out_msh, show_cylinder=True):
    """Write a .msh v4.1 with the SF design's closed CONTOUR FILAMENTS only.

    Each contour is its own 1D entity + Physical Group ("contour_NN") so
    GMSH can colour them per-loop and the user can toggle them
    individually.  NO connection arcs (those are the single-stroke
    afterthought; this view is the pure SF design output).
    """
    polylines, dI, a, L = _sf_design()
    # Convert each (phi, z) polyline to 3D and drop duplicates.
    contours_3d = []
    for poly in polylines:
        p3 = _phi_z_to_3d(poly, a)
        diffs = np.linalg.norm(np.diff(p3, axis=0), axis=1)
        keep = np.concatenate([[True], diffs > 1.0e-9])
        p3 = p3[keep]
        if len(p3) >= 2:
            contours_3d.append(p3)
    print(f"SF design: {len(contours_3d)} closed contour filaments, "
          f"dI = {dI:.4g} (Hz units)")

    cyl_nodes, cyl_n_phi, cyl_n_z = (
        _cylinder_mesh(a, L) if show_cylinder else ([], 0, 0)
    )
    n_cyl_nodes = len(cyl_nodes)

    lines = []
    lines.append("$MeshFormat")
    lines.append("4.1 0 8")
    lines.append("$EndMeshFormat")
    n_phys = len(contours_3d) + (1 if show_cylinder else 0)
    lines.append("$PhysicalNames")
    lines.append(str(n_phys))
    for k, _ in enumerate(contours_3d):
        # one physical group per contour for individual visibility / colour
        lines.append(f'1 {k+1} "contour_{k+1:02d}"')
    cyl_pg = len(contours_3d) + 1
    if show_cylinder:
        lines.append(f'2 {cyl_pg} "cylinder_host"')
    lines.append("$EndPhysicalNames")

    # Entities: one 1D per contour + one 2D for cylinder
    n_1d_ent = len(contours_3d)
    n_2d_ent = 1 if show_cylinder else 0
    lines.append("$Entities")
    lines.append(f"0 {n_1d_ent} {n_2d_ent} 0")
    for k, p in enumerate(contours_3d):
        lines.append(
            f"{k+1} {p[:,0].min()} {p[:,1].min()} {p[:,2].min()} "
            f"{p[:,0].max()} {p[:,1].max()} {p[:,2].max()} "
            f"1 {k+1} 0"
        )
    if show_cylinder:
        c = np.array(cyl_nodes)
        lines.append(
            f"1 {c[:,0].min()} {c[:,1].min()} {c[:,2].min()} "
            f"{c[:,0].max()} {c[:,1].max()} {c[:,2].max()} "
            f"1 {cyl_pg} 0"
        )
    lines.append("$EndEntities")

    # Nodes
    n_contour_nodes = sum(len(p) for p in contours_3d)
    n_node_total = n_contour_nodes + n_cyl_nodes
    n_blocks = n_1d_ent + n_2d_ent
    lines.append("$Nodes")
    lines.append(f"{n_blocks} {n_node_total} 1 {n_node_total}")
    node_id = 0
    contour_node_starts = []
    for k, p in enumerate(contours_3d):
        contour_node_starts.append(node_id + 1)
        # dim=1, tag=k+1, parametric=0, n_nodes
        lines.append(f"1 {k+1} 0 {len(p)}")
        for i in range(len(p)):
            node_id += 1
            lines.append(str(node_id))
        for x, y, z in p:
            lines.append(f"{x} {y} {z}")
    if show_cylinder:
        cyl_node_start = node_id + 1
        lines.append(f"2 1 0 {n_cyl_nodes}")
        for i in range(n_cyl_nodes):
            node_id += 1
            lines.append(str(node_id))
        for x, y, z in cyl_nodes:
            lines.append(f"{x} {y} {z}")
    lines.append("$EndNodes")

    # Elements: line2 per contour + triangles for cylinder
    n_contour_segs = sum(len(p) - 1 for p in contours_3d)
    n_cyl_tri = 2 * cyl_n_phi * (cyl_n_z - 1) if show_cylinder else 0
    n_elem_total = n_contour_segs + n_cyl_tri
    lines.append("$Elements")
    lines.append(f"{n_blocks} {n_elem_total} 1 {n_elem_total}")
    eid = 0
    for k, p in enumerate(contours_3d):
        # dim=1, tag=k+1, type=1 (2-node line), n_elements
        lines.append(f"1 {k+1} 1 {len(p) - 1}")
        for i in range(len(p) - 1):
            eid += 1
            lines.append(
                f"{eid} {contour_node_starts[k] + i} "
                f"{contour_node_starts[k] + i + 1}"
            )
    if show_cylinder:
        lines.append(f"2 1 2 {n_cyl_tri}")
        for k in range(cyl_n_z - 1):
            for p_i in range(cyl_n_phi):
                p_next = (p_i + 1) % cyl_n_phi
                n00 = cyl_node_start + k * cyl_n_phi + p_i
                n10 = cyl_node_start + k * cyl_n_phi + p_next
                n01 = cyl_node_start + (k + 1) * cyl_n_phi + p_i
                n11 = cyl_node_start + (k + 1) * cyl_n_phi + p_next
                eid += 1
                lines.append(f"{eid} {n00} {n10} {n11}")
                eid += 1
                lines.append(f"{eid} {n00} {n11} {n01}")
    lines.append("$EndElements")

    Path(out_msh).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out_msh} ({n_node_total} nodes, {n_elem_total} elements,"
          f" {len(contours_3d)} contour filaments"
          f"{' + 2D cylinder' if show_cylinder else ''})")


def build_chain_msh(out_msh, show_cylinder=True):
    """Regenerate the SF chain and write it as a .msh v4.1 file.

    The chain centerline becomes a chain of 1D line2 (GMSH type 1) elements
    in Physical Group "coil_centerline".  Optionally meshes the host
    cylinder surface (radius=0.15m, length=0.5m) as a Physical Group
    "cylinder_host" with type 2 (triangle) elements -- a visual reference
    so you can see where the coil sits on the cylinder.
    """
    from demo_sf_to_peec_gx import (
        loop_corners, make_dsv, _loop_Hz, contour_polylines_phi_z,
        single_stroke_chain_phi_z, single_stroke_lobe_chain_phi_z,
        single_stroke_kuijpers_chain_phi_z, chain_phi_z_to_3d,
    )
    from radia.stream_function import aca_tsvd, pseudo_inverse_solve

    a, L, Gx, dsv = 0.15, 0.50, 1.0, 0.10
    nphi, nz, nlevels = 24, 40, 12
    phi_grid = np.linspace(0.0, 2 * np.pi, nphi, endpoint=False)
    z_grid = np.linspace(-L / 2, L / 2, nz)
    dphi = 2 * np.pi / nphi
    dz = L / (nz - 1)
    corners = [loop_corners(a, pp, zq, dphi, dz)
               for zq in z_grid for pp in phi_grid]
    N = len(corners)
    obs = make_dsv(dsv, 7)
    M = obs.shape[0]
    B = Gx * obs[:, 0]

    def entry(i, j):
        return _loop_Hz(obs[i], corners[j])

    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1e-8)
    psi = pseudo_inverse_solve(res, B, k_mode=res.modes).reshape(nz, nphi)
    polys, _, _ = contour_polylines_phi_z(psi, phi_grid, z_grid, nlevels)
    # default = Kuijpers 2023 Method-1 inspired chain (per-lobe cut line,
    # straight 'rung' blends between contours -- the cleanest of the three
    # methods and the lowest DSV RMS in field comparison).
    chain_phi_z = single_stroke_kuijpers_chain_phi_z(polys, a)
    path = chain_phi_z_to_3d(chain_phi_z, a)
    print(f"chain regenerated (kuijpers method): {len(path)} points, "
          f"{len(path) - 1} segments")

    # Build a .msh v4.1 file from scratch.  Format reference:
    # https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    lines = []
    lines.append("$MeshFormat")
    lines.append("4.1 0 8")
    lines.append("$EndMeshFormat")
    # Physical groups: 1D coil_centerline, 2D cylinder_host
    pg_coil = 1
    pg_cyl = 2
    n_phys = 2 if show_cylinder else 1
    lines.append("$PhysicalNames")
    lines.append(str(n_phys))
    lines.append(f'1 {pg_coil} "coil_centerline"')
    if show_cylinder:
        lines.append(f'2 {pg_cyl} "cylinder_host"')
    lines.append("$EndPhysicalNames")

    # Cylinder host mesh (only if requested): triangulate (phi, z) -> 3D
    cyl_n_phi, cyl_n_z = 60, 20
    cyl_phi = np.linspace(0.0, 2 * np.pi, cyl_n_phi, endpoint=False)
    cyl_z = np.linspace(-L / 2, L / 2, cyl_n_z)
    cyl_nodes = []
    if show_cylinder:
        for zz in cyl_z:
            for pp in cyl_phi:
                cyl_nodes.append([a * np.cos(pp), a * np.sin(pp), zz])
    n_cyl_nodes = len(cyl_nodes)

    # Entities: 1D entity 1 = coil, 2D entity 1 = cylinder
    n_ent = 1 + (1 if show_cylinder else 0)
    lines.append("$Entities")
    lines.append(f"0 1 {1 if show_cylinder else 0} 0")
    # 1D entity: tag=1, bbox, physical_tag
    p = path
    lines.append(
        f"1 {p[:,0].min()} {p[:,1].min()} {p[:,2].min()} "
        f"{p[:,0].max()} {p[:,1].max()} {p[:,2].max()} "
        f"1 {pg_coil} 0"
    )
    if show_cylinder:
        c = np.array(cyl_nodes)
        lines.append(
            f"1 {c[:,0].min()} {c[:,1].min()} {c[:,2].min()} "
            f"{c[:,0].max()} {c[:,1].max()} {c[:,2].max()} "
            f"1 {pg_cyl} 0"
        )
    lines.append("$EndEntities")

    # Nodes: one block for coil, one for cylinder
    n_path = len(path)
    n_node_total = n_path + n_cyl_nodes
    n_blocks = 1 + (1 if show_cylinder else 0)
    lines.append("$Nodes")
    lines.append(f"{n_blocks} {n_node_total} 1 {n_node_total}")
    # Coil block: dim=1, tag=1, parametric=0, n_nodes
    lines.append(f"1 1 0 {n_path}")
    for i in range(n_path):
        lines.append(str(i + 1))
    for x, y, z in path:
        lines.append(f"{x} {y} {z}")
    # Cylinder block: dim=2, tag=1
    if show_cylinder:
        lines.append(f"2 1 0 {n_cyl_nodes}")
        for i in range(n_cyl_nodes):
            lines.append(str(n_path + i + 1))
        for x, y, z in cyl_nodes:
            lines.append(f"{x} {y} {z}")
    lines.append("$EndNodes")

    # Elements: 1D line2 for coil + 2D triangle for cylinder
    n_seg = n_path - 1
    n_cyl_tri = 2 * cyl_n_phi * (cyl_n_z - 1) if show_cylinder else 0
    n_elem_total = n_seg + n_cyl_tri
    lines.append("$Elements")
    lines.append(f"{n_blocks} {n_elem_total} 1 {n_elem_total}")
    # Coil block: dim=1, tag=1, type=1 (2-node line), n_elements
    lines.append(f"1 1 1 {n_seg}")
    for i in range(n_seg):
        lines.append(f"{i + 1} {i + 1} {i + 2}")
    if show_cylinder:
        # Cylinder block: dim=2, tag=1, type=2 (3-node triangle), n_elements
        lines.append(f"2 1 2 {n_cyl_tri}")
        eid = n_seg
        for k in range(cyl_n_z - 1):
            for p_i in range(cyl_n_phi):
                p_next = (p_i + 1) % cyl_n_phi
                n00 = n_path + 1 + k * cyl_n_phi + p_i
                n10 = n_path + 1 + k * cyl_n_phi + p_next
                n01 = n_path + 1 + (k + 1) * cyl_n_phi + p_i
                n11 = n_path + 1 + (k + 1) * cyl_n_phi + p_next
                eid += 1
                lines.append(f"{eid} {n00} {n10} {n11}")
                eid += 1
                lines.append(f"{eid} {n00} {n11} {n01}")
    lines.append("$EndElements")

    Path(out_msh).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out_msh} ({n_node_total} nodes, {n_elem_total} elements,"
          f" 1D coil + {'2D cyl' if show_cylinder else 'no cyl'})")


def view_msh(msh_path, show_gui=True):
    """Write a .msh display companion and optionally open standalone GMSH."""
    msh_path = Path(msh_path)
    display_geo = _display_geo_path(msh_path)
    _write_display_geo(
        display_geo,
        [msh_path],
        [
            "Mesh.Lines = 1;",
            "Mesh.LineWidth = 2.0;",
            "Mesh.SurfaceFaces = 1;",
            "Mesh.SurfaceEdges = 0;",
            "Mesh.ColorCarousel = 0;",
        ],
    )
    if show_gui:
        _open_standalone_gmsh(display_geo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["contours", "chain", "step"],
                    default="contours",
                    help="contours (default) = SF design's separate "
                         "closed fingerprint filaments (no connection "
                         "arcs); chain = single-stroke chain with the "
                         "cylinder-geodesic connection arcs (= what PEEC "
                         "sees); step = merge the loft-chain STEP file")
    ap.add_argument("--step-path", type=Path,
                    default=HERE / "sf_coil_gx.step",
                    help="STEP file to view in --mode step")
    ap.add_argument("--msh-out", type=Path, default=None,
                    help=".msh output (default sf_coil_gx_<mode>.msh)")
    ap.add_argument("--no-cylinder", action="store_true",
                    help="omit the host cylinder reference surface")
    ap.add_argument("--no-gui", action="store_true",
                    help="don't open the GMSH GUI (write file only)")
    args = ap.parse_args()

    if args.mode == "step":
        if not args.step_path.exists():
            print(f"ERROR: {args.step_path} does not exist. Run "
                  f"`python demo_sf_to_peec_gx.py --with-peec` first.")
            sys.exit(1)
        view_step(args.step_path, show_gui=not args.no_gui)
        return

    msh_out = args.msh_out
    if msh_out is None:
        msh_out = HERE / f"sf_coil_gx_{args.mode}.msh"
    if args.mode == "contours":
        build_contours_msh(msh_out, show_cylinder=not args.no_cylinder)
    else:  # chain
        build_chain_msh(msh_out, show_cylinder=not args.no_cylinder)
    view_msh(msh_out, show_gui=not args.no_gui)


if __name__ == "__main__":
    main()
