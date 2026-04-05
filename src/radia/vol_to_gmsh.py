"""
Convert Netgen .vol to GMSH .msh v4.1.

Usage:
  python -m radia.vol_to_gmsh model.vol                    # -> model.msh
  python -m radia.vol_to_gmsh model.vol -o output.msh      # explicit output
  python -m radia.vol_to_gmsh model.vol --boundary         # surface elements only

The .vol file must contain curvedelements for high-order output.
Physical groups are created from material/boundary labels in the .vol.

This is a thin wrapper around GmshPostExport (mesh-only, no field data).
"""

import argparse
import os
import sys


def vol_to_gmsh(vol_path, output=None, boundary=False):
    """Convert .vol to .msh v4.1.

    Args:
        vol_path: Path to Netgen .vol file.
        output: Output .msh path. Default: same base name + .msh.
        boundary: If True, export boundary surface elements only.

    Returns:
        Output .msh path.
    """
    from ngsolve import Mesh
    from radia.gmsh_post_export import GmshPostExport

    if output is None:
        output = os.path.splitext(vol_path)[0] + ".msh"

    mesh = Mesh(vol_path)
    post = GmshPostExport(mesh, boundary=boundary)
    post.write_mesh(output)

    # Write companion .opt for proper GMSH display
    opt_path = output + ".opt"
    with open(opt_path, "w") as f:
        f.write("Mesh.SurfaceFaces = 1;\n")
        f.write("Mesh.VolumeEdges = 0;\n")
        f.write("Mesh.SurfaceEdges = 1;\n")
        f.write("Mesh.NumSubEdges = 4;\n")

    return output


def main():
    parser = argparse.ArgumentParser(
        description="Convert Netgen .vol to GMSH .msh v4.1")
    parser.add_argument("vol", help="Input .vol file")
    parser.add_argument("-o", "--output", help="Output .msh file")
    parser.add_argument("--boundary", action="store_true",
                        help="Export boundary surface elements only")
    args = parser.parse_args()

    if not os.path.exists(args.vol):
        print(f"ERROR: {args.vol} not found", file=sys.stderr)
        sys.exit(1)

    out = vol_to_gmsh(args.vol, args.output, args.boundary)
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
