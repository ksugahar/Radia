"""Fix inverted TET node ordering in GMSH .msh v4.1 files (docs-local).

Background
----------
The docs-local animation meshes were exported from NGSolve with the GMSH
tet corners written in ``el.vertices`` order.  NGSolve orders tet vertices
so that ``ref(0,0,0) -> el.vertices[3]``, ``ref(1,0,0) -> el.vertices[0]``,
``ref(0,1,0) -> el.vertices[1]``, ``ref(0,0,1) -> el.vertices[2]`` (see
CLAUDE.md "GMSH API Node Ordering Verification Policy").  In that listing
order ``det[v1-v0, v2-v0, v3-v0] = -det(trafo) < 0``, so every exported
tet is an odd (inverted) permutation of a positive GMSH tet: gmsh
``getJacobians`` reports det <= 0 at every integration point while the
|volume| stays exact.

Fix
---
Rewrite each inverted TET connectivity with the canonical corner mapping
(GMSH node i sits at the same reference coordinates as the NGSolve
reference corner):

    gmsh corner 0 = old corner 3   (ref (0,0,0))
    gmsh corner 1 = old corner 0   (ref (1,0,0))
    gmsh corner 2 = old corner 1   (ref (0,1,0))
    gmsh corner 3 = old corner 2   (ref (0,0,1))

For TET10 the 6 mid-edge nodes are re-associated to the GMSH edge order
(0,1), (1,2), (2,0), (3,0), (3,2), (3,1) over the NEW corners:

    new edge (0,1) = mid(v3,v0) -> old node 7
    new edge (1,2) = mid(v0,v1) -> old node 4
    new edge (2,0) = mid(v1,v3) -> old node 9
    new edge (3,0) = mid(v2,v3) -> old node 8
    new edge (3,2) = mid(v2,v1) -> old node 5
    new edge (3,1) = mid(v2,v0) -> old node 6

Node coordinates, node tags, and every ``$NodeData`` block are untouched:
this is a pure element-connectivity (node order) rewrite.

Only elements whose linear corner determinant is negative are permuted,
so the fix is idempotent and safe on already-correct files.

Usage:
    python tet10_orientation_fix.py mesh1.msh [mesh2.msh ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

# New-node-position -> old-node-position (0-indexed into the 10 node ids)
TET10_PERM = [3, 0, 1, 2, 7, 4, 9, 8, 5, 6]
TET4_PERM = [3, 0, 1, 2]

_TET_PERMS = {4: TET4_PERM, 11: TET10_PERM}


def _corner_det(p0, p1, p2, p3):
    """det[p1-p0, p2-p0, p3-p0] (columns), pure python (no numpy needed)."""
    a = [p1[k] - p0[k] for k in range(3)]
    b = [p2[k] - p0[k] for k in range(3)]
    c = [p3[k] - p0[k] for k in range(3)]
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def _parse_nodes(lines):
    """Return {node_tag: (x, y, z)} from the $Nodes section (v4.1)."""
    coords = {}
    i = lines.index("$Nodes")
    header = lines[i + 1].split()
    n_blocks = int(header[0])
    j = i + 2
    for _ in range(n_blocks):
        blk = lines[j].split()
        n_in_block = int(blk[3])
        j += 1
        tags = [int(lines[j + k]) for k in range(n_in_block)]
        j += n_in_block
        for k in range(n_in_block):
            xyz = lines[j + k].split()
            coords[tags[k]] = (float(xyz[0]), float(xyz[1]), float(xyz[2]))
        j += n_in_block
    return coords


def fix_msh_tet_orientation(path: str | Path) -> dict:
    """Fix inverted TET4/TET10 elements of a .msh v4.1 file in place.

    Returns a dict with counts: {"elements": N, "fixed": M, "kept": N-M}.
    Raises ValueError if the file has no $Elements section.
    """
    p = Path(path)
    data = p.read_bytes()
    # Detect EOL from raw bytes: Path.read_text() universal-newline
    # translation would hide "\r\n" and silently rewrite CRLF files as LF.
    newline = "\r\n" if b"\r\n" in data else "\n"
    text = data.decode("utf-8")
    lines = text.splitlines()

    coords = _parse_nodes(lines)

    try:
        i = lines.index("$Elements")
    except ValueError:
        raise ValueError(f"{p}: no $Elements section")
    header = lines[i + 1].split()
    n_blocks = int(header[0])
    j = i + 2
    n_elements = 0
    n_fixed = 0
    for _ in range(n_blocks):
        blk = lines[j].split()
        el_type = int(blk[2])
        n_in_block = int(blk[3])
        j += 1
        perm = _TET_PERMS.get(el_type)
        if perm is None:
            j += n_in_block
            continue
        for k in range(n_in_block):
            parts = lines[j + k].split()
            tag, nodes = parts[0], [int(v) for v in parts[1:]]
            if len(nodes) != len(perm):
                raise ValueError(
                    f"{p}: element {tag}: expected {len(perm)} nodes "
                    f"for gmsh type {el_type}, got {len(nodes)}")
            det = _corner_det(*(coords[n] for n in nodes[:4]))
            n_elements += 1
            if det < 0.0:
                nodes = [nodes[q] for q in perm]
                lines[j + k] = " ".join([tag] + [str(n) for n in nodes])
                n_fixed += 1
        j += n_in_block

    if n_fixed:
        ends_with_newline = text.endswith("\n")
        out = newline.join(lines) + (newline if ends_with_newline else "")
        p.write_text(out, encoding="utf-8", newline="")
    return {"elements": n_elements, "fixed": n_fixed,
            "kept": n_elements - n_fixed}


def main(argv: list[str]) -> int:
    if not argv:
        print("usage: python tet10_orientation_fix.py <file.msh> [...]")
        return 2
    for f in argv:
        stats = fix_msh_tet_orientation(f)
        print(f"{f}: {stats['fixed']} of {stats['elements']} "
              f"tet elements reoriented ({stats['kept']} already positive)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
