"""
Shared utilities for GmshBuilder package.

Internal helper functions used across modules.
"""

import math


def classify_curves_by_direction(vol_tag):
    """Classify curves of a volume by dominant direction (x, y, z).

    Returns dict {direction_idx: [(curve_tag, length), ...]}
    where direction_idx is 0=x, 1=y, 2=z.
    """
    import gmsh

    surfs = gmsh.model.getBoundary([(3, vol_tag)],
                                    oriented=False, recursive=False)
    curves_set = set()
    for s in surfs:
        edges = gmsh.model.getBoundary([s], oriented=False, recursive=False)
        for e in edges:
            curves_set.add(abs(e[1]))

    classified = {0: [], 1: [], 2: []}

    for ctag in curves_set:
        bounds = gmsh.model.getParametrizationBounds(1, ctag)
        t_min, t_max = bounds[0][0], bounds[1][0]
        p0 = gmsh.model.getValue(1, ctag, [t_min])
        p1 = gmsh.model.getValue(1, ctag, [t_max])

        dx = abs(p1[0] - p0[0])
        dy = abs(p1[1] - p0[1])
        dz = abs(p1[2] - p0[2])
        length = math.sqrt(dx * dx + dy * dy + dz * dz)

        if dx >= dy and dx >= dz:
            classified[0].append((ctag, length))
        elif dy >= dx and dy >= dz:
            classified[1].append((ctag, length))
        else:
            classified[2].append((ctag, length))

    return classified


def apply_transfinite(vol_tag, nx, ny, nz):
    """Apply transfinite constraints to a single volume."""
    import gmsh

    classified = classify_curves_by_direction(vol_tag)

    divisions = {0: nx, 1: ny, 2: nz}
    for direction, curves in classified.items():
        n = divisions[direction]
        for ctag, length in curves:
            gmsh.model.mesh.setTransfiniteCurve(ctag, n + 1)

    surfs = gmsh.model.getBoundary([(3, vol_tag)],
                                    oriented=False, recursive=False)
    for s in surfs:
        try:
            gmsh.model.mesh.setTransfiniteSurface(abs(s[1]))
        except Exception:
            pass

    try:
        gmsh.model.mesh.setTransfiniteVolume(vol_tag)
    except Exception:
        pass

    for s in surfs:
        gmsh.model.mesh.setRecombine(2, abs(s[1]))


def compute_divisions_from_size(vol_tag, mesh_size):
    """Compute nx, ny, nz from mesh size and bounding box."""
    import gmsh

    bb = gmsh.model.occ.getBoundingBox(3, vol_tag)
    lx = bb[3] - bb[0]
    ly = bb[4] - bb[1]
    lz = bb[5] - bb[2]

    nx = max(1, round(lx / mesh_size))
    ny = max(1, round(ly / mesh_size))
    nz = max(1, round(lz / mesh_size))

    return nx, ny, nz


def get_bounding_box_union(vol_tags):
    """Get combined bounding box for multiple GMSH volume tags.

    Returns (xmin, ymin, zmin, xmax, ymax, zmax).
    """
    import gmsh

    xmin, ymin, zmin = float('inf'), float('inf'), float('inf')
    xmax, ymax, zmax = float('-inf'), float('-inf'), float('-inf')
    for vt in vol_tags:
        bb = gmsh.model.occ.getBoundingBox(3, vt)
        xmin = min(xmin, bb[0])
        ymin = min(ymin, bb[1])
        zmin = min(zmin, bb[2])
        xmax = max(xmax, bb[3])
        ymax = max(ymax, bb[4])
        zmax = max(zmax, bb[5])

    return xmin, ymin, zmin, xmax, ymax, zmax


def create_half_space_box(axis, position, bbox, pad):
    """Create a half-space box on the positive side of a cutting plane.

    Returns GMSH OCC volume tag.
    """
    import gmsh

    xmin, ymin, zmin, xmax, ymax, zmax = bbox

    if axis == 'x':
        return gmsh.model.occ.addBox(
            position, ymin - pad, zmin - pad,
            (xmax - position) + pad,
            (ymax - ymin) + 2 * pad,
            (zmax - zmin) + 2 * pad)
    elif axis == 'y':
        return gmsh.model.occ.addBox(
            xmin - pad, position, zmin - pad,
            (xmax - xmin) + 2 * pad,
            (ymax - position) + pad,
            (zmax - zmin) + 2 * pad)
    else:  # z
        return gmsh.model.occ.addBox(
            xmin - pad, ymin - pad, position,
            (xmax - xmin) + 2 * pad,
            (ymax - ymin) + 2 * pad,
            (zmax - position) + pad)


def remove_excess_tools(tool_only_tags):
    """Remove excess tool volumes after fragment operation."""
    import gmsh

    for tag in tool_only_tags:
        try:
            gmsh.model.occ.remove([(3, tag)], recursive=True)
        except Exception:
            pass

    if tool_only_tags:
        gmsh.model.occ.synchronize()


def dimtags_to_tags(dimtags, dim=3):
    """Extract tags of a given dimension from (dim, tag) list."""
    return [dt[1] for dt in dimtags if dt[0] == dim]


def ensure_synchronized():
    """Call gmsh.model.occ.synchronize()."""
    import gmsh
    gmsh.model.occ.synchronize()
