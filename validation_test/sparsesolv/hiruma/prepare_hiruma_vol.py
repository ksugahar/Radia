"""Normalize a Hiruma GMSH/Netgen mesh into a checked Netgen ``.vol``.

The historical GMSH input contains the correct tetrahedral volume topology but
does not emit surface descriptors for the enclosed material interfaces.  This
converter preserves every point and volume element and adds those missing
interface triangles so ``check-vol`` can audit domain ownership.
"""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
from netgen.meshing import Element2D, FaceDescriptor, Mesh
from netgen.read_gmsh import ReadGmsh


def _load(path):
    if path.suffix.lower() == ".msh":
        return ReadGmsh(str(path))
    if path.suffix.lower() == ".vol":
        mesh = Mesh(dim=3)
        mesh.Load(str(path))
        return mesh
    raise ValueError(f"expected .msh or .vol input, got {path}")


def add_missing_material_interfaces(mesh):
    """Add one named surface element for each cross-material tetra face."""
    existing_faces = {
        tuple(sorted(point.nr for point in element.vertices))
        for element in mesh.Elements2D()
    }
    open_faces = {}
    interfaces = []
    for element in mesh.Elements3D():
        vertices = list(element.vertices)
        center = np.mean(
            [np.asarray(mesh[point].p, dtype=float) for point in vertices], axis=0
        )
        domain = int(element.index)
        for face in combinations(vertices, 3):
            key = tuple(sorted(point.nr for point in face))
            previous = open_faces.pop(key, None)
            if previous is None:
                open_faces[key] = (domain, face, center)
                continue
            other_domain, other_face, other_center = previous
            if domain == other_domain or key in existing_faces:
                continue
            if other_domain < domain:
                interfaces.append((other_domain, domain, list(other_face), other_center))
            else:
                interfaces.append((domain, other_domain, list(face), center))

    domain_pairs = sorted({(inside, outside) for inside, outside, _, _ in interfaces})
    next_boundary = mesh.GetNFaceDescriptors() + 1
    descriptor = {}
    for offset, pair in enumerate(domain_pairs):
        boundary = next_boundary + offset
        descriptor[pair] = mesh.Add(
            FaceDescriptor(
                surfnr=-1, domin=pair[0], domout=pair[1], bc=boundary
            )
        )
        mesh.SetBCName(boundary - 1, f"interface_{pair[0]}_{pair[1]}")

    for inside, outside, points, element_center in interfaces:
        xyz = [np.asarray(mesh[point].p, dtype=float) for point in points]
        normal = np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])
        face_center = sum(xyz) / 3.0
        if float(np.dot(normal, element_center - face_center)) < 0.0:
            points[1], points[2] = points[2], points[1]
        mesh.Add(Element2D(descriptor[(inside, outside)], points))

    return {
        "interface_elements_added": len(interfaces),
        "domain_pairs": [list(pair) for pair in domain_pairs],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    if args.source.resolve() == args.output.resolve():
        raise ValueError("source and output must be different files")
    mesh = _load(args.source)
    summary = add_missing_material_interfaces(mesh)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mesh.Save(str(args.output))
    serialized = args.output.read_text(encoding="utf-8")
    lines = [line.rstrip() for line in serialized.splitlines()]
    while lines and not lines[-1]:
        lines.pop()
    args.output.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"saved {args.output}")
    print(f"interface elements added: {summary['interface_elements_added']}")
    print(f"domain pairs: {summary['domain_pairs']}")


if __name__ == "__main__":
    main()
