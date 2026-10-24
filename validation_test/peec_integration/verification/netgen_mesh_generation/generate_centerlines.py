"""Generate native Netgen 1D .vol fixtures from the declared analytic models.

No Gmsh data is read. Coordinates, segment counts and SI parameters follow
these cases' existing public analytic definitions; no volume elements are
fabricated for a conductor centerline.
"""
from pathlib import Path
import hashlib
import json
import math


def parameters(path):
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("#", 1)[0].split()
        if parts:
            result[parts[0]] = float(parts[1])
    return result


def generate(directory=None):
    from netgen.meshing import Mesh, MeshPoint, Element1D, Element0D
    from netgen.csg import Pnt
    directory = Path(directory) if directory else Path(__file__).resolve().parent
    for name, closed in [("wire", False), ("coil", True)]:
        params = parameters(directory / f"circular_{name}_params.txt")
        count = int(params["n_segments"])
        mesh = Mesh(dim=1)
        if closed:
            radius = params["coil_radius"]
            coords = [(radius*math.cos(2*math.pi*i/count),
                       radius*math.sin(2*math.pi*i/count), 0.) for i in range(count)]
        else:
            coords = [(params["length"]*i/count, 0., 0.) for i in range(count+1)]
        points = [mesh.Add(MeshPoint(Pnt(*p))) for p in coords]
        for i in range(count):
            mesh.Add(Element1D([points[i], points[(i+1) % len(points)]], index=1))
        label = f"{name}_centerline"
        mesh.SetMaterial(1, label)
        # The historical analytic fixtures classify every point as centerline.
        # Preserve that label; do not invent electrical port conditions.
        for point in points:
            mesh.Add(Element0D(point, index=1))
        mesh.SetBCName(0, label)
        path = directory / f"circular_{name}_centerline.vol"
        mesh.Save(str(path))
        contract = {"schema": "radia.vol-label-contract.v1",
                    "application": "peec_centerline", "strict_labels": True,
                    "required": {"materials": [label],
                                 "boundaries": [label]},
                    "centerline": {"units": "m", "nodes": len(points),
                                   "segments": count, "closed": closed}}
        path.with_suffix(".labels.json").write_text(json.dumps(contract, indent=2)+"\n", encoding="utf-8", newline="\n")
        print(path.name, hashlib.sha256(path.read_bytes()).hexdigest())


if __name__ == "__main__":
    generate()
