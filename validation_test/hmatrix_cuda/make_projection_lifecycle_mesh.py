"""Give the public tetrahedron's single exterior group a semantic label."""
from pathlib import Path
import hashlib
import importlib.metadata
import json


def generate(source, destination, receipt):
    from netgen.meshing import Mesh
    source, destination, receipt = map(Path, (source, destination, receipt))
    if destination.exists() or source.resolve() == destination.resolve():
        raise ValueError("A new output mesh is required")
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    original_hash = sha(source)
    if original_hash != "6ad1b48b0132e55ccc0fe2930a06fef4f1367c41ee0f9e2bb7c27d48985f8deb":
        raise ValueError("Unexpected public source fixture")
    def inventory(mesh):
        return dict(points=[list(p.p) for p in mesh.Points()],
                    volumes=[(e.index, [v.nr for v in e.vertices]) for e in mesh.Elements3D()],
                    surfaces=[(e.index, [v.nr for v in e.vertices]) for e in mesh.Elements2D()],
                    material=mesh.GetMaterial(1))
    mesh = Mesh(dim=3)
    mesh.Load(str(source))
    before = inventory(mesh)
    if len(before["volumes"]) != 1 or len(before["surfaces"]) != 4 or before["material"] != "air":
        raise ValueError("Expected one public air tetrahedron with four exterior faces")
    if mesh.GetBCName(0) != "boundary" or {e[0] for e in before["surfaces"]} != {1}:
        raise ValueError("Expected one exterior boundary group")
    mesh.SetBCName(0, "exterior")
    mesh.Save(str(destination))
    loaded = Mesh(dim=3)
    loaded.Load(str(destination))
    if inventory(loaded) != before or loaded.GetBCName(0) != "exterior" or sha(source) != original_hash:
        raise ValueError("Mesh geometry, connectivity, material or source changed")
    record = dict(accepted=True, scope="Semantic exterior label only; no solver acceptance",
                  source_sha256=original_hash, output_sha256=sha(destination),
                  length_unit="m", unit_source="Caller declaration; VOL has no unit tag",
                  netgen_version=importlib.metadata.version("netgen-mesher"),
                  geometry_connectivity_material_unchanged=True)
    receipt.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return record


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    generate(args.source, args.destination, args.receipt)
