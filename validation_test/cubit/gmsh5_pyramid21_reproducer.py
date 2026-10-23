"""Native-only Pyramid21 constructor diagnostic; no Radia/Cubit dependencies.

Guards: a metadata-only type lookup cannot masquerade as native element support.
Assumes: one fresh, resource-admitted Gmsh 5 process. The result is an observed
upstream capability, not an xfail, accuracy test or authorization to change scope.
"""
import argparse
import hashlib
from importlib import metadata
import json
from pathlib import Path
import traceback

import gmsh


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    result = {"constructible": False, "distribution_version": metadata.version("gmsh"),
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    gmsh.initialize(["-noconfig"])
    try:
        result["native_version"] = gmsh.option.getString("General.Version")
        result["native_sha256"] = hashlib.sha256(Path(gmsh.lib._name).read_bytes()).hexdigest()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.2)
        name, dimension, order, count, coordinates, primary = gmsh.model.mesh.getElementProperties(125)
        result.update(type=125, name=name, dimension=dimension, order=order, nodes=count)
        gmsh.model.add("pyramid21")
        gmsh.model.addDiscreteEntity(3, 1)
        tags = list(range(1, count + 1))
        gmsh.model.mesh.addNodes(3, 1, tags, coordinates)
        gmsh.model.mesh.addElementsByType(1, 125, [1], tags)
        path = out / "pyramid21.msh"
        gmsh.write(str(path))
        gmsh.clear()
        gmsh.open(str(path))
        element_tags, connectivity = gmsh.model.mesh.getElementsByType(125)
        if list(element_tags) != [1] or list(connectivity) != tags:
            raise RuntimeError("Pyramid21 roundtrip changed connectivity")
        result["constructible"] = True
    except Exception:
        result["error"] = traceback.format_exc()
    finally:
        gmsh.finalize()
        (out / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0 if result["constructible"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
