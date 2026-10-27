"""Guards: native export is readable with exact IDs/type/name in NGSolve.

Assumes: CTest generated the public synthetic fixture via the production C++
matcher. Missing input fails; no runtime Identify or repair is performed.
"""
from pathlib import Path
import sys


def verify(path):
    from ngsolve import Mesh
    from netgen.meshing import IdentificationType

    mesh = Mesh(str(path))
    actual = {tuple(map(int, row)) for row in mesh.ngmesh.GetIdentifications()}
    expected = {(i, i + 6, 1) for i in range(1, 7)}
    if mesh.ngmesh.GetNrIdentifications() != 1 or actual != expected:
        raise AssertionError(f"actual group/pairs={actual}, expected={expected}")
    if mesh.ngmesh.GetIdentificationType(1) != IdentificationType.PERIODIC:
        raise AssertionError("Native periodic type differs")
    lines = path.read_text(encoding="utf-8").splitlines()
    sections = [i for i, value in enumerate(lines) if value == "identificationnames"]
    if len(sections) != 1 or lines[sections[0]+1:sections[0]+3] != ["1", "kt"]:
        raise AssertionError("Native persisted name must be the sole kt group")


if __name__ == "__main__":
    verify(Path(sys.argv[1]))
