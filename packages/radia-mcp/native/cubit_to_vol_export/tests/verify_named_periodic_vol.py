"""Guards: native export is readable with exact IDs/type/name in NGSolve.

Assumes: CTest generated the public synthetic fixture via the production C++
matcher. Missing input fails; no runtime Identify or repair is performed.
"""
from pathlib import Path
import sys


def identification_tuple(row):
    if len(row) != 3:
        raise AssertionError("Expected two PointIds and one identification number")
    return tuple(int(p.nr) if hasattr(p, "nr") else int(p) for p in row[:2]) + (int(row[2]),)


def verify_pointid_conversion():
    class PointIdLike:
        def __init__(self, nr):
            self.nr = nr

        def __int__(self):
            raise AssertionError("PointId must be read through nr")

    expected = (1, 7, 1)
    actual = identification_tuple((PointIdLike(1), PointIdLike(7), 1))
    if actual != expected or identification_tuple(expected) != expected:
        raise AssertionError(f"actual={actual}, expected={expected}")


def verify(path):
    from ngsolve import Mesh
    from netgen.meshing import IdentificationType

    mesh = Mesh(str(path))
    actual = {identification_tuple(row) for row in mesh.ngmesh.GetIdentifications()}
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
    verify_pointid_conversion()
    verify(Path(sys.argv[1]))
