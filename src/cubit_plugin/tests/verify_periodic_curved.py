"""Guards: the compact writer's paired q1/q2/q3 faces share a translated map.

Assumes: the native synthetic fixture uses dyadic coefficients and reference
points, making these exact arithmetic checks independent of CAD tolerances.
Full Netgen reads the file; no Curve, Save, Identify or coefficient repair.
"""
from pathlib import Path
import sys


def verify(root):
    import ngsolve as ng
    from netgen.meshing import IdentificationType

    ng.SetNumThreads(1)
    xyz = ng.CF((ng.x, ng.y, ng.z))
    # The first reference point is on vertices 1/2/3; the second is on
    # vertices 2/3/4. Both exercise their shared edge and q3 face bubbles.
    samples = [(0.25, 0.25, 0.5), (0.0, 0.25, 0.25)]
    for suffix, mapping in (("", (8, 7, 6, 5)), ("-mixed", (5, 7, 6, 8))):
      for order in (1, 2, 3):
        mesh = ng.Mesh(str(root / f"periodic-curved-q{order}{suffix}.vol"))
        actual = {(int(a.nr), int(b.nr), int(nr))
                  for a, b, nr in mesh.ngmesh.GetIdentifications()}
        expected = {(i + 1, mapping[i], 1) for i in range(4)}
        if actual != expected or mesh.ngmesh.GetNrIdentifications() != 1:
            raise AssertionError(f"actual IDs={actual}, expected={expected}")
        if mesh.ngmesh.GetIdentificationType(1) != IdentificationType.PERIODIC:
            raise AssertionError("Expected PERIODIC")
        elements = list(mesh.Elements(ng.VOL))
        if len(elements) != 2:
            raise AssertionError(f"actual tets={len(elements)}, expected=2")
        for point in samples:
            ip = next(iter(ng.IntegrationRule(points=[point], weights=[1.0])))
            a = tuple(xyz(mesh.GetTrafo(elements[0])(ip)))
            b = tuple(xyz(mesh.GetTrafo(elements[1])(ip)))
            expected_b = (a[0] + 3.0, a[1], a[2])
            if b != expected_b:
                raise AssertionError(f"q{order} point={point}: actual={b}, expected={expected_b}")


if __name__ == "__main__":
    verify(Path(sys.argv[1]))
