# -*- coding: utf-8 -*-
"""High-order CURVED HEX mesh in NGSolve -- the cubit-mesh-export capability.

Coreform Cubit meshes a sphere with HEXES ('volume scheme sphere') and exports a
high-order Netgen .vol (order 1/2/3) via the plugin command

    block 1 add hex all
    export netgen "hexsph_oN.vol" order N overwrite

NGSolve loads the .vol directly and integrates the volume; as the element order rises
the curved hex boundary captures the sphere and the volume converges to 4/3 pi r^3:

    order 1 (straight hex) :  ~ -23 %    (a coarse 56-hex faceted sphere)
    order 2 (curved hex)   :  ~ -0.2 %
    order 3 (curved hex)   :  ~ +0.1 %

That is high-order *hexahedral* FEM geometry inside NGSolve -- which the Netgen mesher
alone does not provide (it curves tetrahedra).  See packages/cubit-mesh-export.

GOTCHA (important): a high-order .vol already carries its curved mid-side nodes.  Load
it with plain ``Mesh(path)`` and DO NOT call ``mesh.Curve()`` -- ``mesh.Curve(p)``
re-curves from the underlying CAD geometry, which a loaded .vol does NOT have, so it
RESETS every element to straight-sided (the volume jumps back to the -23 % linear value).
(``mesh.Curve`` is the right call only when the mesh was built from an in-memory geometry,
e.g. a netgen.occ / SplineGeometry CAD object.)

Pure NGSolve -- no Cubit needed to RUN this (the three .vol files are committed beside
it).  Regenerate them with Cubit if you change the sphere.  Headless:

    python hex_sphere_curved_ngsolve.py
"""
import math
import os
from ngsolve import Mesh, Integrate, CoefficientFunction, VOL

R = 0.05
V_EXACT = 4.0 / 3.0 * math.pi * R ** 3
HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    print(f"Curved-hex sphere in NGSolve  (R={R} m, V_exact={V_EXACT:.6e} m^3)")
    worst_hi = 0.0
    for order in (1, 2, 3):
        mesh = Mesh(os.path.join(HERE, f"hexsph_o{order}.vol"))   # load AS-IS (no .Curve)
        types = {}
        for el in mesh.Elements(VOL):
            k = str(el.type).replace("ET.", "")
            types[k] = types.get(k, 0) + 1
        V = Integrate(CoefficientFunction(1), mesh, VOL)
        err = (V - V_EXACT) / V_EXACT * 100.0
        print(f"  order {order}  {types}  V={V:.6e}  err={err:+.3f}%")
        if order >= 2:
            worst_hi = max(worst_hi, abs(err))
    assert worst_hi < 1.0, f"high-order hex volume err {worst_hi:.3f}% should be < 1%"
    print("PASS: curved high-order hex converges to the analytic sphere (order>=2 < 1%)")


if __name__ == "__main__":
    main()
