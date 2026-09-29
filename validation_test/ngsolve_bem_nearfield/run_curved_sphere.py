"""Check curved-sphere Laplace single-layer values and gradients near the surface.

The unit sphere with unit density has potential 1 inside and 1/r outside.
Errors are absolute RMS, normalized by the unit surface potential/field scale;
the interior field is zero, so an interior relative field error is undefined.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
from pathlib import Path

import ngsolve as ng
from ngsolve.bem import LaplaceSL
from netgen.occ import Axes, OCCGeometry, Pnt, Sphere, WorkPlane, X, Z


LIMITS = {"potential_rms": 1e-6, "gradient_rms": 1e-4}


def evaluate(maxh, bonus):
    source = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1)).GenerateMesh(maxh=maxh))
    source.Curve(5)
    fes = ng.SurfaceL2(source, order=0)
    density = ng.GridFunction(fes)
    density.Set(1, definedon=source.Boundaries(".*"))
    op = LaplaceSL(fes.TrialFunction() * ng.ds(bonus_intorder=bonus))
    radius = ng.sqrt(ng.x**2 + ng.y**2 + ng.z**2)
    area = ng.Integrate(1, source, definedon=source.Boundaries(".*"))
    geometry_rms = math.sqrt(ng.Integrate((radius-1)**2, source,
                                        definedon=source.Boundaries(".*"), order=16) / area)
    rows = []
    for side in (-1, 1):
        for distance in (.1, .01, .001, 1e-5):
            face = WorkPlane(Axes((0, 0, 1 + side*distance), Z, X)).RectangleC(.002, .002).Face()
            face.faces.name = "target"
            target = ng.Mesh(OCCGeometry(face).GenerateMesh(maxh=.001))
            region = target.Boundaries("target")
            target_area = ng.Integrate(1, target, definedon=region)
            exact = 1/radius if side == 1 else ng.CF(1)
            exact_grad = -ng.CF((ng.x, ng.y, ng.z))/radius**3 if side == 1 else ng.CF((0, 0, 0))
            row = {"side": "outside" if side == 1 else "inside", "distance": distance}
            for name, potential in (("direct", op(density)), ("local_expansion", op(density, region))):
                delta = ng.grad(potential) - exact_grad
                row[name] = {
                    "potential_rms": math.sqrt(abs(ng.Integrate((potential-exact)**2, target,
                                                              definedon=region) / target_area)),
                    "gradient_rms": math.sqrt(abs(ng.Integrate(ng.InnerProduct(delta, delta), target,
                                                             definedon=region) / target_area)),
                }
            rows.append(row)
    passed = all(math.isfinite(row[route][key]) and row[route][key] < limit
                 for row in rows for route in ("direct", "local_expansion")
                 for key, limit in LIMITS.items())
    return {"maxh": maxh, "bonus_intorder": bonus, "curving_order": 5,
            "source_elements": source.ne, "surface_dofs": fes.ndof,
            "geometry_radial_rms": geometry_rms, "passed": passed, "results": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ng.SetNumThreads(1)
    with ng.TaskManager():
        cases = [evaluate(h, bonus) for h, bonus in ((.35, 12), (.35, 24), (.175, 12), (.175, 24))]
    out = {"schema": "radia.bem-curved-sphere/v1", "ngsolve": ng.__version__,
           "python": platform.python_version(), "host": platform.node(),
           "thresholds": LIMITS, "accepted_case": {"maxh": .175, "bonus_intorder": 24},
           "passed": cases[-1]["passed"], "cases": cases}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": out["passed"], "cases": [
        {"maxh": c["maxh"], "bonus": c["bonus_intorder"], "passed": c["passed"]} for c in cases]}))
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
