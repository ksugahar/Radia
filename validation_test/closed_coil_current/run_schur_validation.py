"""Validate RT0 pressure Schur solves at larger sizes and two length scales.

This explicitly uses Netgen's OCC mesher for a synthetic closed-loop test.
It is a convergence/conservation check, not a timing comparison.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import ngsolve as ng
from netgen.occ import Box, Pnt, Glue, OCCGeometry
import radia.meshed_current as current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maxh", type=float, default=0.18)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ng.SetNumThreads(args.threads)
    cases = []
    for scale in (1.0, 1e-3):
        corners = [((1, -1, 0), (2, 1, 1)), ((-2, -1, 0), (-1, 1, 1)),
                   ((-2, 1, 0), (2, 2, 1)), ((-2, -2, 0), (2, -1, 1))]
        parts = [Box(Pnt(*(scale * v for v in low)),
                     Pnt(*(scale * v for v in high))) for low, high in corners]
        for i, part in enumerate(parts):
            part.mat("drive" if i == 0 else "return")
        mesh = ng.Mesh(OCCGeometry(Glue(parts)).GenerateMesh(maxh=args.maxh * scale))
        drive = mesh.MaterialCF({"drive": ng.CF((0, 1, 0))}, default=ng.CF((0, 0, 0)))
        with ng.TaskManager():
            result = current.solve_closed_coil_current(
                mesh, drive, current_A=2, drive_length_m=2 * scale,
                materials=("drive", "return"))
        stats = dict(result["stats"])
        # Runtime is not accepted as performance evidence on this lane.
        stats.pop("runtime_s")
        assert stats["linear_relative_residual"] <= 1e-10
        assert stats["relative_divergence"] <= 1e-8
        assert abs(stats["section_current_A"] - 2) <= 1e-10
        cases.append(dict(scale_m=scale, mesh_elements=mesh.ne, **stats))
    record = dict(schema="radia.closed-current-schur-validation.v1",
                  platform_class=platform.system(), ngsolve=ng.__version__, threads=args.threads,
                  mesher="Netgen OCC synthetic loop", relative_maxh=args.maxh,
                  source_sha256=hashlib.sha256(Path(current.__file__).read_bytes()).hexdigest(),
                  cases=cases, accepted=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
