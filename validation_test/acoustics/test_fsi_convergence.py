"""FSI convergence to the analytic Faran elastic sphere (radia.acoustics.fsi).

The DtN-coupled elastic-sphere scattered field converges to the analytic Faran
series (elastic_sphere_scattering) under mesh refinement; the NGSolve P2 interior
converges ~O(h^2), much faster than P1 ~O(h).  Run as a script to (re)write
fsi_convergence_results.json.
"""
import json
import platform
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from ngsolve import TaskManager

from radia.acoustics import elastic_sphere_scattering, fsi

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "fsi_convergence_results.json"

K, R = 2.0, 1.0
MAT = {"cL": 2.0, "cT": 1.0, "rho_s": 1.5}
OBS = np.array([[0.0, 0.0, -3.0], [3.0, 0.0, 0.0], [0.0, 0.0, 3.0],
                [0.0, 0.0, 2.0], [2.0, 0.0, 1.0]])
MAXH = [0.4, 0.28, 0.20]


def _rel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def _sweep():
    with TaskManager():
        far = elastic_sphere_scattering(K, R, OBS, longitudinal_speed=MAT["cL"],
                                        shear_speed=MAT["cT"], density_ratio=MAT["rho_s"])["scattered"]
        errs = {1: [], 2: []}
        ndof = {1: [], 2: []}
        for h in MAXH:
            mesh = fsi.sphere_mesh(R, maxh=h)
            for p in (1, 2):
                s = fsi.fsi_dtn_solve(mesh, K, cL=MAT["cL"], cT=MAT["cT"],
                                      rho_s=MAT["rho_s"], order=p, obs=OBS)
                errs[p].append(_rel(s["scattered"], far))
                ndof[p].append(s["ndof_u"])
    return errs, ndof


def test_fsi_converges_to_faran_and_p2_beats_p1():
    errs, _ = _sweep()
    # both interior orders converge (monotone decreasing) to the analytic Faran field
    assert all(np.diff(errs[1]) < 0), errs[1]
    assert all(np.diff(errs[2]) < 0), errs[2]
    # P2 is more accurate than P1 at every refinement level
    assert all(e2 < e1 for e1, e2 in zip(errs[1], errs[2])), (errs[1], errs[2])
    # P2 reaches engineering accuracy on the finest mesh here
    assert errs[2][-1] < 0.06, errs[2]


def _write_results():
    import ngsolve

    import radia
    errs, ndof = _sweep()
    out = {
        "schema": "radia.acoustics.fsi-convergence.v1",
        "generated_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "problem": {"radius": R, "wavenumber": K, "exterior": "spherical_dtn", **MAT},
        "maxh": MAXH,
        "rel_error_vs_faran": {"P1": errs[1], "P2": errs[2]},
        "ndof_interior": {"P1": ndof[1], "P2": ndof[2]},
        "versions": {"radia": getattr(radia, "__version__", "?"),
                     "ngsolve": ngsolve.__version__,
                     "python": platform.python_version(),
                     "platform": platform.platform()},
    }
    RESULTS.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    _write_results()
    sys.exit(0)
