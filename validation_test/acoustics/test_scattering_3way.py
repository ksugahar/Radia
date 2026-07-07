"""3-way validation of radia.acoustics analytic sphere scattering.

  (1) Python analytic == MATLAB analytic  -- machine precision (~1e-14).
      Uses the committed golden matlab_scattering_golden.json (dumped once from
      the matlab-acoustic-fembem reference), so MATLAB is not needed at run time.
  (2) analytic == ngsolve.bem numerical   -- sound-soft sphere, Brakhage-Werner
      combined-field BEM (NGSolve tutorial 11.3); discretization-limited (~2e-5).

Run as a script to (re)write scattering_3way_results.json (data-persistence).
"""
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

import radia.acoustics as ac

HERE = Path(__file__).resolve().parent
GOLDEN = HERE / "matlab_scattering_golden.json"
RESULTS = HERE / "scattering_3way_results.json"


def _cx(a):
    a = np.array(a, float)
    return a[:, 0] + 1j * a[:, 1]


def _load_golden():
    return json.loads(GOLDEN.read_text(encoding="utf-8"))


def _python_vs_matlab(d):
    R, k = d["R"], d["k"]
    P = np.array(d["points"], float)
    Pall = np.vstack([P, np.array(d["points_in"], float)])
    cases = {
        "rigid": (_cx(d["rigid_total"]),
                  ac.rigid_sphere_scattering(k, R, P)["total"]),
        "soft": (_cx(d["soft_total"]),
                 ac.soft_sphere_scattering(k, R, P)["total"]),
        "fluid": (_cx(d["fluid_total"]),
                  ac.fluid_sphere_scattering(k, R, Pall,
                      interior_wavenumber=d["fluid_k1"],
                      density_ratio=d["fluid_rho"])["total"]),
        "elastic": (_cx(d["elastic_total"])[:len(P)],
                    ac.elastic_sphere_scattering(k, R, P,
                        longitudinal_speed=d["elastic_cL"],
                        shear_speed=d["elastic_cT"],
                        density_ratio=d["elastic_rho"])["total"]),
    }
    return {name: float(np.max(np.abs(py - ml)) / max(np.max(np.abs(ml)), 1e-30))
            for name, (ml, py) in cases.items()}


def _analytic_vs_ngsbem(k, R):
    """Sound-soft Brakhage-Werner BEM (NGSolve tutorial 11.3) vs the analytic series.

    Returns (relerr, meta).  The tutorial's potential i*k*SL - DL equals the
    NEGATIVE scattered field (its total = uin - uscat), so both signs are compared.
    """
    from netgen.occ import WorkPlane, Axes, Sphere, Fuse, Compound, X, Y
    from ngsolve import (SurfaceL2, Compress, GridFunction, BilinearForm,
                         LinearForm, TaskManager, ds, exp, z, solvers)
    from ngsolve.bem import HelmholtzCF, HelmholtzSL, HelmholtzDL

    order, maxh = 3, 0.3
    screen = WorkPlane(Axes((0, 0, 0), Y, X)).RectangleC(8, 8).Face()
    sphere = Sphere((0, 0, 0), R)
    screen = screen - sphere
    sp = Fuse(sphere.faces)
    screen.faces.name = "screen"
    sp.faces.name = "sphere"
    mesh = Compound([screen, sp]).GenerateMesh(maxh=maxh).Curve(order)

    fes = Compress(SurfaceL2(mesh, order=order, complex=True,
                             definedon=mesh.Boundaries("sphere")))
    u, v = fes.TnT()
    with TaskManager():
        C = HelmholtzCF(u * ds("sphere"), k) * v * ds
        Id = BilinearForm(u * v * ds).Assemble()
        lhs = 0.5 * Id.mat + C.mat
        rhs = LinearForm(-exp(1j * k * z) * v * ds).Assemble()
        gfu = GridFunction(fes)
        pre = BilinearForm(u * v * ds, diagonal=True).Assemble().mat.Inverse()
        gfu.vec[:] = solvers.GMRes(A=lhs, b=rhs.vec, pre=pre, maxsteps=400, tol=1e-9)
        pot = 1j * k * HelmholtzSL(u * ds("sphere"), k) + (-1) * HelmholtzDL(u * ds("sphere"), k)
        pot_cf = pot(gfu, mesh.Boundaries("screen"))

    probes = np.array([[r * np.sin(np.deg2rad(t)), 0.0, r * np.cos(np.deg2rad(t))]
                       for r in (1.5, 2.0, 3.0) for t in (30, 60, 90, 120, 150)])
    bem = np.array([complex(pot_cf(mesh(px, py, pz))) for px, py, pz in probes])
    ana = ac.soft_sphere_scattering(k, R, probes)["scattered"]
    den = max(np.max(np.abs(ana)), 1e-30)
    relerr = min(float(np.max(np.abs(bem - ana)) / den),
                 float(np.max(np.abs(-bem - ana)) / den))
    return relerr, {"formulation": "brakhage_werner_combined_field",
                    "order": order, "maxh": maxh, "ndof_sphere": int(fes.ndof),
                    "n_probes": int(len(probes))}


def test_python_matches_matlab_golden():
    errs = _python_vs_matlab(_load_golden())
    for name, err in errs.items():
        assert err < 1e-12, f"{name}: python vs matlab = {err:.3e}"


def test_analytic_matches_ngsbem_numerical():
    pytest.importorskip("ngsolve.bem")
    d = _load_golden()
    relerr, _ = _analytic_vs_ngsbem(d["k"], d["R"])
    assert relerr < 5e-3, f"analytic vs ngsolve.bem = {relerr:.3e}"


def _write_results():
    import radia
    import ngsolve
    d = _load_golden()
    two_way = _python_vs_matlab(d)
    try:
        relerr, meta = _analytic_vs_ngsbem(d["k"], d["R"])
        three_way = {"soft_analytic_vs_ngsbem_relerr": relerr, **meta}
    except Exception as exc:  # ngsolve/netgen missing or solve failed
        three_way = {"skipped": repr(exc)}
    out = {
        "schema": "radia.acoustics.scattering-3way.v1",
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "problem": {"radius": d["R"], "wavenumber": d["k"], "convention": "e^{+ikr}",
                    "fluid_k1": d["fluid_k1"], "fluid_rho": d["fluid_rho"],
                    "elastic_cL": d["elastic_cL"], "elastic_cT": d["elastic_cT"],
                    "elastic_rho": d["elastic_rho"]},
        "two_way_python_vs_matlab_relerr": two_way,
        "three_way_analytic_vs_ngsbem": three_way,
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
