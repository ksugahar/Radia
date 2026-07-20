"""Independent pure-double-layer validation for sound-soft sphere scattering.

This deliberately excludes the single-layer operator and the Brakhage--Werner
combined field.  With the outward normal used by NGSolve, solve

    (1/2 I + K) mu = -u_inc

and evaluate u_scat = D mu at exterior probes.  The chosen kR=2 case is away
from an irregular frequency; combined fields remain the production choice near
interior resonances.
"""
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

import radia.acoustics as ac

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "double_layer_bem_results.json"


def _pure_double_layer_vs_analytic(k=2.0, radius=1.0, order=3, maxh=0.3):
    from netgen.occ import WorkPlane, Axes, Sphere, Fuse, Compound, X, Y
    from ngsolve import (SurfaceL2, Compress, GridFunction, BilinearForm,
                         LinearForm, TaskManager, ds, exp, z, solvers)
    from ngsolve.bem import HelmholtzDL

    screen = WorkPlane(Axes((0, 0, 0), Y, X)).RectangleC(8, 8).Face()
    sphere = Sphere((0, 0, 0), radius)
    screen = screen - sphere
    boundary = Fuse(sphere.faces)
    screen.faces.name = "screen"
    boundary.faces.name = "sphere"
    mesh = Compound([screen, boundary]).GenerateMesh(maxh=maxh).Curve(order)
    fes = Compress(SurfaceL2(mesh, order=order, complex=True,
                             definedon=mesh.Boundaries("sphere")))
    u, v = fes.TnT()

    with TaskManager():
        identity = BilinearForm(u * v * ds).Assemble()
        double_layer = HelmholtzDL(u * ds("sphere"), k) * v * ds
        lhs = 0.5 * identity.mat + double_layer.mat
        incident_trace = LinearForm(exp(1j * k * z) * v * ds).Assemble()
        rhs = incident_trace.vec.CreateVector()
        rhs.data = -incident_trace.vec
        pre = BilinearForm(u * v * ds, diagonal=True).Assemble().mat.Inverse()
        density = GridFunction(fes)
        density.vec[:] = solvers.GMRes(
            A=lhs, b=rhs, pre=pre, maxsteps=800, tol=1e-10,
            printrates=False,
        )
        potential = HelmholtzDL(u * ds("sphere"), k)
        potential_cf = potential(density, mesh.Boundaries("screen"))

    probes = np.array([
        [r * np.sin(np.deg2rad(theta)), 0.0,
         r * np.cos(np.deg2rad(theta))]
        for r in (1.5, 2.0, 3.0)
        for theta in (30, 60, 90, 120, 150)
    ])
    numerical = np.array([
        complex(potential_cf(mesh(x, y, z0))) for x, y, z0 in probes
    ])
    analytic = ac.soft_sphere_scattering(k, radius, probes)["scattered"]
    rel_error = float(
        np.max(np.abs(numerical - analytic)) /
        max(np.max(np.abs(analytic)), 1e-30)
    )
    return rel_error, {
        "formulation": "pure_double_layer",
        "boundary_equation": "(1/2 I + K) mu = -u_inc",
        "field_representation": "u_scat = D mu",
        "wavenumber": k,
        "radius": radius,
        "order": order,
        "maxh": maxh,
        "ndof_sphere": int(fes.ndof),
        "n_probes": len(probes),
    }


def test_pure_double_layer_matches_soft_sphere_analytic():
    pytest.importorskip("ngsolve.bem")
    rel_error, _ = _pure_double_layer_vs_analytic()
    assert rel_error < 5e-3, rel_error


def _write_results():
    import ngsolve
    import radia

    rel_error, metadata = _pure_double_layer_vs_analytic()
    result = {
        "schema": "radia.acoustics.double-layer-bem.v1",
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "validation_host": platform.node(),
        **metadata,
        "relative_error_vs_soft_sphere_analytic": rel_error,
        "tolerance": 5e-3,
        "passed": rel_error < 5e-3,
        "versions": {
            "radia": getattr(radia, "__version__", "?"),
            "ngsolve": ngsolve.__version__,
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
    }
    RESULTS.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    _write_results()
