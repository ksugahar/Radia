"""Stage-2/3 gates at the research configuration (golden bands + record JSON).

Sector-pole surrogate (annular sector, maxh 0.02, ~194 tets): arc-orbit
radial dB_z/dr objective, two mean-B_z arc constraints held in 0.5 % bands,
50 % iron volume budget, Helmholtz filter.  Golden bands from the measured
2026-07-28 values: J strictly monotone with a total gain of +16.1 %
(band [8 %, 30 %]); violations bounded by 1.25 band (peak 1.06 band);
matched-0/1 exact-void ersatz band +0.69 % (band |.| < 3 %); the gray-design
gap -96 % (detection lock < -50 %).  The aggregated record is written to
``results_design_loop_lane.json`` next to this file (committed).
"""
import json
import platform
from datetime import datetime, timezone
from math import cos, pi, sin
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia")

from netgen.occ import Cylinder, HalfSpace, OCCGeometry, Pnt, Vec, Z  # noqa: E402
from ngsolve import HDiv, InnerProduct, Mesh, SetNumThreads, TaskManager  # noqa: E402

import radia  # noqa: E402
from radia.isochronous_topopt import (  # noqa: E402
    MU0, DensityAdjointVIM, HelmholtzFilter, density_to_s,
    field_functional_load, gradient_pair_points, optimize_density,
    orbit_arc_points, uniform_field_load, verify_design_iron_only,
)

RESULTS = Path(__file__).with_name("results_design_loop_lane.json")


def make_sector_mesh(R1=0.05, R2=0.15, angle_deg=60.0, z0=0.02, thick=0.03,
                     maxh=0.02):
    ang = angle_deg * pi / 180.0
    ring = (Cylinder(Pnt(0, 0, z0), Z, r=R2, h=thick)
            - Cylinder(Pnt(0, 0, z0 - 0.01), Z, r=R1, h=thick + 0.02))
    hs1 = HalfSpace(Pnt(0, 0, 0), Vec(0, -1, 0))
    hs2 = HalfSpace(Pnt(0, 0, 0), Vec(-sin(ang), cos(ang), 0))
    return Mesh(OCCGeometry(ring * hs1 * hs2).GenerateMesh(maxh=maxh))


@pytest.fixture(scope="module")
def study():
    SetNumThreads(4)
    chi_iron = 1000.0
    span = (pi / 12, pi / 2 - pi / 12)
    obj_pts, obj_radial = orbit_arc_points(0.115, 0.0, 7, span=span)
    pair_pts, pair_wts = gradient_pair_points(
        obj_pts, np.full(len(obj_pts), 1.0 / len(obj_pts)), delta=0.01,
        direction=obj_radial)
    con_radii = (0.08, 0.10)

    def state_builder(fes):
        return uniform_field_load(fes, (0.0, 0.0, 1.0e5))

    def objective_builder(fes):
        return field_functional_load(fes, pair_pts, pair_wts, axis=2,
                                     scale=MU0, bonus_intorder=10)

    def constraint_builder(radius):
        def build(fes):
            cpts, _ = orbit_arc_points(radius, 0.0, 7, span=span)
            return field_functional_load(
                fes, cpts, np.full(len(cpts), 1.0 / len(cpts)), axis=2,
                scale=MU0, bonus_intorder=10)
        return build

    con_builders = [constraint_builder(r) for r in con_radii]
    with TaskManager():
        mesh = make_sector_mesh()
        fes = HDiv(mesh, order=1)
        prob = DensityAdjointVIM(fes, eps=1e-7)
        f_state = state_builder(fes)
        f_obj = objective_builder(fes)
        cons = [b(fes) for b in con_builders]
        filt = HelmholtzFilter(mesh, radius=0.012)
        lin0 = prob.linearize(
            density_to_s(filt.apply(np.full(prob.n_el, 0.5)), chi_iron),
            f_state, [f_obj] + cons)
        targets = [float(v) for v in lin0.values[1:]]
        result = optimize_density(prob, f_state, f_obj, cons, targets,
                                  chi_iron=chi_iron, volume_fraction=0.5,
                                  density_filter=filt, move_limit=0.1,
                                  max_iterations=30)
        verification = verify_design_iron_only(
            prob, result.density, state_builder,
            [objective_builder] + con_builders, chi_iron=chi_iron,
            density_filter=filt, gram_kwargs=dict(eps=1e-7))
    return SimpleNamespace(mesh=mesh, prob=prob, lin0=lin0, targets=targets,
                           result=result, verification=verification,
                           chi_iron=chi_iron)


def test_monotone_constrained_ascent(study):
    hist = study.result.history
    assert len(hist) >= 10
    objectives = [h["objective"] for h in hist]
    assert all(b >= a * (1.0 - 1e-6) for a, b in zip(objectives,
                                                     objectives[1:]))
    gain = objectives[-1] / objectives[0] - 1.0
    assert 0.08 < gain < 0.30, gain          # measured +16.1 %
    for h in hist:
        ratio = max(np.array(h["violation"]) / np.array(h["band"]))
        assert ratio <= 1.25 + 1e-9, (h["iteration"], ratio)


def test_exact_void_ersatz_band(study):
    band = float(study.verification.bands[0])   # objective functional
    assert abs(band) < 0.03, band               # measured +0.69 %
    assert study.verification.iron_mesh.ne == int(study.verification.keep.sum())


def test_gray_design_gap_detected(study):
    J_cont = study.result.history[-1]["objective"]
    J_bin = float(study.verification.values_embedded[0])
    gap = (J_cont - J_bin) / abs(J_bin)
    assert gap < -0.5, gap                      # measured -96 %


def test_write_record_json(study):
    hist = study.result.history
    record = dict(
        schema="radia.isochronous-topopt-lane/v1",
        generated_at_utc=datetime.now(timezone.utc)
        .strftime("%Y-%m-%dT%H:%M:%SZ"),
        hostname=platform.node(), radia_version=radia.__version__,
        python_version=platform.python_version(),
        ne=int(study.mesh.ne), ndof=int(study.prob.fes.ndof),
        targets=study.targets,
        J_first=hist[0]["objective"], J_last=hist[-1]["objective"],
        gain=hist[-1]["objective"] / hist[0]["objective"] - 1.0,
        viol_over_band_max=float(max(
            max(np.array(h["violation"]) / np.array(h["band"]))
            for h in hist)),
        t_iter_median_s=float(np.median([h["t_iter_s"] for h in hist])),
        solves=study.result.solves,
        n_iron=int(study.verification.keep.sum()),
        values_embedded=study.verification.values_embedded.tolist(),
        values_iron_only=study.verification.values_iron_only.tolist(),
        bands=study.verification.bands.tolist(),
        history=list(hist),
    )
    RESULTS.write_text(json.dumps(record, indent=1), encoding="utf-8")
    assert RESULTS.exists()
