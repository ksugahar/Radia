"""Regression checks for the PEEC surface-current Helmholtz split."""

import numpy as np
import pytest


ngsolve = pytest.importorskip("ngsolve")


def test_loop_star_keeps_global_cohomology_modes():
    from netgen.occ import Axes, Dir, OCCGeometry, Pnt, WorkPlane

    from radia.ngsbem_interface import LoopStarTransform
    from radia.ngsbem_peec import create_plate_mesh

    with ngsolve.TaskManager():
        plate = create_plate_mesh(0.1, 0.05, 0.025, label="conductor")
        workplane = WorkPlane(
            Axes(p=Pnt(0, 0, 0), n=Dir(0, 0, 1), h=Dir(1, 0, 0))
        )
        annulus = workplane.Circle(0.05).Face() - workplane.Circle(0.02).Face()
        annulus.name = "conductor"
        ring = ngsolve.Mesh(OCCGeometry(annulus).GenerateMesh(maxh=0.015))

    plate_split = LoopStarTransform(plate)
    ring_split = LoopStarTransform(ring)

    assert plate_split.n_local_loop == len(list(plate.vertices)) - 1
    assert plate_split.n_star == len(list(plate.Elements(ngsolve.BND)))
    assert plate_split.n_harmonic == 0
    assert ring_split.n_harmonic == 1

    for split in (plate_split, ring_split):
        assert split.n_loop + split.n_star == split.n_edge
        assert np.linalg.matrix_rank(split.T) == split.n_edge
        assert np.max(np.abs(split.T_star.T @ split.T_loop)) < 1e-12
