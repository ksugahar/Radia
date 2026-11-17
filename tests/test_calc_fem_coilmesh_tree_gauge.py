"""Tree gauge of the coil A-V solve: spans exactly the coupled gradient kernel."""

import functools
import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, Glue, OCCGeometry, X, Z  # noqa: E402

PANELS = Path(__file__).parents[1] / "src" / "radia" / "panels"
sys.path.insert(0, str(PANELS))
import calc_fem_coilmesh as solver  # noqa: E402


def _coil_vol(tmp_path, split_sibc):
    """Air box with a conducting bar; source/sink are the bar end faces and
    the SIBC trace is the top face, or the top and bottom faces (two
    disconnected components), physically disjoint from the ports."""
    air = Box((0, 0, 0), (1, 1, 1))
    bar = Box((0.3, 0.4, 0.4), (0.7, 0.6, 0.6))
    bar.mat("coil")
    bar.faces.Min(X).name = "source"
    bar.faces.Max(X).name = "sink"
    air = air - bar
    air.mat("air")
    air.faces.Max(Z).name = "sibc"
    if split_sibc:
        air.faces.Min(Z).name = "sibc"
    mesh = OCCGeometry(Glue([air, bar])).GenerateMesh(maxh=0.25)
    path = tmp_path / f"coil_split{int(split_sibc)}.vol"
    mesh.Save(str(path))
    bnds = set(ng.Mesh(str(path)).GetBoundaries())
    assert {"source", "sink", "sibc"} <= bnds, f"boundaries {sorted(bnds)}"
    return str(path)


def _solve(vol, reverse_edges, masks):
    orig = solver._tree_gauged_free_dofs

    def recording(*args, **kwargs):
        free = orig(*args, reverse_edges=reverse_edges, **kwargs)
        masks.append(np.array(list(free), dtype=bool))
        return free

    solver._tree_gauged_free_dofs = recording
    try:
        return solver.solve_fem_coilmesh(
            vol, frequency=7000.0, I_target=1.0, coil_sigma=5.8e7,
            wp_sigma=5.8e7, wp_mu_r=1.0, half_thickness=0.0125,
            impedance_model="sibc")
    finally:
        solver._tree_gauged_free_dofs = orig


@pytest.mark.parametrize("split_sibc", [False, True])
def test_port_driven_observables_do_not_depend_on_the_tree(tmp_path,
                                                          split_sibc):
    """Guards: the tree fixes the coupled kernel (psi = 0 on the ports,
    constant per SIBC trace component) plus the gap circulation, so two
    different trees give the same
    current-normalized L and losses with a direct true residual < 1e-10, also
    when the SIBC trace has two disconnected components."""
    vol = _coil_vol(tmp_path, split_sibc)
    masks = []
    r1 = _solve(vol, False, masks)
    r2 = _solve(vol, True, masks)
    m1, m2 = masks
    assert m1.sum() == m2.sum(), (m1.sum(), m2.sum())
    n_diff = int((m1 != m2).sum())
    assert n_diff > 0, "the two tree masks are identical"
    for r in (r1, r2):
        assert r["I_out_pre_scale_abs"] > 1e-6, r["I_out_pre_scale_abs"]
        assert max(r["linear_true_relative_residuals"]) < 1e-10, r["linear_true_relative_residuals"]
    for key in ("L_total_nH", "P_coil_W", "P_wp_W"):
        rel = abs(r1[key] - r2[key]) / abs(r1[key])
        assert rel < 1e-8, (f"{key}: {r1[key]!r} vs {r2[key]!r} "
                            f"({n_diff} mask entries differ)")


def test_tree_gauge_rejects_missing_ground_boundary(tmp_path):
    """Guards: a port label without vertices fails loudly."""
    mesh = ng.Mesh(_coil_vol(tmp_path, False))
    fesA = ng.HCurl(mesh, order=1, nograds=True)
    with pytest.raises(ValueError, match="no vertices"):
        solver._tree_gauged_free_dofs(mesh, fesA, fesA, ("absent", "sink"),
                                      "sibc")
