"""Permanent-magnet (fixed-M) regions mixed with LINEAR soft iron in radia.vim.hdiv_demag_solve
(pm_M = {material: [Mx,My,Mz]}).  A PM region is the M = M_pm (no H-response, chi = 0) special case of
the form-1 projected statement: it contributes a source b_pm = INT M_pm.v dx and its field reaches the
iron through the full-m demag, while its own M stays pinned to M_pm.

Scope + the verified boundary (productionization, docs/hdiv_vim/PRODUCTIONIZATION.md "per-region /
mixed"):
  - PM-only -> M = M_pm exactly (the RT0 projection of a constant vector is exact).
  - PM + LINEAR iron SEPARATED BY AN AIR GAP (the usual magnetic-circuit / PM-motor case) matches
    rad.Solve (yano) to ~1e-3 on the external field and CONVERGES under refinement.
  - PM directly TOUCHING soft iron is NOT supported and FAILS LOUD: the conforming HDiv (RT0) field
    cannot represent the magnetization discontinuity across a PM-iron boundary, so a shared-facet
    interface DIVERGES under refinement (~20% external-field error).  The touching-interface formulation
    is the next step.

This is a transitional gate: it compares against yano-type while it still ships (M5 retires it)."""
import warnings

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import radia as rad  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt, Glue  # noqa: E402

from radia.vim import hdiv_demag_solve  # noqa: E402
from radia.netgen_mesh_import import extract_elements, create_radia_tetrahedron  # noqa: E402

M0 = 9.0e5          # PM magnetization (A/m), +z (~Br 1.13 T)
MU = 1000.0         # soft-iron mu_r
L = 0.02            # block edge (m)
ZERO = ng.CoefficientFunction((0.0, 0.0, 0.0))


def _pm_iron_mesh(gap, maxh):
    """PM block (left) + iron block (right), separated by `gap` in x (gap > 0 -> no shared facets)."""
    pm = Box(Pnt(-L - gap / 2, -L / 2, -L / 2), Pnt(-gap / 2, L / 2, L / 2)).mat("pm")
    iron = Box(Pnt(gap / 2, -L / 2, -L / 2), Pnt(L + gap / 2, L / 2, L / 2)).mat("iron")
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(Glue([pm, iron])).GenerateMesh(maxh=maxh))


def _ext_points():
    return np.array([[3 * L * np.cos(a), 3 * L * np.sin(a), 0.5 * L]
                     for a in np.linspace(0, 2 * np.pi, 24, endpoint=False)])


def _yano_external_B(mesh, ext):
    """rad.Solve (yano/MMM) on the SAME tets: pm -> fixed M0 z, iron -> soft MatLin(MU); external B."""
    rad.UtiDelAll()
    els, _ = extract_elements(mesh)
    matlin = rad.MatLin(MU)
    objs = []
    for e in els:
        reg = mesh[ng.ElementId(ng.VOL, e["element_index"])].mat
        if reg == "pm":
            objs.append(create_radia_tetrahedron(e["vertices"], [0.0, 0.0, M0]))
        else:
            oid = create_radia_tetrahedron(e["vertices"], [0.0, 0.0, 0.0])
            rad.MatApl(oid, matlin)
            objs.append(oid)
    cont = rad.ObjCnt(objs)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)   # yano backend is deprecated; expected here
        rad.Solve(cont, 1e-6, 3000, 0)                        # LU
    return np.array([rad.Fld(cont, "b", p.tolist()) for p in ext])


def _hdiv_external_B(mesh, M_el, ext):
    """External B of the HDiv-solved magnetization state (per-element M -> Radia tets -> rad.Fld)."""
    rad.UtiDelAll()
    els, _ = extract_elements(mesh)
    objs = [create_radia_tetrahedron(e["vertices"], M_el[e["element_index"]].tolist()) for e in els]
    cont = rad.ObjCnt(objs)
    return np.array([rad.Fld(cont, "b", p.tolist()) for p in ext])


def test_pm_only_returns_M_pm():
    """A pure-PM body (no iron) returns M = M_pm exactly (constant-vector RT0 projection is exact)."""
    b = Box(Pnt(-L / 2, -L / 2, -L / 2), Pnt(L / 2, L / 2, L / 2)).mat("pm")
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(b).GenerateMesh(maxh=L / 3))
        r = hdiv_demag_solve(mesh, pm_M={"pm": [0.0, 0.0, M0]}, H_ext=ZERO)
    assert np.allclose(r["M_avg"], [0.0, 0.0, M0], rtol=1e-6, atol=1.0), r["M_avg"]


def test_pm_iron_airgap_parity_vs_yano():
    """PM + LINEAR iron separated by an air gap: external B matches rad.Solve (yano) to < 1% on the SAME
    tet mesh (the method-isolating parity; measured ~3e-4)."""
    mesh = _pm_iron_mesh(gap=0.5 * L, maxh=L / 4)
    ext = _ext_points()
    By = _yano_external_B(mesh, ext)
    with ng.TaskManager():
        r = hdiv_demag_solve(mesh, mu_r={"iron": MU}, pm_M={"pm": [0.0, 0.0, M0]}, H_ext=ZERO)
    Bh = _hdiv_external_B(mesh, r["M"], ext)
    rel = np.linalg.norm(Bh - By) / np.linalg.norm(By)
    assert rel < 1e-2, f"PM+iron(air gap) HDiv vs yano external-B rel {rel:.2e}"


def test_pm_iron_touching_fail_loud():
    """PM directly TOUCHING soft iron (shared RT0 facets) must FAIL LOUD -- the conforming field cannot
    represent the PM-iron magnetization discontinuity (diverges under refinement)."""
    mesh = _pm_iron_mesh(gap=0.0, maxh=L / 4)
    with pytest.raises(NotImplementedError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r={"iron": MU}, pm_M={"pm": [0.0, 0.0, M0]}, H_ext=ZERO)


def test_pm_mixing_classification_fail_loud():
    """A region that is both iron and PM, an unclassified region, or a missing PM region must raise."""
    mesh = _pm_iron_mesh(gap=0.5 * L, maxh=L / 3)
    with pytest.raises(ValueError):                            # 'pm' declared both iron and PM
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r={"pm": MU, "iron": MU}, pm_M={"pm": [0, 0, M0]}, H_ext=ZERO)
    with pytest.raises(ValueError):                            # 'pm' unclassified (iron dict omits it), 'nope' missing
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r={"iron": MU}, pm_M={"nope": [0, 0, M0]}, H_ext=ZERO)
