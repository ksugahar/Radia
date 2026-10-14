"""ChargeGram fail-loud gate: a 3D mesh with no boundary faces must raise.

A bounded 3D body always has a mesh skin; a ``.vol`` whose
``surfaceelements`` section is empty is a broken export.  Without BND
faces the surface charge ``sigma = M.n`` cannot be represented, so the
demag operator silently acts volume-charge-only and uniform magnetization
sees ``N = 0``.  ``build_charge_gram`` now raises instead; this test locks
the gate and its positive control.
"""

import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

from netgen.occ import Box, OCCGeometry, Pnt

from radia import vim

_SURFACE_TOKENS = ("surfaceelements", "surfaceelementsgi",
                   "surfaceelementsuv")


def _strip_surface_elements(vol_path, out_path):
    lines = vol_path.read_text(encoding="utf-8").splitlines(keepends=True)
    out, i, stripped = [], 0, False
    while i < len(lines):
        if lines[i].strip() in _SURFACE_TOKENS:
            count = int(lines[i + 1])
            out.append(lines[i])
            out.append("0\n")
            i += 2 + count
            stripped = True
            continue
        out.append(lines[i])
        i += 1
    assert stripped, "saved .vol had no surfaceelements section to strip"
    out_path.write_text("".join(out), encoding="utf-8")


def test_chargegram_raises_on_3d_mesh_without_boundary_faces(tmp_path):
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(
            Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.6))
        vol = tmp_path / "cube.vol"
        mesh.ngmesh.Save(str(vol))

        bald = tmp_path / "cube_nobnd.vol"
        _strip_surface_elements(vol, bald)

        broken = ng.Mesh(str(bald))
        assert broken.GetNE(ng.BND) == 0, "strip failed to remove BND faces"
        with pytest.raises(ValueError, match="ZERO boundary"):
            vim.DemagOperator(ng.HDiv(broken, order=1), eps=1e-4)

        intact = ng.Mesh(str(vol))
        assert intact.GetNE(ng.BND) > 0
        operator = vim.DemagOperator(ng.HDiv(intact, order=1), eps=1e-4)
        factor = operator.DemagFactor(
            ng.CoefficientFunction((0.0, 0.0, 1.0)))
        # cube demag factor ~ 1/3; the point is only that surface charges
        # are alive, so lock a generous band
        assert 0.15 < factor < 0.6
