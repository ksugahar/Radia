"""Checked-VOL handoff of the prolate HDiv case before Solve.

Guards: run_prolate_case refuses a generated named prolate mesh, or any VOL whose
sha256 is not pinned with its strict label-check report, before ng.Mesh or Solve;
the pinned named fixture still reaches Solve with its labels and element count.
Assumes: Solve is replaced by a recorder, so no numerical solve runs; the maxh 0.6
mesh is the one the validation entry generated before the checked-VOL handoff.
"""
import importlib.util
from pathlib import Path

import pytest

ng = pytest.importorskip("ngsolve")
from netgen.csg import CSGeometry, Ellipsoid, Pnt, Vec  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = Path(__file__).parent / "data" / "prolate_named" / "prolate_named.vol"


def _case():
    spec = importlib.util.spec_from_file_location(
        "prolate_checked_input_case", ROOT / "validation_test/feec/test_hdiv_vim_ellipsoid.py")
    case = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(case)
    return case


def test_generated_or_unpinned_vol_stops_before_solve(tmp_path, monkeypatch):
    case = _case()
    geo = CSGeometry()
    geo.Add(Ellipsoid(Pnt(0, 0, 0), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 2.0)).mat("magnetic_body").bc("exterior"))
    generated = tmp_path / "generated_named.vol"
    with ng.TaskManager():
        geo.GenerateMesh(maxh=0.6).Save(str(generated))
    calls = []
    monkeypatch.setattr(case, "Solve", lambda *args, **kwargs: calls.append("Solve"))
    monkeypatch.setattr(case.ng, "Mesh", lambda *args, **kwargs: calls.append("Mesh"))
    with pytest.raises(ValueError, match="checked VOL handoff required"):
        case.run_prolate_case(generated)
    with pytest.raises(FileNotFoundError):
        case.run_prolate_case(tmp_path / "missing.vol")
    assert calls == [], f"native calls before refusal: actual={calls}; expected=[]"


def test_pinned_fixture_reaches_solve(monkeypatch):
    case = _case()
    seen = []

    def record(mesh, **kwargs):
        labels = (tuple(mesh.GetMaterials()), tuple(sorted(set(mesh.GetBoundaries()))))
        seen.append((labels, mesh.ne, kwargs["mu_r"]))
        return {"demag": None}

    monkeypatch.setattr(case, "Solve", record)
    case.run_prolate_case(FIXTURE)
    expected = [(case.PROLATE_LABELS, 1338, 1000.0)]
    assert seen == expected, f"Solve inputs: actual={seen}; expected={expected}"
