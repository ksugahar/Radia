"""Guards: metadata-only charge staging and unchanged CPU quadrature dispatch.

Assumes: stub FE preparation is not numerical or GPU acceptance; the real
NGSolve and native construction lanes remain separately required.
"""
import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_inputs():
    source = ast.parse((ROOT / "src/radia/vim/_vim.py").read_text(encoding="utf-8"))
    names = {"_flat_tet_charge_inputs", "_prepare_gpu_charge_inputs"}
    selected = ast.Module(body=[n for n in source.body if isinstance(n, ast.FunctionDef)
                               and n.name in names], type_ignores=[])
    calls = []
    basis = dict(B=SimpleNamespace(shape=(6, 12)), M_mass=None,
                 M_mass_ngsolve=object(), host=[0, 1], kind=[0, 1],
                 expo=[[0, 0, 0], [1, 0, 0]], cell_verts=object(),
                 face_verts=object(), n_el=1)

    def prepare(fes, quad, **options):
        calls.append(("basis", quad, options))
        return basis

    def rule(kind, quad):
        calls.append((kind, quad))
        return np.array([[quad, 0., 0.]]), np.array([quad + .5])

    ns = dict(np=np, ng=SimpleNamespace(BND="boundary", VOL="volume"),
              _volume_vertex_counts=lambda mesh: mesh.types,
              validate_hdiv_configuration=lambda *args: calls.append(("validate", args)),
              _charge_basis=prepare,
              _outer_tet=lambda q: rule("tet", q),
              _outer_tri=lambda q: rule("tri", q),
              _i32_buffer=lambda a: np.asarray(a, dtype=np.int32).reshape(-1),
              _f64_buffer=lambda a: np.asarray(a, dtype=np.float64).reshape(-1))
    exec(compile(selected, "device-input-functions", "exec"), ns)
    return ns, calls, basis


def space(order=1, curve=1, types=None, boundary=4, dofs=12, dim=3, elements=1):
    mesh = SimpleNamespace(dim=dim, types={4} if types is None else types,
                           GetCurveOrder=lambda: curve,
                           GetNE=lambda region: boundary if region == "boundary" else elements)
    return SimpleNamespace(mesh=mesh, globalorder=order, ndof=dofs)


@pytest.mark.parametrize("order,outer,inner", [(1, 3, 2), (2, 6, 4)])
def test_device_inputs_keep_production_near_inner_far_rules_without_gram(order, outer, inner):
    ns, calls, basis = load_inputs()
    B, payload = ns["_prepare_gpu_charge_inputs"](space(order), backend="gpu")
    assert B is basis["B"]
    assert calls == [("validate", (3, {4}, order, 1)),
                     ("basis", outer, dict(materialize_mass=False, internal_interfaces=False)),
                     ("tet", outer), ("tri", outer), ("tet", 3), ("tri", 3),
                     ("tet", inner), ("tri", inner)]
    assert payload["eps"] == 1e-10
    assert payload["leaf"] == 64 and payload["eta"] == 2.
    assert payload["ho_far_factor"] == 2.
    assert "build" not in payload
    np.testing.assert_array_equal(payload["ref_tet_pts_in"], [inner, 0., 0.])
    np.testing.assert_array_equal(payload["charge_expo"], [0, 0, 0, 1, 0, 0])
    assert payload["charge_expo"].dtype == np.int32
    # No _rp/native constructor exists in ns: reaching either would fail.


@pytest.mark.parametrize("options,error", [
    ({"types": {8}}, NotImplementedError), ({"types": {6}}, NotImplementedError),
    ({"types": {4, 8}}, NotImplementedError), ({"curve": 2}, NotImplementedError),
    ({"order": 0}, NotImplementedError), ({"order": 3}, NotImplementedError),
    ({"dim": 2}, NotImplementedError), ({"boundary": 0}, ValueError),
    ({"elements": 0}, ValueError), ({"elements": 4097}, ValueError),
    ({"dofs": 0}, ValueError), ({"dofs": 16385}, ValueError)])
def test_unsupported_input_rejected_before_fe_preparation(options, error):
    ns, calls, _ = load_inputs()
    with pytest.raises(error):
        ns["_prepare_gpu_charge_inputs"](space(**options), backend="gpu")
    assert calls == []


def test_backend_rejected_before_touching_space():
    ns, calls, _ = load_inputs()
    for backend in (None, "cpu", "auto"):
        with pytest.raises(ValueError, match="explicit backend"):
            ns["_prepare_gpu_charge_inputs"](None, backend=backend)
    assert calls == []


@pytest.mark.parametrize("shape", [(0, 12), (4609, 12), (6, 11)])
def test_charge_dimensions_fail_before_native_construction(shape):
    ns, _, basis = load_inputs()
    basis["B"] = SimpleNamespace(shape=shape)
    with pytest.raises(ValueError, match="map dimensions"):
        ns["_prepare_gpu_charge_inputs"](space(), backend="gpu")


def test_cpu_input_options_preserved_for_explicit_rules_and_images():
    ns, calls, basis = load_inputs()
    cb, payload = ns["_flat_tet_charge_inputs"](
        space(2), 7, eps=1e-9, leafsize=32, eta=3., far_quad=4,
        ho_far_factor=float("inf"), inner_quad=7, image_masks=[1, 2], image_signs=[1., -1.],
        materialize_mass=True, internal_interfaces=True)
    assert cb is basis
    assert calls == [("basis", 7, dict(materialize_mass=True, internal_interfaces=True)),
                     ("tet", 7), ("tri", 7)]
    assert "ref_tet_pts_in" not in payload and "ref_tet_pts_lo" not in payload
    np.testing.assert_array_equal(payload["image_masks"], [1, 2])
    np.testing.assert_array_equal(payload["image_signs"], [1., -1.])
    assert (payload["eps"], payload["leaf"], payload["eta"]) == (1e-9, 32, 3.)
