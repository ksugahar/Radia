"""Fast ABI and functional contract for the pinned NGSolve/Netgen runtime."""

from __future__ import annotations

import importlib.metadata
import importlib.util
import sys
import tomllib
from pathlib import Path

import numpy as np
import pytest
from netgen.meshing import (
    EdgeDescriptor,
    Element1D,
    Element2D,
    FaceDescriptor,
    Mesh as NetgenMesh,
    MeshPoint,
    Pnt,
)

ng = pytest.importorskip("ngsolve")
MakeStructured3DMesh = pytest.importorskip(
    "ngsolve.meshes"
).MakeStructured3DMesh


ROOT = Path(__file__).resolve().parents[1]
PINNED_VERSION = "6.2.2606"


def _dependencies(path: Path) -> set[str]:
    project = tomllib.loads(path.read_text(encoding="utf-8"))["project"]
    return set(project["dependencies"])


def test_ngsolve_netgen_runtime_matches_both_exact_package_pins():
    assert ng.__version__ == PINNED_VERSION
    assert importlib.metadata.version("netgen-mesher") == PINNED_VERSION

    radia_dependencies = _dependencies(ROOT / "pyproject.toml")
    cubit_dependencies = _dependencies(
        ROOT / "packages" / "cubit-mesh-export" / "pyproject.toml"
    )
    expected = {
        f"ngsolve=={PINNED_VERSION}",
        f"netgen-mesher=={PINNED_VERSION}",
    }
    assert expected <= radia_dependencies
    assert expected <= cubit_dependencies


@pytest.mark.parametrize(
    ("family", "mesh_kwargs", "expected_types"),
    [
        ("TET", {"hexes": False, "prism": False}, {"ET.TET"}),
        ("HEX", {"hexes": True, "prism": False}, {"ET.HEX"}),
        ("WEDGE", {"hexes": False, "prism": True}, {"ET.PRISM"}),
    ],
)
def test_supported_hdiv_order2_spaces_assemble_and_reproduce_constants(
    family, mesh_kwargs, expected_types
):
    mesh = MakeStructured3DMesh(nx=1, ny=1, nz=1, **mesh_kwargs)
    assert {str(element.type) for element in mesh.Elements(ng.VOL)} == expected_types

    with ng.TaskManager():
        space = ng.HDiv(mesh, order=2)
        trial, test = space.TnT()
        mass = ng.BilinearForm(ng.InnerProduct(trial, test) * ng.dx).Assemble()
        assert mass.mat.height == space.ndof > 0, family

        target = ng.CoefficientFunction((1.0, 2.0, 3.0))
        field = ng.GridFunction(space)
        field.Set(target)
        error2 = float(
            ng.Integrate(ng.InnerProduct(field - target, field - target), mesh)
        )
    assert np.isfinite(error2)
    assert error2 < 1e-24, (family, error2)


def test_hdiv_pyramid_remains_an_explicit_upstream_tripwire():
    module_path = ROOT / "tools" / "probe_hdiv_pyramid.py"
    spec = importlib.util.spec_from_file_location("probe_hdiv_pyramid", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    result = module.probe()
    assert result["ngsolve_version"] == PINNED_VERSION
    assert result["detail"]["mesh_valid_h1"] is True
    assert result["verdict"] == "NOT_IMPLEMENTED", result


def test_programmatic_2d_mesh_has_2606_edge_descriptors_before_vol_save(tmp_path):
    mesh = NetgenMesh(dim=2)
    mesh.SetMaterial(1, "domain")
    mesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    mesh.SetBCName(0, "outer")
    edge = EdgeDescriptor()
    edge.edgenr = 1
    edge.surfnr = (1, -1)
    edge.domin = 1
    edge.domout = 0
    edge.name = "outer"
    assert mesh.Add(edge) == 1

    points = [
        mesh.Add(MeshPoint(Pnt(x, y, 0.0)))
        for x, y in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    ]
    mesh.Add(Element2D(1, points))
    for first, second in ((0, 1), (1, 2), (2, 3), (3, 0)):
        mesh.Add(Element1D([points[first], points[second]], index=1))

    path = tmp_path / "ngsolve_2606_edge_descriptor.vol"
    mesh.Save(str(path))
    loaded = ng.Mesh(str(path))
    assert loaded.GetMaterials() == ("domain",)
    assert loaded.GetBoundaries() == ("outer",)
