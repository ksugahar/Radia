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
_mesh_helpers = pytest.importorskip("ngsolve.meshes")
MakeStructured2DMesh = _mesh_helpers.MakeStructured2DMesh
MakeStructured3DMesh = _mesh_helpers.MakeStructured3DMesh

from radia.vim._vim import _curve_mesh  # noqa: E402


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


def test_self_hosted_ci_keeps_the_ngsolve_abi_in_a_run_local_environment():
    workflow = (ROOT / ".github" / "workflows" / "build-test.yml").read_text(
        encoding="utf-8"
    )
    assert "Create isolated NGSolve environment" in workflow
    assert "-m venv --system-site-packages" in workflow
    assert "RADIA_CI_VENV=" in workflow
    assert "NGSOLVE_DIR=" in workflow
    assert "Netgen_DIR=" in workflow
    assert "ngsolve.__version__==want['ngsolve']" in workflow
    assert "Get-Command git -ErrorAction Stop" in workflow
    assert "RADIA_CI_GIT_DIR=$gitDir" in workflow
    assert '$env:PATH = "$env:RADIA_CI_GIT_DIR;$env:PATH"' in workflow
    assert "git is absent from pytest PATH" in workflow
    assert "git --version" in workflow
    assert "shutil.which('git')" in workflow
    assert "[Text.UTF8Encoding]::new($false)" in workflow


def test_build_scripts_resolve_netgen_from_the_active_python_environment():
    for relative in ("Build.ps1", "tools/_build_cubit_plugin.ps1"):
        script = (ROOT / relative).read_text(encoding="utf-8")
        assert "Get-Command python -ErrorAction Stop" in script, relative
        assert "import netgen,os" in script, relative
        assert "-DPython3_EXECUTABLE=" in script, relative
        assert "-Dpybind11_DIR=" in script, relative
        assert (
            "C:\\Program Files\\Python312\\Lib\\site-packages\\netgen"
            not in script
        ), relative

    build_script = (ROOT / "Build.ps1").read_text(encoding="utf-8")
    assert build_script.count("-DPython3_EXECUTABLE=") >= 3

    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "pybind11.get_cmake_dir()" in cmake
    assert "set(pybind11_DIR ${PYTHON_SITE_PACKAGES}" not in cmake
    assert "target_include_directories(radia_mex BEFORE PRIVATE" in cmake
    mex_includes = cmake[
        cmake.index("target_include_directories(radia_mex BEFORE PRIVATE") :
        cmake.index("target_compile_definitions(radia_mex PRIVATE")
    ]
    assert "${Python3_INCLUDE_DIRS}" not in mex_includes


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


@pytest.mark.parametrize(
    ("family", "factory"),
    [
        ("QUAD", lambda: MakeStructured2DMesh(quads=True, nx=1, ny=1)),
        ("HEX", lambda: MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)),
        (
            "WEDGE",
            lambda: MakeStructured3DMesh(
                hexes=False, prism=True, nx=1, ny=1, nz=1
            ),
        ),
    ],
)
def test_programmatic_structured_meshes_are_curve_safe(family, factory):
    mesh = factory()
    assert mesh.ngmesh.EdgeDescriptors() == [], family
    _curve_mesh(mesh, 2)
    assert mesh.GetCurveOrder() == 2, family
    assert mesh.ngmesh.EdgeDescriptors(), family


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
