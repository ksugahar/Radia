"""Distribution contract for the Netgen-backed native curver."""

from __future__ import annotations

import json
import tomllib
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
PINNED_NETGEN = "6.2.2607"


def test_distribution_version_and_netgen_abi_contract_are_aligned():
    project = tomllib.loads(
        (PACKAGE / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    dependencies = set(project["dependencies"])
    package_source = (
        PACKAGE / "src" / "cubit_mesh_export" / "__init__.py"
    ).read_text(encoding="utf-8")
    manifest = json.loads(
        (
            PACKAGE
            / "src"
            / "cubit_mesh_export"
            / "native_payloads.json"
        ).read_text(encoding="utf-8")
    )

    assert project["version"] == "2.1.4"
    assert '__version__ = "2.1.4"' in package_source
    assert f"netgen-mesher=={PINNED_NETGEN}" in dependencies
    assert f"ngsolve=={PINNED_NETGEN}" in dependencies
    assert (
        manifest["payloads"]["cubit_mesh_curver.pyd"]["netgen_version"]
        == PINNED_NETGEN
    )
    assert "netgen_version" not in manifest["payloads"]["cubit_mesh_export.ccm"]


def test_readme_states_the_native_abi_boundary_and_supported_runtime():
    readme = (PACKAGE / "README.md").read_text(encoding="utf-8")

    assert f"Netgen/NGSolve {PINNED_NETGEN}" in readme
    assert "cubit_mesh_curver.pyd" in readme
    assert "links against the Netgen C++" in readme
    assert "cubit_mesh_export.ccm" in readme
    assert "does not link against the Python solver wheels" in readme
