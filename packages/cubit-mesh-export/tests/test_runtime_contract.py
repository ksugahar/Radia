"""Distribution contract for the Netgen-backed native curver."""

from __future__ import annotations

import json
import tomllib
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
PINNED_NETGEN = "6.2.2608"


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

    assert project["version"] == "2.2.8"
    assert '__version__ = "2.2.8"' in package_source
    assert f"netgen-mesher=={PINNED_NETGEN}" in dependencies
    assert f"ngsolve=={PINNED_NETGEN}" in dependencies
    assert (
        manifest["payloads"]["cubit_mesh_curver.pyd"]["netgen_version"]
        == PINNED_NETGEN
    )
    assert "netgen_version" not in manifest["payloads"]["cubit_mesh_export.ccm"]

    # Import the real mandatory payload in the selected external Netgen runtime.
    # Metadata agreement alone cannot detect an extension built for an old ABI.
    import importlib.metadata as metadata
    from cubit_mesh_export import cubit_mesh_curver

    assert metadata.version("netgen-mesher") == PINNED_NETGEN
    assert metadata.version("ngsolve") == PINNED_NETGEN
    assert cubit_mesh_curver.__version__ == project["version"]
    assert cubit_mesh_curver.has_netgen is True


def test_readme_states_the_native_abi_boundary_and_supported_runtime():
    readme = (PACKAGE / "README.md").read_text(encoding="utf-8")

    assert f"Netgen/NGSolve {PINNED_NETGEN}" in readme
    assert "cubit_mesh_curver.pyd" in readme
    assert "links against the Netgen C++" in readme
    assert "cubit_mesh_export.ccm" in readme
    assert "does not link against the Python solver wheels" in readme


def test_native_curver_splits_all_quad_faces_with_netgen_typed_indices():
    import netgen.meshing
    import ngsolve
    import numpy as np
    from cubit_mesh_export import cubit_mesh_curver

    points = np.array([
        [0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
        [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.],
    ])
    faces = [([1, 4, 3, 2], (0, 0, -1)), ([5, 6, 7, 8], (0, 0, 1)),
             ([1, 2, 6, 5], (0, -1, 0)), ([2, 3, 7, 6], (1, 0, 0)),
             ([3, 4, 8, 7], (0, 1, 0)), ([4, 1, 5, 8], (-1, 0, 0))]
    surfaces = [dict(id=i, tris=[], quads=quad, uvs={})
                for i, (quad, _) in enumerate(faces, 1)]

    def project(surface, x, y, z, u, v, has_hint):
        return x, y, z, 0., 0.

    def normal(surface, x, y, z):
        return faces[surface - 1][1]

    with ngsolve.TaskManager():
        mesh = cubit_mesh_curver.build_curved_mesh(
            points, np.arange(1, 9, dtype=np.int32),
            [(1, list(range(1, 9)), 1)], surfaces, project, normal,
            order=2, surface_only=True, split_quads=True,
        )
    triangles = list(mesh.Elements2D())
    assert len(list(mesh.Elements3D())) == 0
    assert len(triangles) == 12
    assert all(len(element.vertices) == 3 for element in triangles)
    assert sorted(element.index for element in triangles) == [i for i in range(1, 7) for _ in range(2)]
