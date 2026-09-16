from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path


def test_fem_force_and_motor_servers_import_without_solver_dependencies() -> None:
    script = textwrap.dedent(
        """
        import importlib.abc
        import sys

        sys.path.insert(0, sys.argv[1])
        blocked = {
            "build123d", "gmsh", "matplotlib", "netgen", "ngsolve",
            "numpy", "pandas", "radia", "scipy", "cubit_mesh_export", "cae_mcp_core",
        }

        class OptionalDependencyBlocker(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname.split(".", 1)[0] in blocked:
                    raise ModuleNotFoundError(
                        f"blocked optional dependency: {fullname}", name=fullname
                    )

        sys.meta_path.insert(0, OptionalDependencyBlocker())
        import radia_mcp.fem.server
        import radia_mcp.force.server
        import radia_mcp.motor.server
        from radia_mcp.common import examples, web_docs
        assert set(examples.FAMILIES) == {'build123d', 'gmsh'}
        assert not hasattr(examples, 'refresh_cubit_examples')
        assert examples._resolve_family('cubit') == []
        assert examples._resolve_family('cubit_local') == []
        assert set(web_docs.DOCS_INDEX) == {'build123d'}
        assert not hasattr(web_docs, 'search_forum')
        """
    )
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            script,
            str(Path(__file__).resolve().parents[1] / "src"),
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
