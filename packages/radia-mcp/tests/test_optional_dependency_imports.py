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
            "numpy", "pandas", "radia", "scipy",
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
