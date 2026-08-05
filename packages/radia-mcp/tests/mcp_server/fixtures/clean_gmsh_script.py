"""Lint fixture: GMSH-policy-clean viewer launcher.

NEVER executed -- only read as text. Launches the gmsh viewer through
the PATH launcher (allowed), no gmsh Python API, no mesh generation.
"""
import shutil
import subprocess


def open_in_gmsh(msh_path: str) -> None:
    launcher = shutil.which("gmsh")
    if launcher is None:
        raise RuntimeError("gmsh launcher not found on PATH")
    subprocess.Popen([launcher, str(msh_path), "-numsubedges", "4"])
