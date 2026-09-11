"""Load a Hiruma benchmark mesh from its Netgen .vol file.

The Hiruma meshes (``mesh1_*.vol``) are not distributed with the repository
(``*.vol`` is gitignored).  Place them next to the benchmark scripts; a
missing mesh raises ``FileNotFoundError`` listing the meshes that are present.
"""
from pathlib import Path

from ngsolve import Mesh

MESH_DIR = Path(__file__).resolve().parent


def mesh_path(name):
    """Return the .vol path for a mesh name such as ``mesh1_3.5T``."""
    if name.endswith(".msh"):
        raise ValueError(f"{name}: the Hiruma benchmarks read Netgen .vol meshes; "
                         "pass the mesh name without an extension")
    stem = name[:-len(".vol")] if name.endswith(".vol") else name
    return MESH_DIR / f"{stem}.vol"


def load_hiruma_mesh(name):
    """Load ``<name>.vol`` from this directory as an NGSolve mesh."""
    path = mesh_path(name)
    if not path.is_file():
        available = sorted(p.stem for p in MESH_DIR.glob("mesh1_*.vol"))
        raise FileNotFoundError(
            f"{path.name} not found in {MESH_DIR}.  The Hiruma meshes are not "
            f"distributed with the repository; place mesh1_*.vol next to the "
            f"benchmark scripts.  Available: {available}")
    return Mesh(str(path))
