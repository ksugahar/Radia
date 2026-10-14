"""Load a Hiruma benchmark mesh from its Netgen ``.vol`` fixture."""
from pathlib import Path

from ngsolve import Mesh

MESH_DIR = Path(__file__).resolve().parent / "meshes"


def mesh_path(name, mesh_dir=None):
    """Return the .vol path for a mesh name such as ``mesh1_3.5T``."""
    name = str(name)
    if name.endswith(".msh"):
        raise ValueError(f"{name}: the Hiruma benchmarks read Netgen .vol meshes; "
                         "pass the mesh name without an extension")
    stem = name[:-len(".vol")] if name.endswith(".vol") else name
    return Path(mesh_dir or MESH_DIR) / f"{stem}.vol"


def load_hiruma_mesh(name, mesh_dir=None):
    """Load ``<name>.vol`` from the validation fixture directory."""
    path = mesh_path(name, mesh_dir)
    if not path.is_file():
        available = sorted(p.stem for p in path.parent.glob("mesh1_*.vol"))
        raise FileNotFoundError(
            f"{path.name} not found in {path.parent}.  The tracked validation "
            "fixture is mesh1_2.5T.vol; place optional larger mesh1_*.vol "
            f"fixtures in the same directory.  Available: {available}")
    return Mesh(str(path))
