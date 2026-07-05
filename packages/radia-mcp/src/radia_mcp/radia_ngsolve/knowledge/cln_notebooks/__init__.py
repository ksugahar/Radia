"""Tanimoto's 3D CLN notebook code (converted from ipynb).

Each .py file is the Python script extraction of the corresponding
Jupyter notebook in public-safe curated corpus

Available notebooks:
  - CLN_AT.py       : A-T formulation, primary 3D Kameari (修論)
  - CLN_T_Omega.py  : T-Ω formulation (修論)
  - CLN_APhi.py     : A-Φ formulation (修論)
  - CLN_2D.py       : 2D scalar reference (修論)
  - A_ICCG_production.py : Latest A + ICCG production code (2024-09-17)
"""

from pathlib import Path

_NOTEBOOK_DIR = Path(__file__).parent

NOTEBOOK_FILES = {
    "AT": "CLN_AT.py",
    "T_Omega": "CLN_T_Omega.py",
    "APhi": "CLN_APhi.py",
    "2D": "CLN_2D.py",
    "production": "A_ICCG_production.py",
}


def get_notebook(name: str) -> str:
    """Retrieve a Tanimoto CLN notebook code as text.

    Args:
        name: One of "AT", "T_Omega", "APhi", "2D", "production".

    Returns:
        Python script content as string. Returns error message if not found.
    """
    if name not in NOTEBOOK_FILES:
        return (f"Unknown notebook '{name}'. Available: "
                f"{', '.join(NOTEBOOK_FILES.keys())}")
    path = _NOTEBOOK_DIR / NOTEBOOK_FILES[name]
    if not path.exists():
        return f"File not found: {path}"
    return path.read_text(encoding="utf-8", errors="replace")


def list_notebooks() -> str:
    """List all available Tanimoto CLN notebooks with sizes."""
    lines = ["# Available Tanimoto 3D CLN Notebooks\n"]
    for key, fname in NOTEBOOK_FILES.items():
        path = _NOTEBOOK_DIR / fname
        if path.exists():
            size = path.stat().st_size
            lines.append(f"- **{key}** ({fname}): {size:,} bytes")
        else:
            lines.append(f"- **{key}** ({fname}): MISSING")
    return "\n".join(lines)
