"""
Lint rules for GMSH-related Python scripts.

Each rule receives (filepath, lines) and returns a list of findings.
A finding is a dict with keys: line, severity, rule, message.

Focus: Enforce Radia policy that GMSH is for visualization only,
not mesh generation. Catch common mistakes in .msh handling.
"""

import re
from pathlib import Path
from typing import List, Dict


def check_gmsh_api_mesh_generation(filepath: str, lines: List[str]) -> List[Dict]:
    """CRITICAL: GMSH Python API must not be used for mesh generation."""
    findings = []
    has_gmsh_import = any("import gmsh" in line for line in lines)
    if not has_gmsh_import:
        return findings

    mesh_gen_patterns = [
        (r"gmsh\.model\.occ\.", "gmsh.model.occ.* (OCC geometry creation)"),
        (r"gmsh\.model\.geo\.add", "gmsh.model.geo.add* (geometry creation)"),
        (r"gmsh\.model\.mesh\.generate", "gmsh.model.mesh.generate() (mesh generation)"),
        (r"gmsh\.model\.geo\.synchronize", "gmsh.model.geo.synchronize() (geometry sync)"),
    ]

    for i, line in enumerate(lines, 1):
        stripped = line.split("#")[0]
        for pattern, desc in mesh_gen_patterns:
            if re.search(pattern, stripped):
                findings.append({
                    "line": i,
                    "severity": "CRITICAL",
                    "rule": "gmsh-mesh-generation",
                    "message": (
                        f"GMSH Python API used for mesh generation: {desc}. "
                        "Radia policy: use Netgen (tet) or Cubit (hex) for meshing. "
                        "GMSH is for visualization only."
                    ),
                })
    return findings


def check_gmsh_builder_import(filepath: str, lines: List[str]) -> List[Dict]:
    """CRITICAL: GmshBuilder was removed from the repository."""
    findings = []
    for i, line in enumerate(lines, 1):
        if "GmshBuilder" in line or "gmsh_builder" in line:
            findings.append({
                "line": i,
                "severity": "CRITICAL",
                "rule": "gmsh-builder-removed",
                "message": (
                    "GmshBuilder is removed. Use Netgen or Cubit for mesh generation. "
                    "GMSH is for visualization and post-processing only."
                ),
            })
    return findings


def check_msh_version_mismatch(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: Warn when using ReadGmsh (deprecated, use .vol instead)."""
    findings = []
    for i, line in enumerate(lines, 1):
        stripped = line.split("#")[0]
        if "ReadGmsh" in stripped:
            findings.append({
                "line": i,
                "severity": "HIGH",
                "rule": "readgmsh-deprecated",
                "message": (
                    "ReadGmsh() is deprecated. Use .vol files instead: "
                    "Mesh('model.vol'). Export from Cubit with "
                    "'export netgen \"model.vol\" order N'."
                ),
            })
    return findings


def check_numsubedges_missing(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Suggest NumSubEdges=4 when high-order elements are used."""
    findings = []
    has_curve = any("mesh.Curve(" in line or "Curve(" in line.split("#")[0]
                    for line in lines)
    has_gmsh_post = any("GmshPostExport" in line for line in lines)
    has_numsubedges = (
        any("NumSubEdges" in line for line in lines)
        or _has_numsubedges_companion(filepath)
    )

    if (has_curve or has_gmsh_post) and not has_numsubedges:
        findings.append({
            "line": 0,
            "severity": "MODERATE",
            "rule": "numsubedges-missing",
            "message": (
                "High-order elements detected but Mesh.NumSubEdges not set. "
                "Add NumSubEdges=4 to .geo companion or use "
                "-numsubedges 4 on command line for correct curved display."
            ),
        })
    return findings


def _has_numsubedges_companion(filepath: str) -> bool:
    """Return True if a nearby GMSH display companion sets NumSubEdges."""
    p = Path(filepath)
    candidates = [
        p.with_suffix(".geo"),
        p.with_name(f"{p.stem}_display.geo"),
        p.with_suffix(".msh.opt"),
        p.with_name(f"{p.stem}_display.msh.opt"),
    ]
    for parent in p.parents:
        candidates.extend([
            parent / "_gmsh_display.geo",
            parent / "_gmsh_display.msh.opt",
        ])

    for candidate in dict.fromkeys(candidates):
        try:
            text = candidate.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if "Mesh.NumSubEdges" in text:
            return True
    return False


def check_pip_gmsh_import(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: pip gmsh package should not be used (standalone exe only)."""
    findings = []
    for i, line in enumerate(lines, 1):
        stripped = line.split("#")[0].strip()
        if stripped == "import gmsh" or stripped.startswith("from gmsh "):
            # Allow in test files or post-processing scripts
            if "test_" in filepath or "gmsh_post" in filepath:
                continue
            findings.append({
                "line": i,
                "severity": "HIGH",
                "rule": "pip-gmsh-import",
                "message": (
                    "GMSH Python package (pip install gmsh) should not be used. "
                    "Radia uses standalone gmsh.exe for visualization. "
                    "For mesh input to NGSolve, use .vol files (export netgen)."
                ),
            })
    return findings


def check_meshio_import(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: meshio is removed from the project."""
    findings = []
    for i, line in enumerate(lines, 1):
        stripped = line.split("#")[0].strip()
        if "import meshio" in stripped or "from meshio" in stripped:
            findings.append({
                "line": i,
                "severity": "HIGH",
                "rule": "meshio-removed",
                "message": (
                    "meshio is removed from the project. "
                    "Use .vol files for mesh input to NGSolve, "
                    "or GmshPostExport for .msh output (visualization)."
                ),
            })
    return findings


ALL_RULES = [
    check_gmsh_api_mesh_generation,
    check_gmsh_builder_import,
    check_msh_version_mismatch,
    check_numsubedges_missing,
    check_pip_gmsh_import,
    check_meshio_import,
]
