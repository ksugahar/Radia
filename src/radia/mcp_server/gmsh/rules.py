"""
Lint rules for GMSH-related Python scripts.

Each rule receives (filepath, lines) and returns a list of findings.
A finding is a dict with keys: line, severity, rule, message.

Focus: Enforce Radia policy that GMSH is for visualization only,
not mesh generation. Catch common mistakes in .msh handling.
"""

import re
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
    """HIGH: Warn when writing v4.1 for NGSolve input (expects v2.2)."""
    findings = []
    for i, line in enumerate(lines, 1):
        stripped = line.split("#")[0]
        if "ReadGmsh" in stripped and "4.1" in stripped:
            findings.append({
                "line": i,
                "severity": "HIGH",
                "rule": "msh-version-mismatch",
                "message": (
                    "NGSolve ReadGmsh() expects .msh v2.2, not v4.1. "
                    "Use GmshPostExport.write() which defaults to v2.2."
                ),
            })
    return findings


def check_numsubedges_missing(filepath: str, lines: List[str]) -> List[Dict]:
    """MODERATE: Suggest NumSubEdges=4 when high-order elements are used."""
    findings = []
    has_curve = any("mesh.Curve(" in line or "Curve(" in line.split("#")[0]
                    for line in lines)
    has_gmsh_post = any("GmshPostExport" in line for line in lines)
    has_numsubedges = any("NumSubEdges" in line for line in lines)

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


def check_pip_gmsh_import(filepath: str, lines: List[str]) -> List[Dict]:
    """HIGH: pip gmsh package should not be used (standalone exe only)."""
    findings = []
    for i, line in enumerate(lines, 1):
        stripped = line.split("#")[0].strip()
        if stripped == "import gmsh" or stripped.startswith("from gmsh "):
            # Allow in test files or gmsh_mesh_import.py (pure reader)
            if "gmsh_mesh_import" in filepath or "test_" in filepath:
                continue
            findings.append({
                "line": i,
                "severity": "HIGH",
                "rule": "pip-gmsh-import",
                "message": (
                    "GMSH Python package (pip install gmsh) should not be used. "
                    "Radia uses standalone gmsh.exe for visualization. "
                    "For .msh reading, use gmsh_mesh_import.py (no gmsh dependency)."
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
                    "Use gmsh_mesh_import.py for .msh reading, "
                    "or GmshPostExport for .msh writing."
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
