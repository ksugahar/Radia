"""Native Gmsh readback of an MSH file: does the viewer actually load it?

The pure-Python inspector (`msh_inspect`) checks MSH 4.2 structure; this
check asks native Gmsh 5 itself, in a child process with the qualified
runtime, to open the file and reports what it loaded: node count, element
counts per dimension, physical groups, nodes that no element references and
any Gmsh error or warning. A file that the Python parser accepts but native
Gmsh rejects or reads with warnings fails here instead of in the GUI.
"""
from __future__ import annotations

import hashlib
import shutil
import tempfile
from pathlib import Path
from typing import Any

from radia_mcp.gmsh._gmsh_subprocess import run_gmsh_json_subprocess

_SCRIPT = r'''
import json, sys
from pathlib import Path
import gmsh

path, out = sys.argv[1], sys.argv[-1]
result = {"ok": False, "ran": True, "path": path}
gmsh.initialize(["radia-native-readback"], readConfigFiles=False)
try:
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.logger.start()
    try:
        gmsh.open(path)
    except Exception as exc:
        result["error"] = f"native gmsh.open failed: {exc}"
    else:
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        nodes = {int(tag) for tag in node_tags}
        referenced = set()
        elements = {}
        for dim in range(4):
            _, element_tags, element_nodes = gmsh.model.mesh.getElements(dim)
            count = sum(len(tags) for tags in element_tags)
            if count:
                elements[str(dim)] = count
            for block in element_nodes:
                referenced.update(int(tag) for tag in block)
        groups = sorted([dim, gmsh.model.getPhysicalName(dim, tag)]
                        for dim, tag in gmsh.model.getPhysicalGroups())
        result.update(nodes=len(nodes), elements=elements, physical_groups=groups,
                      unreferenced_nodes=len(nodes - referenced),
                      missing_nodes=len(referenced - nodes),
                      views=len(gmsh.view.getTags()))
        result["ok"] = True
    log = gmsh.logger.get()
    gmsh.logger.stop()
    result["errors"] = [m for m in log if m.startswith("Error")][:20]
    result["warnings"] = [m for m in log if m.startswith("Warning")][:20]
finally:
    gmsh.finalize()
if result["ok"] and (result["errors"] or result["warnings"] or result["unreferenced_nodes"]
                     or result["missing_nodes"]):
    result["ok"] = False
    result["error"] = "native Gmsh loaded the file with errors, warnings or unconnected nodes"
Path(out).write_text(json.dumps(result), encoding="utf-8")
'''


def native_readback(path: str | Path, *, timeout_s: float = 120.0) -> dict[str, Any]:
    """Open ``path`` with native Gmsh 5 in a child process and report what loaded.

    ``ok`` is True only when Gmsh opened the file with no error or warning,
    every node is referenced by an element and every referenced node exists.
    ``ran`` is False when the runtime is missing or unqualified.
    """
    resolved = Path(path).resolve()
    if not resolved.is_file():
        return {"ok": False, "ran": False, "path": str(resolved), "error": f"file not found: {resolved}"}
    # Native Gmsh also reads companions such as FILE.msh.opt beside the file it
    # opens. Read an isolated byte copy so the MSH bytes are the only input, as
    # in the viewer, which never receives the exporter's companions.
    with tempfile.TemporaryDirectory(prefix="radia_mcp_gmsh_readback_input_") as work:
        isolated = Path(work) / "readback.msh"
        shutil.copyfile(resolved, isolated)
        digest = _sha256(isolated)
        if digest != _sha256(resolved):
            return {"ok": False, "ran": False, "path": str(resolved),
                    "error": "MSH changed while it was copied for native readback"}
        result = run_gmsh_json_subprocess(_SCRIPT, [str(isolated)], timeout_s=timeout_s,
                                          prefix="radia_mcp_gmsh_readback_")
    return {**result, "path": str(resolved), "sha256": digest}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_native_readback(path: str | Path, *, expected_groups=None, expected_elements=None,
                            timeout_s: float = 120.0) -> dict[str, Any]:
    """Raise unless native Gmsh loads ``path`` cleanly and matches the expectations.

    ``expected_groups`` is a set of ``(dim, name)``; ``expected_elements`` maps
    a dimension to its element count. Returns the readback summary.
    """
    result = native_readback(path, timeout_s=timeout_s)
    if not result.get("ok"):
        detail = {key: result.get(key) for key in
                  ("error", "errors", "warnings", "unreferenced_nodes", "missing_nodes")}
        raise RuntimeError(f"Native Gmsh readback failed for {path}: {detail}")
    if expected_groups is not None:
        observed = {(int(dim), name) for dim, name in result["physical_groups"]}
        if observed != set(expected_groups):
            raise RuntimeError(f"Native Gmsh physical groups {sorted(observed)} != expected {sorted(expected_groups)}")
    if expected_elements is not None:
        observed = {int(dim): count for dim, count in result["elements"].items()}
        if observed != dict(expected_elements):
            raise RuntimeError(f"Native Gmsh element counts {observed} != expected {dict(expected_elements)}")
    return {key: result[key] for key in ("nodes", "elements", "physical_groups", "unreferenced_nodes", "views",
                                         "sha256")}
