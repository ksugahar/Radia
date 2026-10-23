"""Gmsh 5 / MSH 4.2 boundary for exporter-owned isolated processes."""
from __future__ import annotations

import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import re
import tempfile


def require_msh42(path: str | Path) -> None:
    """Reject other formats before delegating to Gmsh's permissive reader."""
    with Path(path).open("rb") as stream:
        marker = stream.readline().strip()
        row = stream.readline().split()
        end = stream.readline().strip()
    if (marker != b"$MeshFormat" or len(row) != 3 or row[0] != b"4.2"
            or row[1] not in (b"0", b"1") or row[2] != b"8"
            or (row[1] == b"0" and end != b"$EndMeshFormat")):
        raise ValueError("expected Gmsh MSH 4.2 with 8-byte data")


def _wheel_identity() -> dict:
    report = os.environ.get("CUBIT_GMSH_INSTALL_REPORT")
    if not report:
        raise RuntimeError("Gmsh wheel identity missing: install with --report and set CUBIT_GMSH_INSTALL_REPORT")
    document = json.loads(Path(report).read_text(encoding="utf-8"))
    entries = [row for row in document["install"]
               if row["metadata"]["name"].lower() == "gmsh"]
    if len(entries) != 1 or entries[0]["metadata"]["version"] != metadata.version("gmsh"):
        raise RuntimeError("Gmsh install report does not match installed distribution")
    digest = entries[0]["download_info"]["archive_info"]["hashes"]["sha256"]
    if not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise RuntimeError("Gmsh install report lacks a full wheel SHA256")
    return {"wheel_sha256": digest, "wheel_provenance": "pip install report"}


def probe_gmsh5() -> dict:
    """Check actual native version and a real 4.2 roundtrip in an owned child.

    Call before initializing Gmsh. Never attach to or alter a caller's model.
    Archive identity, when retained by the installer, is distinct from loaded
    wrapper/native hashes; installed RECORD bytes are not a wheel hash.
    """
    import gmsh

    version = metadata.version("gmsh")
    if version.split(".", 1)[0] != "5":
        raise RuntimeError(f"Gmsh major 5 required, installed {version}")
    if gmsh.isInitialized():
        raise RuntimeError("Gmsh readiness probe requires an uninitialized owned process")
    identity = {"distribution_version": version, **_wheel_identity()}
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        native_version = gmsh.option.getString("General.Version")
        if not native_version.startswith("5."):
            raise RuntimeError(f"Gmsh native major 5 required, loaded {native_version}")
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.model.add("cubit_msh42_probe")
        gmsh.model.addDiscreteEntity(2, 1)
        gmsh.model.mesh.addNodes(2, 1, [1, 2, 3], [0, 0, 0, 1, 0, 0, 0, 1, 0])
        gmsh.model.mesh.addElementsByType(1, 2, [1], [1, 2, 3])
        with tempfile.TemporaryDirectory(prefix="cubit_msh42_probe_") as work:
            path = Path(work) / "probe.msh"
            gmsh.write(str(path))
            require_msh42(path)
            gmsh.clear()
            gmsh.open(str(path))
            tags, coordinates, _ = gmsh.model.mesh.getNodes()
            etags, connectivity = gmsh.model.mesh.getElementsByType(2)
            if (list(tags) != [1, 2, 3] or list(coordinates) != [0, 0, 0, 1, 0, 0, 0, 1, 0]
                    or list(etags) != [1] or list(connectivity) != [1, 2, 3]):
                raise RuntimeError("Gmsh MSH 4.2 read/write probe changed mesh data")
        wrapper = Path(gmsh.__file__).resolve()
        native = Path(gmsh.lib._name).resolve()
        identity.update(native_version=native_version, msh_format="4.2",
                        msh42_roundtrip=True,
                        wrapper_sha256=hashlib.sha256(wrapper.read_bytes()).hexdigest(),
                        native_sha256=hashlib.sha256(native.read_bytes()).hexdigest())
        return identity
    finally:
        gmsh.finalize()
