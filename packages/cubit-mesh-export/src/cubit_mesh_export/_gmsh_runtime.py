"""Gmsh 5 / MSH 4.2 boundary for exporter-owned isolated processes."""
from __future__ import annotations

from email import message_from_bytes
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import re
import tempfile
import zipfile


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


def _payload_identity(gmsh, identity) -> dict:
    archive = os.environ.get("CUBIT_GMSH_WHEEL")
    if not archive:
        raise RuntimeError("Gmsh original wheel missing: set CUBIT_GMSH_WHEEL")
    wheel = Path(archive).resolve(strict=True)
    if hashlib.sha256(wheel.read_bytes()).hexdigest() != identity["wheel_sha256"]:
        raise RuntimeError("Gmsh original wheel SHA256 differs from install report")
    result = {}
    with zipfile.ZipFile(wheel) as contents:
        names = contents.namelist()
        metadatas = [name for name in names if name.endswith(".dist-info/METADATA")]
        if len(metadatas) != 1:
            raise RuntimeError("Gmsh original wheel metadata is ambiguous")
        version = message_from_bytes(contents.read(metadatas[0]))["Version"]
        if version != identity["distribution_version"]:
            raise RuntimeError("Gmsh original wheel version differs from installation")
        for name, path in (("wrapper", Path(gmsh.__file__).resolve()),
                           ("native", Path(gmsh.lib._name).resolve())):
            matches = [member for member in names
                       if Path(member).name.casefold() == path.name.casefold()]
            if len(matches) != 1:
                raise RuntimeError(f"Gmsh {name} is not a unique original-wheel payload")
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != hashlib.sha256(contents.read(matches[0])).hexdigest():
                raise RuntimeError(f"Gmsh loaded {name} differs from original wheel")
            result[name + "_sha256"] = actual
    return result


def probe_gmsh5() -> dict:
    """Check actual native version and a real 4.2 roundtrip in an owned child.

    Call before initializing Gmsh. Never attach to or alter a caller's model.
    The original archive and installer report are required. Archive identity
    is distinct from loaded
    wrapper/native hashes; installed RECORD bytes are not a wheel hash.
    """
    import gmsh

    version = metadata.version("gmsh")
    if version.split(".", 1)[0] != "5":
        raise RuntimeError(f"Gmsh major 5 required, installed {version}")
    if gmsh.isInitialized():
        raise RuntimeError("Gmsh readiness probe requires an uninitialized owned process")
    identity = {"distribution_version": version, **_wheel_identity()}
    identity.update(_payload_identity(gmsh, identity))
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
        identity.update(native_version=native_version, msh_format="4.2",
                        msh42_roundtrip=True)
        return identity
    finally:
        gmsh.finalize()
