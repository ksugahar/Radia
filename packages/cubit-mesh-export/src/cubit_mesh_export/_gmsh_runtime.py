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


NATIVE_FIELD_TAG_MAX = 2**31 - 1


def validate_native_field_tags(tags, association):
    """Reject unsupported native field indices; never renumber or substitute."""
    for tag in tags:
        if int(tag) > NATIVE_FIELD_TAG_MAX:
            raise ValueError(
                f'Gmsh native {association} field tag {tag} exceeds '
                f'{NATIVE_FIELD_TAG_MAX}; large-tag native field data is '
                'unsupported. Re-export with compact producer numbering; '
                'tags are not remapped and no fallback is used. '
                'Pure-Python readers and mesh-only tags are unaffected.')


def validate_native_msh_fields(path: str | Path) -> None:
    """Check ASCII MSH 4.2 field indices before ANY native initialization.

    Mesh-only APIs accept direct MSH files, never GEO scripts or other formats.
    Only field indices have a signed-int native limit; mesh-only size_t tags
    and the independent pure-Python inventory reader remain unrestricted.
    """
    path = Path(path)
    require_msh42(path)
    text = path.read_text(encoding="utf-8")
    rows = text.splitlines()
    if len(rows) < 3 or rows[0].strip() != '$MeshFormat' or rows[1].split() != ['4.2', '0', '8']:
        raise ValueError(f'{path}: native field hand-off requires ASCII MSH 4.2 (size 8)')
    lines = iter(rows)
    for line in lines:
        kind = line.strip()
        if kind not in ('$NodeData', '$ElementData', '$ElementNodeData'):
            continue
        association = 'node' if kind == '$NodeData' else 'element'
        try:
            strings = int(next(lines))
            if strings < 0:
                raise ValueError('negative string tag count')
            for _ in range(strings):
                next(lines)
            body = []
            end = kind.replace('$', '$End')
            for row in lines:
                if row.strip() == end:
                    break
                body.append(row)
            else:
                raise ValueError(f'missing {end}')
            tokens = iter(' '.join(body).split())
            nreal = int(next(tokens))
            if nreal < 0:
                raise ValueError('negative real tag count')
            for _ in range(nreal):
                next(tokens)
            ninteger = int(next(tokens))
            if ninteger < 3:
                raise ValueError('field requires three integer tags')
            integers = [int(next(tokens)) for _ in range(ninteger)]
            ncomp, declared = integers[1:3]
            blocks, count, minimum, maximum = [int(next(tokens)) for _ in range(4)]
            if ncomp <= 0 or min(blocks, count, declared) < 0 or count != declared:
                raise ValueError('invalid field count/component header')
            if count:
                validate_native_field_tags((maximum,), association)
            consumed = 0
            for _ in range(blocks):
                width = ncomp
                if kind == '$ElementNodeData':
                    int(next(tokens))  # element type
                    width = int(next(tokens))
                n = int(next(tokens))
                if n < 0 or n > count - consumed or width <= 0:
                    raise ValueError('invalid field block width/count')
                for _ in range(n):
                    validate_native_field_tags((int(next(tokens)),), association)
                for _ in range(n * width):
                    next(tokens)
                consumed += n
            if consumed != count or next(tokens, None) is not None:
                raise ValueError('field block count/trailing token mismatch')
        except StopIteration as exc:
            raise ValueError(f'{path}: truncated {kind} before native hand-off') from exc


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
