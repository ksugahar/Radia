"""Named periodic identification for curved Netgen .vol export.

A periodic export is requested by an explicit nine-field contract; nothing is
inferred.  Master and slave boundary points are paired by a given
translation, curved edge/face coefficients of master entities are copied to
their images (``Mesh.CopyPeriodicCurvedCoefficients``), and one PERIODIC
identification group named ``kt`` is added.  The .vol and its companion JSON
are published as one bundle: the ``.vol`` appears only after the JSON record
was written exclusively and verified.

Everything here works on a Netgen mesh and does not import Cubit.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

FIELDS = ("periodic_name", "periodic_master", "periodic_slave", "periodic_unit",
          "periodic_patches", "periodic_tx", "periodic_ty", "periodic_tz", "periodic_match")
PERIODIC_TYPE = 2  # Netgen Identifications::PERIODIC


class PeriodicError(RuntimeError):
    """The periodic contract or the mesh does not allow a faithful identification."""


@dataclass
class Contract:
    name: str
    master_label: str
    slave_label: str
    unit: str
    translation: tuple[float, float, float]
    match: float
    boundary_pairs: list[tuple[int, int]] = field(default_factory=list)


def parse_patches(text: str) -> list[tuple[int, int]]:
    """Parse ``"m1:s1,m2:s2"``; digits only, no empty entries."""
    if not text or text.endswith(","):
        raise PeriodicError("Empty BC patch entry")
    pairs = []
    for patch in text.split(","):
        parts = patch.split(":")
        if len(parts) != 2 or not parts[0] or not parts[1]:
            raise PeriodicError("Expected BC patch list master:slave,...")
        if not (parts[0].isascii() and parts[0].isdigit() and parts[1].isascii() and parts[1].isdigit()):
            raise PeriodicError("BC patch IDs must contain digits only")
        pairs.append((int(parts[0]), int(parts[1])))
    return pairs


def parse_contract(values: dict[str, Any]) -> Contract | None:
    """Build the contract from the nine fields; ``None`` when none is given."""
    present = [k for k in FIELDS if values.get(k) is not None]
    if not present:
        return None
    if len(present) != len(FIELDS):
        raise PeriodicError("Named periodic export requires all nine explicit periodic fields.")
    return Contract(name=str(values["periodic_name"]), master_label=str(values["periodic_master"]),
                    slave_label=str(values["periodic_slave"]), unit=str(values["periodic_unit"]),
                    translation=(float(values["periodic_tx"]), float(values["periodic_ty"]),
                                 float(values["periodic_tz"])),
                    match=float(values["periodic_match"]),
                    boundary_pairs=parse_patches(str(values["periodic_patches"])))


def validate_contract(contract: Contract) -> tuple[dict[int, int], set[int]]:
    if (contract.name != "kt" or contract.unit != "m" or not contract.master_label
            or not contract.slave_label or contract.master_label == contract.slave_label
            or not math.isfinite(contract.match) or contract.match <= 0
            or not contract.boundary_pairs):
        raise PeriodicError("Explicit named periodic contract is incomplete")
    if not all(math.isfinite(t) for t in contract.translation):
        raise PeriodicError("Nonfinite periodic translation")
    master_to_slave: dict[int, int] = {}
    slaves: set[int] = set()
    for m, s in contract.boundary_pairs:
        if m < 1 or s < 1 or m == s or m in master_to_slave or s in slaves:
            raise PeriodicError("Periodic patch map is not one-to-one")
        master_to_slave[m] = s
        slaves.add(s)
    if any(m in slaves for m in master_to_slave):
        raise PeriodicError("Periodic master and slave BC sets overlap")
    return master_to_slave, slaves


def geometry_serialization(text: str) -> str:
    """The saved mesh without its identification sections."""
    lines = text.split("\n")
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        if line in ("identifications", "identificationtypes", "identificationnames"):
            if i + 1 >= len(lines):
                raise PeriodicError("Truncated periodic metadata count")
            try:
                count = int(lines[i + 1].split()[0])
            except (ValueError, IndexError) as error:
                raise PeriodicError("Invalid periodic metadata count") from error
            if count < 0:
                raise PeriodicError("Invalid periodic metadata count")
            rows = count if line == "identifications" else 1
            if line != "identifications" and count != 1:
                raise PeriodicError("Only one named periodic group is supported")
            if i + 1 + rows >= len(lines):
                raise PeriodicError("Truncated periodic metadata rows")
            i += 2 + rows
            continue
        out.append(line)
        i += 1
    return "\n".join(out)


def _labels(mesh) -> tuple[dict[int, int], dict[int, str]]:
    descriptor_to_bc: dict[int, int] = {}
    labels: dict[int, str] = {}
    for fd in range(1, mesh.GetNFaceDescriptors() + 1):
        bc = mesh.FaceDescriptor(fd).bc
        if bc < 1:
            raise PeriodicError("Invalid boundary-condition property")
        descriptor_to_bc[fd] = bc
        name = mesh.GetBCName(fd - 1)
        if bc in labels and labels[bc] != name:
            raise PeriodicError(f"Conflicting descriptor labels for BC {bc}")
        labels[bc] = name
    return descriptor_to_bc, labels


def triangle_plan(mesh, master_to_slave: dict[int, int], points: dict[int, int],
                  descriptor_to_bc: dict[int, int]) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Pair every sorted master triangle with its mapped image, both directions checked."""
    selected = set(master_to_slave) | set(master_to_slave.values())
    patches: dict[int, set[tuple[int, ...]]] = {}
    for el in mesh.Elements2D():
        bc = descriptor_to_bc[el.index]
        if bc not in selected:
            continue
        if str(el.type).rsplit(".", 1)[-1] != "TRIG" or len(el.vertices) != 3:
            raise PeriodicError("Named periodic geometry requires triangular linear connectivity")
        tri = tuple(sorted(v.nr for v in el.vertices))
        if tri[0] == tri[1] or tri[1] == tri[2] or tri in patches.setdefault(bc, set()):
            raise PeriodicError("Duplicate or degenerate periodic triangle")
        patches[bc].add(tri)
    plan = []
    for mbc, sbc in master_to_slave.items():
        if mbc not in patches or sbc not in patches:
            raise PeriodicError(f"Missing periodic BC patch {mbc if mbc not in patches else sbc}")
        covered: set[tuple[int, ...]] = set()
        for master in sorted(patches[mbc]):
            try:
                mapped = tuple(points[p] for p in master)
            except KeyError as error:
                raise PeriodicError(f"Missing periodic point {error.args[0]}") from None
            key = tuple(sorted(mapped))
            if key not in patches[sbc] or key in covered:
                raise PeriodicError("Periodic triangle correspondence is incomplete or ambiguous")
            covered.add(key)
            plan.append((master, mapped))
        if covered != patches[sbc]:
            raise PeriodicError("Incomplete reverse periodic triangle coverage")
    return plan


def apply_identification(mesh, contract: Contract,
                         save_text: Callable[[Any], str]) -> dict[str, Any]:
    """Pair points, synchronize curved seams and add the ``kt`` group.

    ``save_text(mesh)`` must return the mesh serialization; it is used to
    prove that adding the identification leaves the geometry unchanged.
    """
    master_to_slave, slave_bcs = validate_contract(contract)
    if mesh.GetNrIdentifications() != 0:
        raise PeriodicError("Pre-existing periodic metadata is unsupported")
    descriptor_to_bc, labels = _labels(mesh)
    for bc, label in labels.items():
        if ((label == contract.master_label and bc not in master_to_slave)
                or (label == contract.slave_label and bc not in slave_bcs)):
            raise PeriodicError("Periodic patch map omits a labeled boundary")

    coords = {i + 1: tuple(p.p) for i, p in enumerate(mesh.Points())}
    boundary_points: dict[int, set[int]] = {}
    master_incidence: dict[int, set[int]] = {}
    for el in mesh.Elements2D():
        bc = descriptor_to_bc.get(el.index)
        if bc is None:
            raise PeriodicError("Surface element references an unknown descriptor")
        is_master = bc in master_to_slave
        if not is_master and bc not in slave_bcs:
            continue
        if labels[bc] != (contract.master_label if is_master else contract.slave_label):
            raise PeriodicError(f"Periodic BC label mismatch at BC {bc}")
        for v in el.vertices:
            pi = v.nr
            if pi < 1 or pi > len(coords):
                raise PeriodicError("Periodic surface point index out of range")
            if not all(math.isfinite(c) for c in coords[pi]):
                raise PeriodicError("Nonfinite periodic point coordinate")
            boundary_points.setdefault(bc, set()).add(pi)
            if is_master:
                master_incidence.setdefault(pi, set()).add(bc)
    masters: set[int] = set()
    slaves: set[int] = set()
    for mbc, sbc in master_to_slave.items():
        if not boundary_points.get(mbc) or not boundary_points.get(sbc):
            raise PeriodicError(f"Empty periodic patch at master BC {mbc}")
        masters |= boundary_points[mbc]
        slaves |= boundary_points[sbc]
    if len(masters) != len(slaves):
        raise PeriodicError("Periodic point-set cardinality mismatch")
    if masters & slaves:
        raise PeriodicError("Periodic point index belongs to both sides")

    tx, ty, tz = contract.translation
    pairs: list[tuple[int, int]] = []
    assigned: set[int] = set()
    max_distance = 0.0
    slave_list = sorted(slaves)
    for m in sorted(masters):
        mx, my, mz = coords[m]
        candidates = []
        for s in slave_list:
            if not all(s in boundary_points[master_to_slave[bc]] for bc in master_incidence[m]):
                continue
            sx, sy, sz = coords[s]
            d = math.hypot(sx - (mx + tx), sy - (my + ty), sz - (mz + tz))
            if d <= contract.match:
                candidates.append((s, d))
        if len(candidates) != 1:
            raise PeriodicError(f"Periodic master point {m} has {len(candidates)} candidates")
        s, d = candidates[0]
        if s in assigned:
            raise PeriodicError("Periodic target point assigned more than once")
        assigned.add(s)
        pairs.append((m, s))
        max_distance = max(max_distance, d)
    if assigned != slaves:
        raise PeriodicError("Incomplete periodic slave coverage")
    mapping = dict(pairs)
    for mbc, sbc in master_to_slave.items():
        if {mapping[p] for p in boundary_points[mbc]} != boundary_points[sbc]:
            raise PeriodicError(f"Incomplete periodic patch coverage at master BC {mbc}")

    plan = triangle_plan(mesh, master_to_slave, mapping, descriptor_to_bc)
    if any(str(el.type).rsplit(".", 1)[-1] != "TET" or len(el.vertices) != 4
           for el in mesh.Elements3D()):
        raise PeriodicError("Periodic export requires linear tetrahedral connectivity")
    try:
        transfers = mesh.CopyPeriodicCurvedCoefficients([list(m) for m, _ in plan],
                                                        [list(s) for _, s in plan])
    except Exception as error:  # Netgen raises a plain exception type
        raise PeriodicError(f"Curved periodic synchronization failed: {error}") from error

    before = geometry_serialization(save_text(mesh))
    nr = 0
    for m, s in pairs:
        nr = mesh.AddNamedPointIdentification(m, s, contract.name)
    if (mesh.GetNrIdentifications() != 1 or nr != 1
            or "PERIODIC" not in str(mesh.GetIdentificationType(1)).upper()):
        raise PeriodicError("Named periodic metadata insertion failed")
    if geometry_serialization(save_text(mesh)) != before:
        raise PeriodicError("Periodic insertion changed serialized geometry")
    return {
        "pairs": pairs,
        "descriptor_to_bc": descriptor_to_bc,
        "max_distance": max_distance,
        "curved_transfers": [dict(t) for t in transfers],
    }


def companion_record(result: dict[str, Any], descriptor_surface: dict[int, int]) -> dict[str, Any]:
    return {
        "name": "kt", "nr": 1, "type": PERIODIC_TYPE,
        "pair_count": len(result["pairs"]),
        "max_distance": result["max_distance"],
        "descriptor_mapping": [{"fd": fd, "bc": bc, "cubit_surface": descriptor_surface[fd]}
                               for fd, bc in sorted(result["descriptor_to_bc"].items())],
        "curved_transfers": [{"entity": t["entity"], "master_topology_index": t["master"],
                              "slave_topology_index": t["slave"], "mode": t["mode"],
                              "sign": t["sign"], "original": list(t["original"]),
                              "expected": list(t["expected"])}
                             for t in result["curved_transfers"]],
    }


def output_paths(final: Path) -> tuple[Path, Path, Path]:
    if final.suffix != ".vol":
        raise PeriodicError("Named periodic output requires an explicit .vol suffix")
    return final, final.with_name(final.name + ".partial.vol"), final.with_name(final.name + ".json")


def require_fresh(final: Path) -> None:
    for path in output_paths(final):
        if os.path.lexists(path):
            raise PeriodicError(
                f"Named periodic output already exists: {path}. Inspect the earlier failure; "
                "remove only your own unqualified leftovers before retrying.")


def publish_bundle(final: Path, vol_bytes: bytes, record_bytes: bytes) -> None:
    """Stage the .vol, exclusively write the JSON, then publish the .vol without replacing."""
    final, partial, sidecar = output_paths(final)
    require_fresh(final)
    if not record_bytes:
        raise PeriodicError("Empty periodic sidecar record")
    with open(partial, "xb") as handle:
        handle.write(vol_bytes)
        handle.flush()
        os.fsync(handle.fileno())
    if partial.read_bytes() != vol_bytes:
        raise PeriodicError("Periodic temporary serialization is incomplete or differs")
    with open(sidecar, "xb") as handle:
        handle.write(record_bytes)
        handle.flush()
        os.fsync(handle.fileno())
        identity = os.fstat(handle.fileno())
    # Keep a read handle on the record while publishing the .vol marker.
    with open(sidecar, "rb") as check:
        info = os.fstat(check.fileno())
        if (info.st_ino, info.st_dev) != (identity.st_ino, identity.st_dev) or check.read() != record_bytes:
            raise PeriodicError("Periodic bundle record content changed")
        if os.name == "nt":
            os.rename(partial, final)  # fails when the target exists
        else:
            os.link(partial, final)
            os.unlink(partial)
