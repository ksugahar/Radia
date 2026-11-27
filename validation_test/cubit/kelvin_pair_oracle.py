"""Independent oracle for the Kelvin pairs the Cubit export writes into a .vol.

Strips the identification sections from a copy of the .vol, pairs that copy
after load with ``radia.kelvin_identify_ngsolve.add_kelvin_identification``
(offset detected from the mesh, assignment by distance), and compares the
result with the pairs stored in the original file.
"""
from __future__ import annotations

import os

_SECTIONS = ("identifications", "identificationtypes", "identificationnames")


def without_identifications(vol: str, out_path: str) -> str:
    """Write ``vol`` without its identification sections to ``out_path``."""
    with open(vol, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    kept, i = [], 0
    while i < len(lines):
        if lines[i].strip() in _SECTIONS:
            i += 2 + int(lines[i + 1].split()[0])
            continue
        kept.append(lines[i])
        i += 1
    with open(out_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(kept) + "\n")
    return out_path


def pair_set(mesh) -> set[frozenset[int]]:
    return {frozenset((int(row[0].nr), int(row[1].nr)))
            for row in mesh.ngmesh.GetIdentifications()}


def compare(vol: str, work_dir: str) -> dict:
    """Exported pairs of ``vol`` against the post-hoc pairing of its bare copy."""
    from ngsolve import Mesh
    from radia.kelvin_identify_ngsolve import (
        add_kelvin_identification, has_kelvin_identification)

    mesh = Mesh(vol)
    bare_path = without_identifications(
        vol, os.path.join(work_dir, os.path.basename(vol) + ".bare.vol"))
    bare = Mesh(bare_path)
    bare_groups = bare.ngmesh.GetNrIdentifications()
    info = add_kelvin_identification(bare, skip_if_existing=False)
    exported, oracle = pair_set(mesh), pair_set(bare)
    return {
        "groups": mesh.ngmesh.GetNrIdentifications(),
        "bare_groups_before_pairing": bare_groups,
        "covers_kelvin": has_kelvin_identification(mesh),
        "exported_pairs": len(exported),
        "posthoc_pairs": len(oracle),
        "posthoc_unmatched": info["n_unmatched"],
        "posthoc_max_dist": info["max_dist"],
        "posthoc_offset": tuple(info["kelvin_offset"]),
        "differing_pairs": len(exported ^ oracle),
    }
