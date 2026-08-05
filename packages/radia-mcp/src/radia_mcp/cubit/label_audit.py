"""
label_audit.py — pure block/sideset convention audit (no cubit dependency).

Shared between the Cubit daemon (`probe("labels")`, running under Cubit's
bundled Python 3.10) and cubit-free tests / callers in the MCP server's
Python 3.12.  Keep this module free of any import beyond the stdlib and
free of side effects: the daemon imports it by path, pytest imports it as
``radia_mcp.cubit.label_audit``.

Encodes the lab's Block/Sideset Label Convention (CLAUDE.md):
  - a block mixing volume elements and surface elements loses its
    boundary label at Netgen export time (ERROR),
  - unnamed blocks/sidesets export as generated labels (WARNING),
  - casefold name collisions break the label contract (ERROR),
  - names outside check-vol strict naming fail --strict-labels (WARNING).
"""

import re

_STRICT_LABEL_RE = re.compile(
    r"^(?:[a-z0-9]+(?:_[a-z0-9]+)*"
    r"|sym_(?:bn|ht)=0_[xyz]"
    r"|GND)$"
)


def is_strict_label(name):
    """check-vol strict naming: lower snake case, sym_bn=0_<axis> /
    sym_ht=0_<axis>, and the reserved GND point label."""
    return bool(_STRICT_LABEL_RE.match(name))


def audit_label_records(blocks, sidesets):
    """Audit block/sideset records against the lab label convention.

    Args:
        blocks: list of dicts with keys ``id``, ``name``, and the element
            counts ``volume_elems`` / ``surface_elems`` (extra keys such as
            ``volumes`` / ``surfaces`` are carried by the caller, ignored
            here).
        sidesets: list of dicts with keys ``id`` and ``name``.

    Returns:
        {"errors": [...], "warnings": [...], "passed": bool} where
        ``passed`` is False iff any error was found.
    """
    errors, warnings = [], []
    names = []
    for b in blocks:
        label = f'block {b["id"]}'
        name = (b.get("name") or "").strip()
        if b.get("volume_elems") and b.get("surface_elems"):
            errors.append(
                f'{label} ("{name}") mixes volume and surface elements -- '
                "the boundary label is LOST on Netgen export; use a separate "
                "block or a sideset for the boundary"
            )
        if not name:
            warnings.append(f"{label} is unnamed -- exports as a generated label")
        else:
            names.append((label, name))
            if not is_strict_label(name):
                warnings.append(
                    f'{label} name "{name}" is outside strict naming '
                    "(lower snake case / sym_bn=0_<axis> / GND)"
                )
    for s in sidesets:
        label = f'sideset {s["id"]}'
        name = (s.get("name") or "").strip()
        if not name:
            warnings.append(f"{label} is unnamed -- exports as a generated label")
        else:
            names.append((label, name))
            if not is_strict_label(name):
                warnings.append(
                    f'{label} name "{name}" is outside strict naming '
                    "(lower snake case / sym_bn=0_<axis> / GND)"
                )
    seen = {}
    for label, name in names:
        key = name.casefold()
        if key in seen and seen[key][1] != name:
            errors.append(
                f'casefold collision: {seen[key][0]} "{seen[key][1]}" vs '
                f'{label} "{name}"'
            )
        seen.setdefault(key, (label, name))
    return {"errors": errors, "warnings": warnings, "passed": not errors}
