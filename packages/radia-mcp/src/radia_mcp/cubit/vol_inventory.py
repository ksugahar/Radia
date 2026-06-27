"""Semantic inventory for Netgen ``.vol`` exports from Cubit/Coreform.

This helper is intentionally broader than
``radia_mcp.radia_ngsolve.netgen_vol``.  The radia-ngsolve parser is the
first-order education path and rejects anything except boundary triangles and
volume tetrahedra.  Cubit, however, is the lab's hex-led and mixed-mesh lane,
so the MCP server also needs a light preflight that can *recognize* hex,
pyramid, wedge, tet, quad, and tri records before routing the file.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable


VOLUME_KIND_BY_NP = {
    4: "tet",
    5: "pyramid",
    6: "wedge",
    8: "hex",
}

SURFACE_KIND_BY_NP = {
    3: "triangle",
    4: "quad",
}


def read_netgen_vol_inventory(path: str | Path) -> dict[str, object]:
    """Read a Netgen ``.vol`` file and return a semantic element inventory."""

    p = Path(path)
    return summarize_netgen_vol_inventory(p.read_text(encoding="utf-8"), source=str(p))


def summarize_netgen_vol_inventory(text: str, source: str | None = None) -> dict[str, object]:
    """Return mixed-element inventory for a Netgen ``.vol`` text.

    The result is a routing preflight, not a solver parser.  It detects whether
    a file belongs to the Netgen tri/tet-only teaching path or the Cubit
    hex/mixed path, and it deliberately refuses to split or reinterpret element
    types.
    """

    lines = text.splitlines()
    surface_rows = _read_counted_section(lines, "surfaceelements", required=False)
    volume_rows = _read_counted_section(lines, "volumeelements", required=False)
    point_rows = _read_counted_section(lines, "points", required=False)
    material_rows = _read_counted_section(lines, "materials", required=False)

    surface_kind_counts = _count_by_np(surface_rows, 4, SURFACE_KIND_BY_NP)
    volume_kind_counts = _count_by_np(volume_rows, 1, VOLUME_KIND_BY_NP)
    materials = _parse_materials(material_rows)
    is_tri_tet_only = (
        set(surface_kind_counts).issubset({"triangle"})
        and set(volume_kind_counts).issubset({"tet"})
        and bool(volume_rows)
    )
    has_mixed_hex_transition = any(
        volume_kind_counts.get(kind, 0) > 0 for kind in ("hex", "pyramid", "wedge")
    )

    routing_hint = (
        "netgen_tri_tet_path"
        if is_tri_tet_only
        else "cubit_hex_or_mixed_path"
        if has_mixed_hex_transition
        else "inspect_before_solver_import"
    )

    return {
        "source": source,
        "surface_elements": len(surface_rows),
        "surface_kind_counts": surface_kind_counts,
        "volume_elements": len(volume_rows),
        "volume_kind_counts": volume_kind_counts,
        "points": len(point_rows),
        "materials": materials,
        "is_tri_tet_only": is_tri_tet_only,
        "has_mixed_hex_transition": has_mixed_hex_transition,
        "routing_hint": routing_hint,
        "policy": (
            "Cubit/Coreform owns hex-led and mixed hex+pyramid+tet inventory; "
            "Netgen/OCC owns tet-only generation for the first-order education path."
        ),
    }


def _next_payload_line(lines: list[str], start: int) -> tuple[int, str]:
    for i in range(start, len(lines)):
        line = lines[i].strip()
        if line and not line.startswith("#"):
            return i, line
    raise ValueError("unexpected end of .vol file")


def _read_counted_section(lines: list[str], name: str, *, required: bool) -> list[str]:
    for i, line in enumerate(lines):
        if line.strip().lower() == name.lower():
            count_i, count_line = _next_payload_line(lines, i + 1)
            count = int(count_line.split()[0])
            rows: list[str] = []
            j = count_i + 1
            while len(rows) < count:
                if j >= len(lines):
                    raise ValueError(f"section {name!r} ended early")
                row = lines[j].strip()
                if row and not row.startswith("#"):
                    rows.append(row)
                j += 1
            return rows
    if required:
        raise ValueError(f"section {name!r} not found")
    return []


def _count_by_np(rows: Iterable[str], np_column: int, kind_by_np: dict[int, str]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in rows:
        parts = row.split()
        np_value = int(parts[np_column])
        counts[kind_by_np.get(np_value, f"np{np_value}")] += 1
    return dict(sorted(counts.items()))


def _parse_materials(rows: Iterable[str]) -> dict[int, str]:
    materials: dict[int, str] = {}
    for row in rows:
        parts = row.split(maxsplit=1)
        if not parts:
            continue
        material_id = int(parts[0])
        materials[material_id] = parts[1] if len(parts) > 1 else f"material_{material_id}"
    return materials
