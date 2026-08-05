"""Pure-Python GMSH MSH v4.1 inspection and validation.

The structural layer (``inspect_msh`` / ``validate_msh`` / ``validate_geo``)
has no gmsh dependency, so it works on any machine including minimal-dep CI.
The optional high-order Jacobian check runs the gmsh Python API in a
subprocess -- never in-process -- so a crashing gmsh cannot take the MCP
server down and the stdio transport cannot dead-lock on inherited pipes
(temp-file redirect + ``stdin=DEVNULL``, per the lab subprocess pattern).

Scope: MSH **v4.1 ASCII** only, matching the repository-wide format policy.
Binary MSH and v2.2 are rejected loudly, never half-parsed.
"""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from importlib import util as _importlib_util
from pathlib import Path
from typing import Any, Iterator

# MSH element type registry: code -> (name, nodes_per_element, dim, order).
# Covers every type emitted by cubit-mesh-export (order 1-3 gmsh export),
# GmshPostExport (Tri3..Tri21), and the standard linear/quadratic families.
ELEMENT_TYPES: dict[int, tuple[str, int, int, int]] = {
    1: ("line2", 2, 1, 1),
    2: ("tri3", 3, 2, 1),
    3: ("quad4", 4, 2, 1),
    4: ("tet4", 4, 3, 1),
    5: ("hex8", 8, 3, 1),
    6: ("prism6", 6, 3, 1),
    7: ("pyr5", 5, 3, 1),
    8: ("line3", 3, 1, 2),
    9: ("tri6", 6, 2, 2),
    10: ("quad9", 9, 2, 2),
    11: ("tet10", 10, 3, 2),
    12: ("hex27", 27, 3, 2),
    13: ("prism18", 18, 3, 2),
    14: ("pyr14", 14, 3, 2),
    15: ("point1", 1, 0, 1),
    16: ("quad8", 8, 2, 2),
    17: ("hex20", 20, 3, 2),
    18: ("prism15", 15, 3, 2),
    19: ("pyr13", 13, 3, 2),
    20: ("tri9", 9, 2, 3),
    21: ("tri10", 10, 2, 3),
    22: ("tri12", 12, 2, 4),
    23: ("tri15", 15, 2, 4),
    24: ("tri15i", 15, 2, 5),
    25: ("tri21", 21, 2, 5),
    26: ("line4", 4, 1, 3),
    27: ("line5", 5, 1, 4),
    28: ("line6", 6, 1, 5),
    29: ("tet20", 20, 3, 3),
    30: ("tet35", 35, 3, 4),
    31: ("tet56", 56, 3, 5),
    36: ("quad16", 16, 2, 3),
    37: ("quad25", 25, 2, 4),
    38: ("quad36", 36, 2, 5),
    92: ("hex64", 64, 3, 3),
    93: ("hex125", 125, 3, 4),
}

_DATA_SECTIONS = ("NodeData", "ElementData", "ElementNodeData")
_MAX_DETAIL = 20  # cap per-issue detail strings so huge files stay readable

# Known-bad options for .geo validation.  These do not exist in GMSH 4.x and
# crash or error on merge; each maps to the correct replacement.
INVALID_GEO_OPTIONS: dict[str, str] = {
    "Mesh.Volumes": "use Mesh.VolumeEdges + Mesh.VolumeFaces",
    "Mesh.Surfaces": "use Mesh.SurfaceEdges + Mesh.SurfaceFaces",
    "General.GraphicsSizeX": "use General.GraphicsWidth",
    "General.GraphicsSizeY": "use General.GraphicsHeight",
}


# ======================================================================
# Section splitting
# ======================================================================

def _split_sections(lines: list[str]) -> tuple[list[tuple[str, list[str], int]], list[str]]:
    """Split raw lines into (name, body_lines, start_line) tuples.

    Returns (sections, errors) where errors report unbalanced
    ``$Section`` / ``$EndSection`` pairs.
    """
    sections: list[tuple[str, list[str], int]] = []
    errors: list[str] = []
    open_name: str | None = None
    open_start = 0
    body: list[str] = []

    for i, raw in enumerate(lines, 1):
        s = raw.strip()
        if s.startswith("$End"):
            end_name = s[4:]
            if open_name is None:
                errors.append(f"L{i}: {s} without an opening section")
            elif end_name != open_name:
                errors.append(f"L{i}: {s} closes ${open_name}")
                sections.append((open_name, body, open_start))
                open_name, body = None, []
            else:
                sections.append((open_name, body, open_start))
                open_name, body = None, []
        elif s.startswith("$"):
            if open_name is not None:
                errors.append(f"L{i}: ${open_name} not closed before {s}")
                sections.append((open_name, body, open_start))
            open_name = s[1:]
            open_start = i
            body = []
        elif open_name is not None:
            body.append(raw)
        # tokens outside any section: blank lines only in well-formed files

    if open_name is not None:
        errors.append(f"${open_name} never closed ($End{open_name} missing)")
        sections.append((open_name, body, open_start))
    return sections, errors


def _token_stream(lines: list[str]) -> Iterator[str]:
    for ln in lines:
        yield from ln.split()


def _cap(items: list[str]) -> list[str]:
    if len(items) > _MAX_DETAIL:
        return items[:_MAX_DETAIL] + [f"... and {len(items) - _MAX_DETAIL} more"]
    return items


# ======================================================================
# Per-section parsers
# ======================================================================

def _parse_physical_names(body: list[str]) -> tuple[list[dict[str, Any]], list[str]]:
    names: list[dict[str, Any]] = []
    errors: list[str] = []
    rows = [ln.strip() for ln in body if ln.strip()]
    if not rows:
        return names, ["$PhysicalNames: empty section"]
    try:
        declared = int(rows[0])
    except ValueError:
        return names, [f"$PhysicalNames: bad count line {rows[0]!r}"]
    for row in rows[1:]:
        m = re.match(r'(\d+)\s+(\d+)\s+"(.*)"\s*$', row)
        if not m:
            errors.append(f"$PhysicalNames: malformed row {row!r}")
            continue
        names.append({"dim": int(m.group(1)), "tag": int(m.group(2)),
                      "name": m.group(3)})
    if declared != len(names):
        errors.append(
            f"$PhysicalNames: declared {declared} but parsed {len(names)}")
    return names, errors


def _parse_entities(body: list[str]) -> tuple[dict[str, Any], list[str]]:
    """Parse $Entities enough to count entities and collect physical tags."""
    toks = _token_stream(body)
    errors: list[str] = []
    result: dict[str, Any] = {
        "counts": None,
        "physical_tags": {0: {}, 1: {}, 2: {}, 3: {}},  # dim -> {tag: [phys]}
    }
    try:
        counts = [int(next(toks)) for _ in range(4)]
        result["counts"] = {"points": counts[0], "curves": counts[1],
                            "surfaces": counts[2], "volumes": counts[3]}
        # points: tag x y z numPhys phys...
        for _ in range(counts[0]):
            tag = int(next(toks))
            for _ in range(3):
                next(toks)
            n_phys = int(next(toks))
            result["physical_tags"][0][tag] = [int(next(toks)) for _ in range(n_phys)]
        # curves / surfaces / volumes: tag bbox(6) numPhys phys... numBnd bnd...
        for dim in (1, 2, 3):
            for _ in range(counts[dim]):
                tag = int(next(toks))
                for _ in range(6):
                    next(toks)
                n_phys = int(next(toks))
                result["physical_tags"][dim][tag] = [int(next(toks)) for _ in range(n_phys)]
                n_bnd = int(next(toks))
                for _ in range(n_bnd):
                    next(toks)
    except (StopIteration, ValueError) as exc:
        errors.append(f"$Entities: truncated or malformed ({exc})")
    return result, errors


def _parse_nodes(body: list[str]) -> tuple[dict[str, Any], list[str]]:
    toks = _token_stream(body)
    errors: list[str] = []
    out: dict[str, Any] = {
        "declared": None, "blocks": None, "min_tag": None, "max_tag": None,
        "tags": [], "bbox_min": None, "bbox_max": None,
    }
    try:
        n_blocks = int(next(toks))
        out["blocks"] = n_blocks
        out["declared"] = int(next(toks))
        out["min_tag"] = int(next(toks))
        out["max_tag"] = int(next(toks))
        lo = [math.inf] * 3
        hi = [-math.inf] * 3
        for _ in range(n_blocks):
            dim = int(next(toks))
            next(toks)  # entity tag
            parametric = int(next(toks))
            n_in_block = int(next(toks))
            block_tags = [int(next(toks)) for _ in range(n_in_block)]
            out["tags"].extend(block_tags)
            n_extra = dim if parametric else 0
            for _ in range(n_in_block):
                x = float(next(toks))
                y = float(next(toks))
                z = float(next(toks))
                for _ in range(n_extra):
                    next(toks)
                for k, v in enumerate((x, y, z)):
                    if v < lo[k]:
                        lo[k] = v
                    if v > hi[k]:
                        hi[k] = v
        if out["tags"]:
            out["bbox_min"] = lo
            out["bbox_max"] = hi
    except (StopIteration, ValueError) as exc:
        errors.append(f"$Nodes: truncated or malformed ({exc})")
    return out, errors


def _parse_elements(body: list[str]) -> tuple[dict[str, Any], list[str]]:
    toks = _token_stream(body)
    errors: list[str] = []
    out: dict[str, Any] = {
        "declared": None, "blocks": None, "min_tag": None, "max_tag": None,
        "tags": [], "by_type": Counter(), "connectivity": {},  # type -> list[list[int]] capped
        "node_refs": set(), "unknown_types": [], "parse_complete": True,
    }
    try:
        n_blocks = int(next(toks))
        out["blocks"] = n_blocks
        out["declared"] = int(next(toks))
        out["min_tag"] = int(next(toks))
        out["max_tag"] = int(next(toks))
        for _ in range(n_blocks):
            next(toks)  # entity dim
            next(toks)  # entity tag
            etype = int(next(toks))
            n_in_block = int(next(toks))
            if etype not in ELEMENT_TYPES:
                # Token width per element is unknown for unrecognized types:
                # stop loudly instead of guessing (fail fast, no silent skip).
                out["unknown_types"].append(etype)
                out["parse_complete"] = False
                errors.append(
                    f"$Elements: unknown element type {etype}; remaining "
                    f"element parsing stopped (extend ELEMENT_TYPES if this "
                    f"type is legitimate)")
                break
            n_per = ELEMENT_TYPES[etype][1]
            out["by_type"][etype] += n_in_block
            for _ in range(n_in_block):
                out["tags"].append(int(next(toks)))
                refs = [int(next(toks)) for _ in range(n_per)]
                out["node_refs"].update(refs)
    except (StopIteration, ValueError) as exc:
        out["parse_complete"] = False
        errors.append(f"$Elements: truncated or malformed ({exc})")
    return out, errors


def _parse_data_section(kind: str, body: list[str], start_line: int) -> tuple[dict[str, Any], list[str]]:
    errors: list[str] = []
    view: dict[str, Any] = {
        "section": kind, "start_line": start_line, "name": "",
        "time": None, "step": None, "components": None,
        "declared": None, "data_rows": 0, "tags": [],
        "header_ok": True,
    }
    rows = [ln.strip() for ln in body if ln.strip()]
    idx = 0

    def _take() -> str:
        nonlocal idx
        if idx >= len(rows):
            raise StopIteration
        val = rows[idx]
        idx += 1
        return val

    try:
        n_str = int(_take())
        strings = [_take().strip('"') for _ in range(n_str)]
        n_real = int(_take())
        reals = [float(_take()) for _ in range(n_real)]
        n_int = int(_take())
        ints = [int(_take()) for _ in range(n_int)]
    except (StopIteration, ValueError):
        view["header_ok"] = False
        errors.append(f"${kind} @L{start_line}: truncated or malformed header")
        return view, errors

    view["name"] = strings[0] if strings else ""
    view["time"] = reals[0] if reals else None
    if len(ints) >= 3:
        view["step"], view["components"], view["declared"] = ints[0], ints[1], ints[2]
    else:
        view["header_ok"] = False
        errors.append(
            f"${kind} \"{view['name']}\" @L{start_line}: integer-tag block "
            f"needs >=3 entries (step, components, count), got {len(ints)}")

    for row in rows[idx:]:
        view["data_rows"] += 1
        tok = row.split(None, 1)[0]
        try:
            view["tags"].append(int(tok))
        except ValueError:
            errors.append(
                f"${kind} \"{view['name']}\" @L{start_line}: non-integer "
                f"tag {tok!r} in data row")
    return view, errors


# ======================================================================
# Whole-file parse
# ======================================================================

def _parse_msh(path: Path) -> dict[str, Any]:
    """Parse an MSH file into an internal structure (ASCII v4.x only)."""
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    sections, section_errors = _split_sections(lines)

    parsed: dict[str, Any] = {
        "version": None, "ascii": None,
        "section_names": [name for name, _, _ in sections],
        "section_errors": section_errors,
        "physical_names": None, "entities": None,
        "nodes": None, "elements": None, "views_raw": [],
        "errors": [], "warnings": [],
    }

    fmt = next(((b, ln) for name, b, ln in sections if name == "MeshFormat"), None)
    if fmt is None:
        parsed["errors"].append("$MeshFormat section missing")
        return parsed
    fmt_tokens = " ".join(fmt[0]).split()
    if len(fmt_tokens) < 2:
        parsed["errors"].append("$MeshFormat: malformed header")
        return parsed
    parsed["version"] = fmt_tokens[0]
    parsed["ascii"] = fmt_tokens[1] == "0"
    if not parsed["version"].startswith("4"):
        parsed["errors"].append(
            f"MSH v{parsed['version']} is not supported; the lab standard is "
            f"v4.1 (re-export, do not downgrade the checker)")
        return parsed
    if not parsed["ascii"]:
        parsed["errors"].append(
            "binary MSH is not supported by this validator; re-export ASCII")
        return parsed

    n_nodes_sections = 0
    n_elems_sections = 0
    for name, body, start_line in sections:
        if name == "PhysicalNames":
            parsed["physical_names"], errs = _parse_physical_names(body)
            parsed["errors"].extend(errs)
        elif name == "Entities":
            parsed["entities"], errs = _parse_entities(body)
            parsed["errors"].extend(errs)
        elif name == "Nodes":
            n_nodes_sections += 1
            if n_nodes_sections > 1:
                parsed["warnings"].append("multiple $Nodes sections; using the first")
                continue
            parsed["nodes"], errs = _parse_nodes(body)
            parsed["errors"].extend(errs)
        elif name == "Elements":
            n_elems_sections += 1
            if n_elems_sections > 1:
                parsed["warnings"].append("multiple $Elements sections; using the first")
                continue
            parsed["elements"], errs = _parse_elements(body)
            parsed["errors"].extend(errs)
        elif name in _DATA_SECTIONS:
            view, errs = _parse_data_section(name, body, start_line)
            parsed["views_raw"].append(view)
            parsed["errors"].extend(errs)
    return parsed


def _group_views(views_raw: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group per-section data blocks into named views with time steps."""
    grouped: dict[tuple[str, str], dict[str, Any]] = {}
    for view in views_raw:
        key = (view["section"], view["name"])
        g = grouped.setdefault(key, {
            "section": view["section"], "name": view["name"],
            "components": view["components"], "steps": 0,
            "times": [], "entries_per_step": [],
            "components_consistent": True,
        })
        g["steps"] += 1
        if view["time"] is not None:
            g["times"].append(view["time"])
        g["entries_per_step"].append(view["data_rows"])
        if view["components"] != g["components"]:
            g["components_consistent"] = False
    out = []
    for g in grouped.values():
        times = g.pop("times")
        g["time_range"] = [min(times), max(times)] if times else None
        entries = g.pop("entries_per_step")
        g["entries_per_step"] = entries[0] if len(set(entries)) == 1 else entries
        out.append(g)
    return out


# ======================================================================
# Public API: inspect
# ======================================================================

def inspect_msh(msh_path: str | Path) -> dict[str, Any]:
    """Summarize the structure of an MSH v4.1 file (no gmsh dependency)."""
    path = Path(msh_path)
    if not path.is_file():
        return {"ok": False, "path": str(path), "error": f"file not found: {path}"}

    parsed = _parse_msh(path)
    result: dict[str, Any] = {
        "ok": True,
        "path": str(path),
        "file_size_bytes": path.stat().st_size,
        "version": parsed["version"],
        "ascii": parsed["ascii"],
        "sections": parsed["section_names"],
        "parse_errors": _cap(parsed["section_errors"] + parsed["errors"]),
        "warnings": parsed["warnings"],
    }
    if parsed["version"] is None or not str(parsed["version"]).startswith("4") \
            or parsed["ascii"] is False:
        result["ok"] = False
        return result

    result["physical_names"] = parsed["physical_names"]
    if parsed["entities"] is not None:
        result["entities"] = parsed["entities"]["counts"]

    nodes = parsed["nodes"]
    if nodes is not None:
        result["nodes"] = {
            "count": len(nodes["tags"]), "declared": nodes["declared"],
            "blocks": nodes["blocks"], "min_tag": nodes["min_tag"],
            "max_tag": nodes["max_tag"],
        }
        if nodes["bbox_min"] is not None:
            result["bbox"] = {"min": nodes["bbox_min"], "max": nodes["bbox_max"]}

    elements = parsed["elements"]
    high_order = False
    max_dim = None
    if elements is not None:
        by_type = []
        for etype, count in sorted(elements["by_type"].items()):
            name, n_per, dim, order = ELEMENT_TYPES[etype]
            by_type.append({
                "type": etype, "name": name, "dim": dim, "order": order,
                "nodes_per_element": n_per, "count": count,
            })
            high_order = high_order or order >= 2
            max_dim = dim if max_dim is None else max(max_dim, dim)
        result["elements"] = {
            "count": len(elements["tags"]), "declared": elements["declared"],
            "blocks": elements["blocks"], "by_type": by_type,
            "max_dim": max_dim, "parse_complete": elements["parse_complete"],
        }
        if elements["unknown_types"]:
            result["elements"]["unknown_types"] = elements["unknown_types"]

    result["views"] = _group_views(parsed["views_raw"])
    result["high_order"] = high_order

    hints = []
    if high_order:
        hints.append(
            "order>=2 elements present: open with Mesh.NumSubEdges = 4 "
            "(companion .geo or `gmsh file.msh -numsubedges 4`)")
    if result["views"]:
        multi_step = [v for v in result["views"] if v["steps"] > 1]
        if multi_step:
            hints.append(
                f"{len(multi_step)} time-stepped view(s): use PostProcessing."
                f"Link = 1 and per-view TimeStep for synchronized animation")
    if not parsed["physical_names"]:
        hints.append(
            "no $PhysicalNames: GMSH shows elements without group names; "
            "Radia/Cubit exporters normally emit physical names")
    result["hints"] = hints
    return result


# ======================================================================
# Public API: validate .msh
# ======================================================================

def validate_msh(msh_path: str | Path,
                 check_jacobians: bool = False,
                 quadrature: str = "Gauss2",
                 timeout_s: float = 300.0) -> dict[str, Any]:
    """Validate MSH v4.1 structural consistency; optional Jacobian check.

    The Jacobian check implements the repository policy that high-order
    meshes are verified through the GMSH API (``getJacobians``): negative
    determinants indicate inverted elements from wrong node ordering.
    """
    path = Path(msh_path)
    if not path.is_file():
        return {"ok": False, "path": str(path), "status": "needs_attention",
                "checks": {}, "errors": [f"file not found: {path}"],
                "warnings": []}

    parsed = _parse_msh(path)
    errors: list[str] = list(parsed["section_errors"]) + list(parsed["errors"])
    warnings: list[str] = list(parsed["warnings"])
    checks: dict[str, bool] = {
        "format_is_v41": str(parsed["version"]) == "4.1",
        "is_ascii": bool(parsed["ascii"]),
        "sections_balanced": not parsed["section_errors"],
    }

    if not checks["format_is_v41"] or not checks["is_ascii"]:
        status = "needs_attention"
        return {"ok": False, "path": str(path), "status": status,
                "checks": checks, "errors": _cap(errors), "warnings": warnings}

    nodes = parsed["nodes"]
    elements = parsed["elements"]
    checks["nodes_section_present"] = nodes is not None
    checks["elements_section_present"] = elements is not None

    node_tag_set: set[int] = set()
    if nodes is not None:
        node_tag_set = set(nodes["tags"])
        checks["node_header_count_matches"] = len(nodes["tags"]) == nodes["declared"]
        if not checks["node_header_count_matches"]:
            errors.append(
                f"$Nodes header declares {nodes['declared']} nodes but "
                f"{len(nodes['tags'])} were parsed")
        checks["node_tags_unique"] = len(node_tag_set) == len(nodes["tags"])
        if not checks["node_tags_unique"]:
            dupes = [t for t, c in Counter(nodes["tags"]).items() if c > 1]
            errors.append(f"duplicate node tags: {_cap([str(d) for d in dupes])}")

    if elements is not None:
        checks["element_parse_complete"] = elements["parse_complete"]
        checks["element_types_known"] = not elements["unknown_types"]
        if elements["parse_complete"]:
            checks["element_header_count_matches"] = (
                len(elements["tags"]) == elements["declared"])
            if not checks["element_header_count_matches"]:
                errors.append(
                    f"$Elements header declares {elements['declared']} elements "
                    f"but {len(elements['tags'])} were parsed")
        tag_counter = Counter(elements["tags"])
        dupes = [t for t, c in tag_counter.items() if c > 1]
        checks["element_tags_unique"] = not dupes
        if dupes:
            errors.append(
                f"duplicate element tags (v4.1 requires globally unique "
                f"element tags): {_cap([str(d) for d in dupes])}")
        if nodes is not None:
            missing = elements["node_refs"] - node_tag_set
            checks["element_node_refs_exist"] = not missing
            if missing:
                errors.append(
                    f"elements reference {len(missing)} undefined node tag(s): "
                    f"{_cap([str(m) for m in sorted(missing)])}")

    # Entity physical tags must be declared when $PhysicalNames exists.
    if parsed["physical_names"] and parsed["entities"] is not None:
        declared_tags = {(p["dim"], p["tag"]) for p in parsed["physical_names"]}
        undeclared = []
        for dim, tag_map in parsed["entities"]["physical_tags"].items():
            for etag, phys_list in tag_map.items():
                for phys in phys_list:
                    if (dim, abs(phys)) not in declared_tags:
                        undeclared.append(f"dim{dim} entity {etag} -> physical {phys}")
        checks["entity_physical_tags_declared"] = not undeclared
        if undeclared:
            errors.append(
                "entities reference physical tags missing from "
                f"$PhysicalNames: {_cap(undeclared)}")

    views_raw = parsed["views_raw"]
    if views_raw:
        checks["data_headers_wellformed"] = all(v["header_ok"] for v in views_raw)
        count_ok = True
        tags_ok = True
        comp_ok = True
        for view in views_raw:
            label = f"${view['section']} \"{view['name']}\" @L{view['start_line']}"
            if view["declared"] is not None and view["declared"] != view["data_rows"]:
                count_ok = False
                errors.append(
                    f"{label}: declares {view['declared']} entries but has "
                    f"{view['data_rows']} data rows")
            if view["components"] is not None and view["components"] not in (1, 3, 9):
                comp_ok = False
                errors.append(
                    f"{label}: numComponents={view['components']} "
                    f"(must be 1, 3, or 9)")
            if view["section"] == "NodeData" and node_tag_set:
                missing = [t for t in view["tags"] if t not in node_tag_set]
                if missing:
                    tags_ok = False
                    errors.append(
                        f"{label}: {len(missing)} node tag(s) not in $Nodes: "
                        f"{_cap([str(m) for m in missing])}")
            elif view["section"] in ("ElementData", "ElementNodeData") \
                    and elements is not None and elements["parse_complete"]:
                elem_set = set(elements["tags"])
                missing = [t for t in view["tags"] if t not in elem_set]
                if missing:
                    tags_ok = False
                    errors.append(
                        f"{label}: {len(missing)} element tag(s) not in "
                        f"$Elements: {_cap([str(m) for m in missing])}")
        checks["data_declared_counts_match"] = count_ok
        checks["data_components_valid"] = comp_ok
        checks["data_tags_exist"] = tags_ok
        grouped = _group_views(views_raw)
        checks["view_components_consistent"] = all(
            g["components_consistent"] for g in grouped)
        for g in grouped:
            if not g["components_consistent"]:
                errors.append(
                    f"view \"{g['name']}\": numComponents changes between "
                    f"time steps")

    result: dict[str, Any] = {
        "path": str(path),
        "checks": checks,
        "errors": _cap(errors),
        "warnings": warnings,
    }

    if check_jacobians:
        jac = _run_jacobian_check(path, quadrature, timeout_s)
        result["jacobian"] = jac
        checks["jacobians_positive"] = bool(jac.get("ok")) and bool(jac.get("ran"))
        if not jac.get("ran"):
            errors.append(f"jacobian check did not run: {jac.get('error', 'unknown')}")
            result["errors"] = _cap(errors)

    result["status"] = "ok" if all(checks.values()) else "needs_attention"
    result["ok"] = result["status"] == "ok"
    return result


# ======================================================================
# Jacobian subprocess (gmsh API, per repo verification policy)
# ======================================================================

_JACOBIAN_CHECK_SCRIPT = r"""
import json
import sys

msh_path, out_path, quadrature = sys.argv[1], sys.argv[2], sys.argv[3]
result = {"ok": False, "ran": False}
try:
    import numpy as np
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(msh_path)
        by_type = []
        total_neg = 0
        total_volume = 0.0
        etypes, _, _ = gmsh.model.mesh.getElements(3)
        for et in etypes:
            name, _dim, order, _nn, _coords, _ = \
                gmsh.model.mesh.getElementProperties(int(et))
            local, weights = gmsh.model.mesh.getIntegrationPoints(
                int(et), quadrature)
            _jac, det, _pts = gmsh.model.mesh.getJacobians(int(et), local)
            det = np.asarray(det, dtype=float)
            w = np.asarray(weights, dtype=float)
            n_q = len(w)
            n_el = det.size // n_q if n_q else 0
            neg = int((det <= 0.0).sum())
            vol = float((det.reshape(n_el, n_q) * w).sum()) if n_el else 0.0
            by_type.append({
                "type": int(et), "name": str(name), "order": int(order),
                "n_elements": int(n_el),
                "min_det": float(det.min()) if det.size else None,
                "negative_count": neg,
                "volume": vol,
            })
            total_neg += neg
            total_volume += vol
        result.update({
            "ran": True,
            "ok": total_neg == 0,
            "applicable": bool(by_type),
            "quadrature": quadrature,
            "by_type": by_type,
            "total_negative": total_neg,
            "total_volume_dim3": total_volume,
        })
        if not by_type:
            result["ok"] = True
            result["note"] = ("no 3D elements; Jacobian sign check not "
                              "applicable (surface dets are norms, always >0)")
    finally:
        gmsh.finalize()
except Exception as exc:
    result["error"] = f"{type(exc).__name__}: {exc}"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def _run_jacobian_check(msh_path: Path, quadrature: str,
                        timeout_s: float) -> dict[str, Any]:
    """Run the gmsh getJacobians check in a subprocess.

    Never imports gmsh in-process: gmsh keeps global state and a hard
    crash on a corrupt file must not kill the MCP server.  Output goes
    through temp files with stdin=DEVNULL (stdio-MCP pipe-hang pattern).
    """
    if _importlib_util.find_spec("gmsh") is None:
        return {"ok": False, "ran": False,
                "error": "gmsh Python package not installed (pip install gmsh)"}

    work = tempfile.mkdtemp(prefix="radia_mcp_gmsh_jac_")
    out_json = str(Path(work) / "result.json")
    log_path = Path(work) / "run.log"
    cmd = [sys.executable, "-c", _JACOBIAN_CHECK_SCRIPT,
           str(msh_path), out_json, quadrature]
    try:
        with open(log_path, "w", encoding="utf-8") as log:
            proc = subprocess.run(
                cmd, stdin=subprocess.DEVNULL, stdout=log,
                stderr=subprocess.STDOUT, timeout=timeout_s, check=False)
        if not Path(out_json).is_file():
            tail = log_path.read_text(encoding="utf-8", errors="replace")[-1500:]
            return {"ok": False, "ran": False,
                    "error": f"gmsh subprocess exited {proc.returncode} "
                             f"without result: {tail}"}
        return json.loads(Path(out_json).read_text(encoding="utf-8"))
    except subprocess.TimeoutExpired:
        return {"ok": False, "ran": False,
                "error": f"gmsh subprocess timed out after {timeout_s}s"}
    finally:
        shutil.rmtree(work, ignore_errors=True)


# ======================================================================
# Public API: validate .geo
# ======================================================================

def validate_geo(geo_path: str | Path) -> dict[str, Any]:
    """Validate a .geo launch/companion file before opening it in GMSH.

    Checks that every ``Merge "..."`` target exists on disk and that no
    known-invalid GMSH 4.x option names are used.  Also reports the max
    ``View[N]`` index and whether ``Mesh.NumSubEdges`` is set.
    """
    path = Path(geo_path)
    if not path.is_file():
        return {"ok": False, "path": str(path), "status": "needs_attention",
                "checks": {}, "errors": [f"file not found: {path}"],
                "warnings": []}

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    errors: list[str] = []
    warnings: list[str] = []
    merge_targets: list[dict[str, Any]] = []
    invalid_options: list[dict[str, Any]] = []
    max_view_index = -1
    numsubedges: int | None = None

    for lineno, raw in enumerate(lines, 1):
        stripped = raw.split("//")[0].strip()
        if not stripped:
            continue

        m = re.match(r'Merge\s+"([^"]+)"', stripped)
        if m:
            target = m.group(1)
            resolved = (path.parent / target) if not Path(target).is_absolute() \
                else Path(target)
            exists = resolved.is_file()
            entry: dict[str, Any] = {"line": lineno, "target": target,
                                     "exists": exists,
                                     "resolved": str(resolved)}
            if exists:
                entry["size_bytes"] = resolved.stat().st_size
            else:
                errors.append(f'L{lineno}: Merge target missing: "{target}"')
            merge_targets.append(entry)

        for opt, fix in INVALID_GEO_OPTIONS.items():
            if re.match(rf"{re.escape(opt)}\s*=", stripped):
                invalid_options.append({"line": lineno, "option": opt, "fix": fix})
                errors.append(
                    f"L{lineno}: invalid GMSH 4.x option {opt} -- {fix}")

        for idx_str in re.findall(r"View\[(\d+)\]", stripped):
            max_view_index = max(max_view_index, int(idx_str))

        m = re.match(r"Mesh\.NumSubEdges\s*=\s*(\d+)", stripped)
        if m:
            numsubedges = int(m.group(1))

    checks = {
        "merge_targets_exist": all(t["exists"] for t in merge_targets),
        "no_invalid_options": not invalid_options,
    }
    if not merge_targets:
        warnings.append(
            "no Merge lines found: this .geo does not load any data file")

    status = "ok" if all(checks.values()) else "needs_attention"
    return {
        "ok": status == "ok",
        "path": str(path),
        "status": status,
        "checks": checks,
        "merge_targets": merge_targets,
        "invalid_options": invalid_options,
        "max_view_index": max_view_index,
        "numsubedges": numsubedges,
        "errors": _cap(errors),
        "warnings": warnings,
    }
