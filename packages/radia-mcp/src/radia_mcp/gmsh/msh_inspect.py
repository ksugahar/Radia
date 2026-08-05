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

import math
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterator

from ._gmsh_subprocess import run_gmsh_json_subprocess

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


def _parse_data_section(kind: str, body: list[str], start_line: int,
                        parse_values: bool = False) -> tuple[dict[str, Any], list[str]]:
    errors: list[str] = []
    view: dict[str, Any] = {
        "section": kind, "start_line": start_line, "name": "",
        "time": None, "step": None, "components": None,
        "declared": None, "data_rows": 0, "tags": [],
        "header_ok": True, "value_stats": None,
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

    ncomp = view["components"]
    stats: dict[str, Any] | None = None
    if parse_values:
        stats = {"nan": 0, "inf": 0, "bad_width_rows": 0,
                 "comp_min": None, "comp_max": None,
                 "min": None, "max": None, "s1": 0.0, "s2": 0.0,
                 "n_finite_rows": 0}

    for row in rows[idx:]:
        view["data_rows"] += 1
        toks = row.split()
        try:
            view["tags"].append(int(toks[0]))
        except (ValueError, IndexError):
            errors.append(
                f"${kind} \"{view['name']}\" @L{start_line}: non-integer "
                f"tag {toks[0] if toks else row!r} in data row")
            continue
        if stats is None:
            continue

        # ElementNodeData rows carry an extra per-element node count.
        vals_start = 2 if kind == "ElementNodeData" else 1
        try:
            vals = [float(t) for t in toks[vals_start:]]
        except ValueError:
            stats["bad_width_rows"] += 1
            continue
        if ncomp:
            expected = ncomp
            if kind == "ElementNodeData" and len(toks) > 1:
                try:
                    expected = ncomp * int(toks[1])
                except ValueError:
                    expected = None
            if expected is not None and len(vals) != expected:
                stats["bad_width_rows"] += 1
        stats["nan"] += sum(1 for v in vals if math.isnan(v))
        stats["inf"] += sum(1 for v in vals if math.isinf(v))
        finite = [v for v in vals if math.isfinite(v)]
        if not finite:
            continue
        lo, hi = min(finite), max(finite)
        stats["comp_min"] = lo if stats["comp_min"] is None else min(stats["comp_min"], lo)
        stats["comp_max"] = hi if stats["comp_max"] is None else max(stats["comp_max"], hi)
        if len(vals) == len(finite):
            # scalar: signed value; vector/tensor: euclidean magnitude
            metric = vals[0] if ncomp == 1 and len(vals) == 1 else \
                math.sqrt(sum(v * v for v in vals))
            stats["min"] = metric if stats["min"] is None else min(stats["min"], metric)
            stats["max"] = metric if stats["max"] is None else max(stats["max"], metric)
            stats["s1"] += metric
            stats["s2"] += metric * metric
            stats["n_finite_rows"] += 1

    view["value_stats"] = stats
    return view, errors


# ======================================================================
# Whole-file parse
# ======================================================================

def _parse_msh(path: Path, parse_values: bool = False) -> dict[str, Any]:
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
            view, errs = _parse_data_section(name, body, start_line,
                                             parse_values=parse_values)
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

    parsed = _parse_msh(path, parse_values=True)
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
        finite_ok = True
        width_ok = True
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
            st = view.get("value_stats")
            if st is not None:
                if st["nan"] or st["inf"]:
                    finite_ok = False
                    errors.append(
                        f"{label}: {st['nan']} NaN / {st['inf']} Inf data "
                        f"value(s) -- GMSH renders these silently wrong")
                if st["bad_width_rows"]:
                    width_ok = False
                    errors.append(
                        f"{label}: {st['bad_width_rows']} data row(s) whose "
                        f"value count does not match numComponents")
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
        checks["data_values_finite"] = finite_ok
        checks["data_row_width_matches"] = width_ok
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

msh_path, quadrature, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
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
    crash on a corrupt file must not kill the MCP server.
    """
    return run_gmsh_json_subprocess(
        _JACOBIAN_CHECK_SCRIPT, [str(msh_path), quadrature],
        timeout_s=timeout_s, prefix="radia_mcp_gmsh_jac_")


# ======================================================================
# Public API: field statistics
# ======================================================================

def field_stats(msh_path: str | Path,
                view_name: str | None = None) -> dict[str, Any]:
    """Per-view, per-time-step field value statistics for an MSH file.

    Scalars report signed min/max/mean/rms; vectors and tensors report
    Euclidean-magnitude statistics plus the pooled per-component
    min/max.  NaN/Inf values are counted (and excluded from the stats).
    """
    path = Path(msh_path)
    if not path.is_file():
        return {"ok": False, "path": str(path),
                "error": f"file not found: {path}"}

    parsed = _parse_msh(path, parse_values=True)
    if parsed["version"] is None or not str(parsed["version"]).startswith("4") \
            or parsed["ascii"] is False:
        return {"ok": False, "path": str(path),
                "error": "not an ASCII MSH v4.x file",
                "parse_errors": _cap(parsed["errors"])}

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for view in parsed["views_raw"]:
        grouped.setdefault((view["section"], view["name"]), []).append(view)

    if view_name is not None:
        names = sorted({name for _, name in grouped})
        grouped = {key: v for key, v in grouped.items() if key[1] == view_name}
        if not grouped:
            return {"ok": False, "path": str(path),
                    "error": f"view {view_name!r} not found; "
                             f"available: {names}"}

    views_out = []
    for (section, name), sections in grouped.items():
        sections.sort(key=lambda v: (v["step"] if v["step"] is not None
                                     else v["start_line"]))
        per_step = []
        overall_nan = overall_inf = 0
        overall_min = overall_max = None
        for view in sections:
            st = view.get("value_stats") or {}
            n = st.get("n_finite_rows", 0)
            entry: dict[str, Any] = {
                "step": view["step"], "time": view["time"],
                "entities": view["data_rows"],
                "nan": st.get("nan", 0), "inf": st.get("inf", 0),
                "min": st.get("min"), "max": st.get("max"),
                "mean": (st.get("s1", 0.0) / n) if n else None,
                "rms": math.sqrt(st.get("s2", 0.0) / n) if n else None,
            }
            if (view["components"] or 1) > 1:
                entry["comp_min"] = st.get("comp_min")
                entry["comp_max"] = st.get("comp_max")
                entry["metric"] = "magnitude"
            else:
                entry["metric"] = "value"
            per_step.append(entry)
            overall_nan += entry["nan"]
            overall_inf += entry["inf"]
            for bound, pick in (("min", min), ("max", max)):
                val = entry[bound]
                if val is None:
                    continue
                current = overall_min if bound == "min" else overall_max
                new = val if current is None else pick(current, val)
                if bound == "min":
                    overall_min = new
                else:
                    overall_max = new
        views_out.append({
            "section": section, "name": name,
            "components": sections[0]["components"],
            "steps": len(sections),
            "per_step": per_step,
            "overall": {"min": overall_min, "max": overall_max,
                        "nan": overall_nan, "inf": overall_inf},
        })

    return {"ok": True, "path": str(path), "views": views_out}


# ======================================================================
# Public API: directory audit
# ======================================================================

def audit_msh_directory(directory: str | Path,
                        check_jacobians: bool = False,
                        pattern: str = "**/*.msh",
                        limit: int = 500) -> dict[str, Any]:
    """Validate every .msh under a directory and summarize the health.

    The .msh companion of the Python-lint ``gmsh_audit_summary``: one
    call answers "are the repository's mesh artifacts structurally
    sound?".  ``check_jacobians=True`` additionally runs the gmsh
    getJacobians gate per file (slower; needs the gmsh package).
    """
    root = Path(directory)
    if not root.is_dir():
        return {"ok": False, "directory": str(root),
                "error": f"directory not found: {root}"}

    files = sorted(root.glob(pattern))[: max(1, int(limit))]
    entries: list[dict[str, Any]] = []
    issues: list[dict[str, Any]] = []
    by_status: Counter[str] = Counter()

    for path in files:
        result = validate_msh(path, check_jacobians=check_jacobians)
        status = result.get("status", "needs_attention")
        by_status[status] += 1
        try:
            rel = str(path.relative_to(root))
        except ValueError:
            rel = str(path)
        info = inspect_msh(path)
        entry: dict[str, Any] = {
            "path": rel,
            "status": status,
            "version": info.get("version"),
            "nodes": info.get("nodes", {}).get("count"),
            "elements": info.get("elements", {}).get("count"),
            "views": len(info.get("views") or []),
            "high_order": info.get("high_order", False),
        }
        if check_jacobians and "jacobian" in result:
            entry["negative_jacobians"] = result["jacobian"].get(
                "total_negative")
        entries.append(entry)
        if status != "ok":
            issues.append({
                "path": rel,
                "failed_checks": [k for k, v in result.get("checks", {}).items()
                                  if not v],
                "errors": result.get("errors", [])[:3],
            })

    return {
        "ok": True,
        "directory": str(root),
        "pattern": pattern,
        "files_scanned": len(files),
        "clean": not issues,
        "by_status": dict(by_status),
        "issues": issues,
        "files": entries,
        "jacobians_checked": bool(check_jacobians),
    }


# ======================================================================
# Public API: structural + field diff
# ======================================================================

def _rel_delta(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None if a is b else math.inf
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


def diff_msh(msh_a: str | Path, msh_b: str | Path,
             rel_tol: float = 1e-9) -> dict[str, Any]:
    """Compare two MSH v4.1 files: structure and field statistics.

    Built for before/after checks (re-export, node-order fix, solver
    change): reports node/element/physical/view structure differences
    and, for common views, the relative drift of per-view min/max/mean.
    Field values are compared through statistics, not tag-by-tag.
    """
    results = {}
    for key, p in (("a", Path(msh_a)), ("b", Path(msh_b))):
        info = inspect_msh(p)
        if not info.get("ok"):
            return {"ok": False,
                    "error": f"cannot inspect {key}: "
                             f"{info.get('error') or info.get('parse_errors')}",
                    key: str(p)}
        stats = field_stats(p)
        results[key] = (info, stats)

    (ia, fa), (ib, fb) = results["a"], results["b"]
    differences: list[str] = []

    nodes_a = ia.get("nodes", {}).get("count")
    nodes_b = ib.get("nodes", {}).get("count")
    if nodes_a != nodes_b:
        differences.append(f"node count: {nodes_a} vs {nodes_b}")

    by_type_a = {t["name"]: t["count"]
                 for t in ia.get("elements", {}).get("by_type", [])}
    by_type_b = {t["name"]: t["count"]
                 for t in ib.get("elements", {}).get("by_type", [])}
    for name in sorted(set(by_type_a) | set(by_type_b)):
        ca, cb = by_type_a.get(name, 0), by_type_b.get(name, 0)
        if ca != cb:
            differences.append(f"element {name}: {ca} vs {cb}")

    phys_a = {(p["dim"], p["name"]) for p in ia.get("physical_names") or []}
    phys_b = {(p["dim"], p["name"]) for p in ib.get("physical_names") or []}
    only_phys_a = sorted(f"dim{d} {n}" for d, n in phys_a - phys_b)
    only_phys_b = sorted(f"dim{d} {n}" for d, n in phys_b - phys_a)
    if only_phys_a:
        differences.append(f"physical names only in a: {only_phys_a}")
    if only_phys_b:
        differences.append(f"physical names only in b: {only_phys_b}")

    bbox_max_delta = None
    box_a, box_b = ia.get("bbox"), ib.get("bbox")
    if box_a and box_b:
        bbox_max_delta = max(
            abs(box_a[k][i] - box_b[k][i])
            for k in ("min", "max") for i in range(3))
        if bbox_max_delta > 0:
            differences.append(f"bbox max coordinate delta: {bbox_max_delta:.3e}")

    views_a = {(v["section"], v["name"]): v for v in fa.get("views", [])}
    views_b = {(v["section"], v["name"]): v for v in fb.get("views", [])}
    only_view_a = sorted(name for _, name in set(views_a) - set(views_b))
    only_view_b = sorted(name for _, name in set(views_b) - set(views_a))
    if only_view_a:
        differences.append(f"views only in a: {only_view_a}")
    if only_view_b:
        differences.append(f"views only in b: {only_view_b}")

    structural_diff_count = len(differences)

    common_views = []
    fields_match: bool | None = None
    if set(views_a) & set(views_b):
        fields_match = True
    for key in sorted(set(views_a) & set(views_b)):
        va, vb = views_a[key], views_b[key]
        entry: dict[str, Any] = {
            "section": key[0], "name": key[1],
            "a": {"steps": va["steps"], "components": va["components"],
                  **va["overall"]},
            "b": {"steps": vb["steps"], "components": vb["components"],
                  **vb["overall"]},
        }
        if va["steps"] != vb["steps"]:
            differences.append(
                f"view \"{key[1]}\": {va['steps']} vs {vb['steps']} steps")
            structural_diff_count += 1
        if va["components"] != vb["components"]:
            differences.append(
                f"view \"{key[1]}\": {va['components']} vs "
                f"{vb['components']} components")
            structural_diff_count += 1
        deltas = [d for d in (_rel_delta(va["overall"]["min"], vb["overall"]["min"]),
                              _rel_delta(va["overall"]["max"], vb["overall"]["max"]))
                  if d is not None]
        max_rel = max(deltas) if deltas else None
        entry["max_rel_delta"] = max_rel
        if max_rel is not None and max_rel > rel_tol:
            fields_match = False
            differences.append(
                f"view \"{key[1]}\": min/max drift rel {max_rel:.3e} "
                f"(a: [{va['overall']['min']}, {va['overall']['max']}], "
                f"b: [{vb['overall']['min']}, {vb['overall']['max']}])")
        common_views.append(entry)

    identical_structure = structural_diff_count == 0
    return {
        "ok": True,
        "a": str(msh_a),
        "b": str(msh_b),
        "identical_structure": identical_structure,
        "fields_match": fields_match,
        "rel_tol": rel_tol,
        "differences": _cap(differences),
        "nodes": {"a": nodes_a, "b": nodes_b},
        "elements": {"a": by_type_a, "b": by_type_b},
        "bbox_max_delta": bbox_max_delta,
        "views": {"only_a": only_view_a, "only_b": only_view_b,
                  "common": common_views},
    }


# ======================================================================
# Public API: validate .geo
# ======================================================================

def _count_views_light(msh: Path) -> int:
    """Count distinct (section, name) view groups without a full parse.

    GMSH folds consecutive same-named data sections of one file into one
    view with multiple time steps, so the group count is the number of
    views that file contributes after a Merge.
    """
    groups: set[tuple[str, str]] = set()
    section: str | None = None
    state = 0  # 0 idle, 1 expect numStringTags, 2 expect view name
    try:
        with open(msh, encoding="utf-8", errors="replace") as f:
            for raw in f:
                s = raw.strip()
                if not s:
                    continue
                if s.startswith("$End"):
                    section, state = None, 0
                    continue
                if s.startswith("$"):
                    name = s[1:]
                    section = name if name in _DATA_SECTIONS else None
                    state = 1 if section else 0
                    continue
                if state == 1 and section:
                    try:
                        n_str = int(s)
                    except ValueError:
                        state = 0
                        continue
                    if n_str <= 0:
                        groups.add((section, ""))
                        state = 0
                    else:
                        state = 2
                    continue
                if state == 2 and section:
                    groups.add((section, s.strip('"')))
                    state = 0
    except OSError:
        return 0
    return len(groups)


def validate_geo(geo_path: str | Path, deep: bool = True) -> dict[str, Any]:
    """Validate a .geo launch/companion file before opening it in GMSH.

    Checks that every ``Merge "..."`` target exists on disk and that no
    known-invalid GMSH 4.x option names are used.  With ``deep=True``
    (default) the merged .msh files are scanned for their view count so
    out-of-range ``View[N]`` references -- the classic "opens black"
    bug -- are caught, and the exact-autoload sidecars are reported.
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

    result: dict[str, Any] = {
        "path": str(path),
        "merge_targets": merge_targets,
        "invalid_options": invalid_options,
        "max_view_index": max_view_index,
        "numsubedges": numsubedges,
    }

    if deep:
        total_views = 0
        for entry in merge_targets:
            if not entry["exists"]:
                continue
            resolved = Path(entry["resolved"])
            if resolved.suffix.lower() in (".msh", ".pos"):
                n = _count_views_light(resolved)
                entry["views"] = n
                total_views += n
        result["merged_views_total"] = total_views
        if max_view_index >= 0:
            in_range = max_view_index < total_views
            checks["view_indices_in_range"] = in_range
            if not in_range:
                errors.append(
                    f"View[{max_view_index}] is referenced but the merged "
                    f"files only contribute {total_views} view(s) -- "
                    f"out-of-range view options are silently ignored and "
                    f"the intended field stays invisible")

        geo_opt = Path(str(path) + ".opt")
        sidecars = {"geo_opt": geo_opt.is_file()}
        msh_targets = [Path(e["resolved"]) for e in merge_targets
                       if e["exists"] and e["resolved"].lower().endswith(".msh")]
        if msh_targets:
            sidecars["msh_opt"] = all(
                Path(str(t) + ".opt").is_file() for t in msh_targets)
        result["sidecars"] = sidecars
        if not geo_opt.is_file():
            warnings.append(
                f"no exact-autoload sidecar {geo_opt.name}: display options "
                f"must live inline in the .geo for double-click launches")

    status = "ok" if all(checks.values()) else "needs_attention"
    result.update({
        "ok": status == "ok",
        "status": status,
        "checks": checks,
        "errors": _cap(errors),
        "warnings": warnings,
    })
    return result


# ======================================================================
# CLI:  python -m radia_mcp.gmsh.msh_inspect <target> [options]
# ======================================================================

def main(argv: list[str] | None = None) -> int:
    """CI-friendly command line over inspect/validate/diff/audit.

    Exit code 0 = ok, 1 = needs_attention/error, 2 = usage error.
    """
    import argparse
    import json as _json

    parser = argparse.ArgumentParser(
        prog="python -m radia_mcp.gmsh.msh_inspect",
        description="Inspect / validate / diff GMSH MSH v4.1 artifacts "
                    "(pure Python; --jacobians shells out to gmsh).")
    parser.add_argument("target", help=".msh or .geo file, or a directory")
    parser.add_argument("--validate", action="store_true",
                        help="run validation (default for .geo and dirs)")
    parser.add_argument("--jacobians", action="store_true",
                        help="also run the gmsh getJacobians gate")
    parser.add_argument("--diff", metavar="OTHER",
                        help="diff target .msh against OTHER .msh")
    parser.add_argument("--stats", action="store_true",
                        help="print per-view field statistics")
    parser.add_argument("--json", action="store_true",
                        help="print the full JSON result")
    args = parser.parse_args(argv)

    target = Path(args.target)
    if args.diff:
        result = diff_msh(target, Path(args.diff))
        ok = bool(result.get("ok")) and result.get("identical_structure") \
            and result.get("fields_match") is not False
    elif target.is_dir():
        result = audit_msh_directory(target, check_jacobians=args.jacobians)
        ok = bool(result.get("ok")) and bool(result.get("clean"))
    elif target.suffix.lower() == ".geo":
        result = validate_geo(target)
        ok = bool(result.get("ok"))
    elif args.stats:
        result = field_stats(target)
        ok = bool(result.get("ok"))
    elif args.validate or args.jacobians:
        result = validate_msh(target, check_jacobians=args.jacobians)
        ok = bool(result.get("ok"))
    else:
        result = inspect_msh(target)
        ok = bool(result.get("ok"))

    if args.json:
        print(_json.dumps(result, indent=2, default=str))
    else:
        _print_cli_summary(result)
    return 0 if ok else 1


def _print_cli_summary(result: dict[str, Any]) -> None:
    status = result.get("status") or ("ok" if result.get("ok") else "error")
    print(f"[{status}] {result.get('path') or result.get('directory') or ''}")
    for key in ("version", "high_order", "merged_views_total",
                "files_scanned", "identical_structure", "fields_match"):
        if key in result:
            print(f"  {key}: {result[key]}")
    if result.get("nodes") is not None and isinstance(result.get("nodes"), dict):
        print(f"  nodes: {result['nodes'].get('count')}")
    if isinstance(result.get("elements"), dict) and "by_type" in result["elements"]:
        for t in result["elements"]["by_type"]:
            print(f"  {t['name']} x {t['count']} (order {t['order']})")
    for view in result.get("views") or []:
        if isinstance(view, dict) and "name" in view and "steps" in view:
            print(f"  view \"{view['name']}\": {view['steps']} step(s), "
                  f"{view.get('components')} comp")
    for check, value in (result.get("checks") or {}).items():
        if not value:
            print(f"  FAILED: {check}")
    for issue in result.get("issues") or []:
        print(f"  ISSUE {issue['path']}: {issue['failed_checks']}")
    for err in result.get("errors") or []:
        print(f"  ERROR: {err}")
    for warning in result.get("warnings") or []:
        print(f"  warning: {warning}")
    if result.get("error"):
        print(f"  ERROR: {result['error']}")
    for diff_line in result.get("differences") or []:
        print(f"  diff: {diff_line}")


if __name__ == "__main__":
    import sys as _sys

    _sys.exit(main())
