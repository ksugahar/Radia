"""Inspect the docs-local GMSH displacement-animation artifact."""

from __future__ import annotations

import datetime as _dt
import importlib.metadata as _metadata
import json
import math
import platform
import sys
from pathlib import Path


def _pkg_version(name: str) -> str | None:
    try:
        return _metadata.version(name)
    except Exception:
        return None


def runtime_versions() -> dict:
    return {
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "radia_version": _pkg_version("radia"),
    }


def _strip_tag(line: str) -> str:
    return line.strip().strip('"')


def _read_geo_merges(path: Path) -> list[str]:
    merges = []
    if not path.exists():
        return merges
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("Merge "):
            merges.append(line.split('"')[1])
    return merges


def inspect_msh(path: str | Path) -> dict:
    """Inspect mesh-format, node/element headers, and NodeData time series."""
    p = Path(path)
    lines = p.read_text(encoding="utf-8").splitlines()
    mesh_format = None
    nodes_header = None
    elements_header = None
    node_data: list[dict] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "$MeshFormat":
            mesh_format = lines[i + 1].strip()
            i += 3
            continue
        if line == "$Nodes":
            parts = lines[i + 1].split()
            nodes_header = {
                "entity_blocks": int(parts[0]),
                "num_nodes": int(parts[1]),
                "min_node_tag": int(parts[2]),
                "max_node_tag": int(parts[3]),
            }
            i += 1
            continue
        if line == "$Elements":
            parts = lines[i + 1].split()
            elements_header = {
                "entity_blocks": int(parts[0]),
                "num_elements": int(parts[1]),
                "min_element_tag": int(parts[2]),
                "max_element_tag": int(parts[3]),
            }
            i += 1
            continue
        if line == "$NodeData":
            j = i + 1
            n_string = int(lines[j]); j += 1
            string_tags = [_strip_tag(lines[j + k]) for k in range(n_string)]
            j += n_string
            n_real = int(lines[j]); j += 1
            real_tags = [float(lines[j + k]) for k in range(n_real)]
            j += n_real
            n_int = int(lines[j]); j += 1
            int_tags = [int(lines[j + k]) for k in range(n_int)]
            j += n_int
            time_value = real_tags[0] if real_tags else None
            time_step = int_tags[0] if len(int_tags) > 0 else None
            num_components = int_tags[1] if len(int_tags) > 1 else None
            num_entries = int_tags[2] if len(int_tags) > 2 else None
            max_norm = 0.0
            min_norm = None
            sum_norm = 0.0
            for _ in range(num_entries or 0):
                vals = lines[j].split()
                comps = [float(v) for v in vals[1:1 + (num_components or 0)]]
                norm = math.sqrt(sum(c * c for c in comps))
                max_norm = max(max_norm, norm)
                min_norm = norm if min_norm is None else min(min_norm, norm)
                sum_norm += norm
                j += 1
            node_data.append({
                "name": string_tags[0] if string_tags else "",
                "time_value": time_value,
                "time_step": time_step,
                "num_components": num_components,
                "num_entries": num_entries,
                "min_norm": min_norm,
                "max_norm": max_norm,
                "mean_norm": sum_norm / num_entries if num_entries else None,
            })
            i = j + 1
            continue
        i += 1

    return {
        "path": p.as_posix(),
        "bytes": p.stat().st_size,
        "mesh_format": mesh_format,
        "nodes": nodes_header,
        "elements": elements_header,
        "node_data_count": len(node_data),
        "node_data": node_data,
    }


def inspect_example(example_dir: str | Path = ".") -> dict:
    example = Path(example_dir)
    if not example.is_absolute():
        example = (Path(__file__).resolve().parent / example).resolve()
    repo_root = Path(__file__).resolve().parents[2]

    def _repo_rel(path: Path) -> str:
        try:
            return path.resolve().relative_to(repo_root).as_posix()
        except ValueError:
            return path.as_posix()

    msh = example / "rotor_animation.msh"
    geo = example / "animation.geo"
    stator = example / "stator.step"
    rotor = example / "rotor.vol"
    msh_info = inspect_msh(msh)
    msh_info["path"] = _repo_rel(msh)
    time_values = [row["time_value"] for row in msh_info["node_data"]]
    max_disp = [row["max_norm"] for row in msh_info["node_data"]]
    monotone_time = all(a <= b for a, b in zip(time_values, time_values[1:]))
    monotone_disp = all(a <= b + 1e-14 for a, b in zip(max_disp, max_disp[1:]))
    final_displacement_m = max_disp[-1] if max_disp else None
    return {
        "schema": "radia.docs.gmsh_animation.inspect.v1",
        "generated_at_utc": _dt.datetime.now(
            _dt.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "versions": runtime_versions(),
        "example_dir": _repo_rel(example),
        "files": {
            "animation_geo": _repo_rel(geo),
            "rotor_animation_msh": _repo_rel(msh),
            "rotor_vol": _repo_rel(rotor),
            "stator_step": _repo_rel(stator),
        },
        "geo_merges": _read_geo_merges(geo),
        "msh": msh_info,
        "checks": {
            "is_msh_v41": msh_info["mesh_format"] == "4.1 0 8",
            "has_vector_node_data": all(
                row["num_components"] == 3 for row in msh_info["node_data"]
            ),
            "time_values_monotone": monotone_time,
            "max_displacement_monotone": monotone_disp,
            "final_displacement_m": final_displacement_m,
            "expected_final_displacement_m": 0.15,
            "final_displacement_close": (
                abs(final_displacement_m - 0.15) < 1e-12
                if final_displacement_m is not None else False
            ),
        },
    }


def write_results_json(result: dict, path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                 encoding="utf-8")
    return p


def format_summary(result: dict) -> str:
    msh = result["msh"]
    checks = result["checks"]
    nodes = msh["nodes"] or {}
    elements = msh["elements"] or {}
    return "\n".join([
        f"mesh_format: {msh['mesh_format']}",
        f"nodes: {nodes.get('num_nodes')}",
        f"elements: {elements.get('num_elements')}",
        f"NodeData frames: {msh['node_data_count']}",
        f"geo merges: {', '.join(result['geo_merges'])}",
        f"final displacement: {checks['final_displacement_m']:.6g} m",
        f"checks: v4.1={checks['is_msh_v41']}, "
        f"vector={checks['has_vector_node_data']}, "
        f"time_monotone={checks['time_values_monotone']}, "
        f"disp_monotone={checks['max_displacement_monotone']}",
    ])
