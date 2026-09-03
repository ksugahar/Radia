"""OCC sphere volume/area p-convergence validation driver.

This is the Cubit-free numerical lane. The full Cubit ACIS-trampoline
presentation context lives beside it in ``acis_volume_area_convergence_demo.py``.
The durable result is ``volume_area_convergence_results.json``.
"""

from __future__ import annotations

import datetime as _dt
import importlib.metadata as _metadata
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Iterable


RADIUS_M = 0.05
V_EXACT_M3 = (4.0 / 3.0) * math.pi * RADIUS_M**3
A_EXACT_M2 = 4.0 * math.pi * RADIUS_M**2


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
        "ngsolve_version": _pkg_version("ngsolve"),
        "netgen_version": _pkg_version("netgen-mesher"),
        "radia_version": _pkg_version("radia"),
    }


def run_occ_convergence(
    maxh_cases: Iterable[tuple[float, str]] = ((0.015, "coarse"), (0.008, "fine")),
    orders: Iterable[int] = range(1, 6),
    radius_m: float = RADIUS_M,
) -> dict:
    """Run OCC native ``mesh.Curve(p)`` volume/area convergence on a sphere."""
    from netgen.occ import OCCGeometry, Sphere
    from ngsolve import BND, CF, Integrate, Mesh as NGMesh, TaskManager

    v_exact = (4.0 / 3.0) * math.pi * radius_m**3
    a_exact = 4.0 * math.pi * radius_m**2
    geo = OCCGeometry(Sphere((0, 0, 0), radius_m))
    rows: list[dict] = []

    for maxh, label in maxh_cases:
        with TaskManager():
            base = geo.GenerateMesh(maxh=maxh)
            base_ne = int(base.ne)

        for order in orders:
            with TaskManager():
                ng = geo.GenerateMesh(maxh=maxh)
                t0 = time.perf_counter()
                if order >= 2:
                    ng.Curve(order)
                curve_time_s = time.perf_counter() - t0
                mesh = NGMesh(ng)
                volume_m3 = float(Integrate(CF(1), mesh))
                area_m2 = float(Integrate(CF(1), mesh, BND))

            rows.append({
                "case": label,
                "maxh_m": float(maxh),
                "base_elements": base_ne,
                "order": int(order),
                "volume_m3": volume_m3,
                "volume_error_percent": (volume_m3 - v_exact) / v_exact * 100.0,
                "area_m2": area_m2,
                "area_error_percent": (area_m2 - a_exact) / a_exact * 100.0,
                "curve_time_s": curve_time_s,
            })

    return {
        "schema": "radia.validation.ngsolve-user-meeting-volume-area-convergence.v1",
        "generated_at_utc": _dt.datetime.now(
            _dt.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "versions": runtime_versions(),
        "parameters": {
            "radius_m": radius_m,
            "v_exact_m3": v_exact,
            "a_exact_m2": a_exact,
            "maxh_cases": [{"maxh_m": float(m), "label": str(label)}
                           for m, label in maxh_cases],
            "orders": [int(p) for p in orders],
        },
        "rows": rows,
        "summary": summarize(rows),
    }


def summarize(rows: list[dict]) -> dict:
    by_case: dict[str, list[dict]] = {}
    for row in rows:
        by_case.setdefault(row["case"], []).append(row)
    out = {}
    for case, case_rows in by_case.items():
        p1 = next((r for r in case_rows if r["order"] == 1), None)
        p5 = next((r for r in case_rows if r["order"] == 5), None)
        out[case] = {
            "base_elements": case_rows[0]["base_elements"],
            "p1_volume_error_percent": p1["volume_error_percent"] if p1 else None,
            "p5_volume_error_percent": p5["volume_error_percent"] if p5 else None,
            "p1_area_error_percent": p1["area_error_percent"] if p1 else None,
            "p5_area_error_percent": p5["area_error_percent"] if p5 else None,
        }
    return out


def write_results_json(result: dict, path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                 encoding="utf-8")
    return p


def format_rows(rows: list[dict]) -> str:
    lines = [
        f"{'case':<8} {'p':>2} {'ne':>7} {'V err (%)':>12} "
        f"{'A err (%)':>12} {'curve (s)':>10}"
    ]
    for row in rows:
        lines.append(
            f"{row['case']:<8} {row['order']:>2d} {row['base_elements']:>7d} "
            f"{row['volume_error_percent']:>+12.6f} "
            f"{row['area_error_percent']:>+12.6f} "
            f"{row['curve_time_s']:>10.4f}"
        )
    return "\n".join(lines)
