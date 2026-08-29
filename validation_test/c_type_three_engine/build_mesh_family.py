"""Build a checked, geometric Cubit mesh family for C-type convergence.

The three levels preserve one exact ACIS C-yoke and one Kelvin construction.
Only the Cubit size controls change.  The default scale sequence
``1.25, 1.0, 0.8`` has a common refinement ratio of 1.25, which permits an
observed-order and Richardson analysis without pretending that the unstructured
meshes are nested.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

from build_cubit_meshes import CUBIT_DEFAULT, build


HERE = Path(__file__).resolve().parent
DEFAULT_LEVELS = (
    ("coarse", 1.25),
    ("medium", 1.00),
    ("fine", 0.80),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_level(value: str) -> tuple[str, float]:
    try:
        name, raw_scale = value.split("=", 1)
        scale = float(raw_scale)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("level must be NAME=SCALE") from exc
    name = name.strip()
    if not name or scale <= 0.0:
        raise argparse.ArgumentTypeError("level name and scale must be positive")
    return name, scale


def _load_matching_level(path: Path, expected: dict[str, float | int]) -> dict:
    result_path = path / "mesh_result.json"
    if not result_path.is_file():
        raise FileNotFoundError(result_path)
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    if not payload.get("passed", False):
        raise RuntimeError(f"existing mesh level is not passing: {result_path}")
    for key, value in expected.items():
        actual = payload.get(key)
        if isinstance(value, float):
            matches = actual is not None and abs(float(actual) - value) <= 1e-14
        else:
            matches = actual == value
        if not matches:
            raise RuntimeError(
                f"existing mesh level {path.name!r} has {key}={actual!r}; "
                f"expected {value!r}"
            )
    for artifact_name, hash_name in (
        ("iron_vol", "iron_vol_sha256"),
        ("kelvin_domain_vol", "kelvin_domain_vol_sha256"),
    ):
        artifact = path / Path(payload["artifacts"][artifact_name]).name
        if not artifact.is_file() or sha256(artifact) != payload["artifacts"][hash_name]:
            raise RuntimeError(f"existing mesh artifact hash mismatch: {artifact}")
    return payload


def build_family(options: argparse.Namespace) -> dict[str, object]:
    output_dir = options.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    levels = tuple(options.level or DEFAULT_LEVELS)
    names = [name for name, _ in levels]
    if len(set(names)) != len(names):
        raise ValueError("mesh level names must be unique")
    scales = [float(scale) for _, scale in levels]
    if len(scales) < 3 or any(
        left <= right for left, right in zip(scales, scales[1:])
    ):
        raise ValueError("provide at least three levels ordered coarse to fine")
    ratios = [left / right for left, right in zip(scales, scales[1:])]
    if max(ratios) - min(ratios) > 1e-12:
        raise ValueError("mesh scales must form one geometric sequence")

    rows = []
    for name, scale in levels:
        level_dir = output_dir / name
        expected = {
            "iron_size_m": float(options.iron_size * scale),
            "air_size_m": float(options.air_size * scale),
            "gap_size_m": float(options.gap_size * scale),
            "kelvin_mesh_size_m": float(options.kelvin_mesh_size * scale),
            "kelvin_radius_m": float(options.kelvin_radius),
            "curve_order": int(options.curve_order),
        }
        if options.reuse_existing and (level_dir / "mesh_result.json").is_file():
            payload = _load_matching_level(level_dir, expected)
        else:
            payload = build(
                SimpleNamespace(
                    output_dir=level_dir,
                    cubit=options.cubit,
                    command_plugin_dir=options.command_plugin_dir,
                    iron_size=expected["iron_size_m"],
                    air_size=expected["air_size_m"],
                    gap_size=expected["gap_size_m"],
                    kelvin_radius=expected["kelvin_radius_m"],
                    kelvin_mesh_size=expected["kelvin_mesh_size_m"],
                    curve_order=expected["curve_order"],
                )
            )
        rows.append(
            {
                "name": name,
                "scale": float(scale),
                "relative_directory": name,
                "mesh_result": f"{name}/mesh_result.json",
                "mesh_result_sha256": sha256(level_dir / "mesh_result.json"),
                "iron_size_m": expected["iron_size_m"],
                "air_size_m": expected["air_size_m"],
                "gap_size_m": expected["gap_size_m"],
                "kelvin_mesh_size_m": expected["kelvin_mesh_size_m"],
                "iron_elements": int(payload["inventory"]["iron"]["elements"]),
                "kelvin_elements": int(
                    payload["inventory"]["kelvin_domain"]["elements"]
                ),
                "iron_vol_sha256": payload["artifacts"]["iron_vol_sha256"],
                "kelvin_domain_vol_sha256": payload["artifacts"][
                    "kelvin_domain_vol_sha256"
                ],
            }
        )

    monotone_elements = all(
        left[family] < right[family]
        for left, right in zip(rows, rows[1:])
        for family in ("iron_elements", "kelvin_elements")
    )
    result = {
        "schema": "radia.validation.c-type-cubit-mesh-family.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "passed": bool(monotone_elements),
        "cad_authority": "cad/c_type_iron.jou",
        "builder": "build_cubit_meshes.py",
        "builder_sha256": sha256(HERE / "build_cubit_meshes.py"),
        "family_builder_sha256": sha256(Path(__file__).resolve()),
        "refinement_ratio": float(sum(ratios) / len(ratios)),
        "curve_order": int(options.curve_order),
        "kelvin_radius_m": float(options.kelvin_radius),
        "levels": rows,
        "checks": {
            "all_mesh_contracts_passed": True,
            "element_counts_strictly_increase": monotone_elements,
            "geometric_size_sequence": True,
        },
    }
    manifest = output_dir / "mesh_family.json"
    manifest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not result["passed"]:
        raise RuntimeError(f"C-type mesh family contract failed; see {manifest}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cubit", type=Path, default=CUBIT_DEFAULT)
    parser.add_argument("--command-plugin-dir", type=Path)
    parser.add_argument("--level", action="append", type=_parse_level)
    parser.add_argument("--iron-size", type=float, default=0.030)
    parser.add_argument("--air-size", type=float, default=0.040)
    parser.add_argument("--gap-size", type=float, default=0.002)
    parser.add_argument("--kelvin-radius", type=float, default=0.22)
    parser.add_argument("--kelvin-mesh-size", type=float, default=0.050)
    parser.add_argument("--curve-order", type=int, choices=range(2, 6), default=2)
    parser.add_argument("--reuse-existing", action="store_true")
    options = parser.parse_args()
    if any(
        value <= 0.0
        for value in (
            options.iron_size,
            options.air_size,
            options.gap_size,
            options.kelvin_radius,
            options.kelvin_mesh_size,
        )
    ):
        raise ValueError("mesh sizes must be positive")
    build_family(options)


if __name__ == "__main__":
    main()
