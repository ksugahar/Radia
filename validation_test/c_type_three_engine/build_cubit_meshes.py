"""Build the canonical Cubit 2025.12 meshes for the C-type comparison.

Two meshes are deliberate:

* ``iron.vol`` contains only the exact ACIS iron and is the HDiv-MMM open-
  boundary input. HDiv-MMM must not acquire an artificial air domain.
* ``kelvin_domain.vol`` contains the same iron, a spherical physical-air
  domain, and the translated Kelvin exterior sphere. Reduced-A and
  Omega-reduced-Omega share this open-boundary mesh.

Both artifacts originate from ``cad/c_type_iron.jou``. Python never authors a
replacement C-yoke with ``netgen.occ``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
CAD_JOURNAL = HERE / "cad" / "c_type_iron.jou"
CUBIT_DEFAULT = Path(
    r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.com"
)
EXACT_IRON_VOLUME_M3 = 1_446_095.333333e-9


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run(command: list[str], cwd: Path, log: Path) -> dict[str, object]:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log.write_text(completed.stdout, encoding="utf-8")
    return {
        "command": command,
        "returncode": int(completed.returncode),
        "runtime_s": float(time.perf_counter() - started),
        "log": str(log),
    }


def _write_journal(
    path: Path,
    output: Path,
    *,
    iron_size_m: float,
    air_size_m: float | None,
    gap_size_m: float,
    kelvin_radius_m: float,
    kelvin_mesh_size_m: float,
    curve_order: int,
) -> None:
    cad_lines = CAD_JOURNAL.read_text(encoding="utf-8").rstrip().splitlines()
    reflect_command = "volume all copy reflect z"
    if not cad_lines or cad_lines[-1].strip().lower() != reflect_command:
        raise RuntimeError(
            f"{CAD_JOURNAL} must end with {reflect_command!r}; the canonical "
            "builder replaces that geometry-only reflection with a reflected "
            "mesh copy"
        )
    positive_cad = "\n".join(cad_lines[:-1]).rstrip()
    physical_helper = HERE / "cubit_reflection_mesh.py"
    physical_script = path.with_suffix(".physical.py")
    physical_script.write_text(
        "import importlib.util\n"
        f'_spec = importlib.util.spec_from_file_location("_radia_c_type_reflection", r"{physical_helper}")\n'
        "_module = importlib.util.module_from_spec(_spec)\n"
        "_spec.loader.exec_module(_module)\n"
        "_physical = _module.build_reflection_invariant_physical_mesh("
        f"iron_size={float(iron_size_m):.17g}, "
        f"air_size={None if air_size_m is None else float(air_size_m)!r}, "
        f"gap_size={float(gap_size_m):.17g}, "
        f"kelvin_radius={float(kelvin_radius_m):.17g})\n"
        "assert _physical is not None\n",
        encoding="utf-8",
    )
    lines = [positive_cad, "", f'play "{physical_script.as_posix()}"']
    if air_size_m is None:
        export = f'export netgen "{output.as_posix()}" order 1 overwrite'
    else:
        helper = (
            HERE.parents[1]
            / "packages"
            / "cubit-mesh-export"
            / "src"
            / "cubit_mesh_export"
            / "cubit_helpers"
            / "add_kelvin.py"
        )
        kelvin_script = path.with_suffix(".kelvin.py")
        kelvin_script.write_text(
            "import importlib.util\n"
            f'_spec = importlib.util.spec_from_file_location("_radia_c_type_kelvin", r"{helper}")\n'
            "_module = importlib.util.module_from_spec(_spec)\n"
            "_spec.loader.exec_module(_module)\n"
            f'_info = _module.add_kelvin_cubit(R={float(kelvin_radius_m):.17g}, air_block="air", symmetry=["z"], mesh_size={float(kelvin_mesh_size_m):.17g}, kelvin_block="kelvin")\n'
            "assert _info is not None\n",
            encoding="utf-8",
        )
        lines.append(f'play "{kelvin_script.as_posix()}"')
        export = f'export netgen "{output.as_posix()}" order {int(curve_order)} overwrite'
    lines.extend([export, "exit"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _check_vol(vol: Path, output_dir: Path, *, kelvin_domain: bool) -> dict[str, object]:
    report = output_dir / f"{vol.stem}.vol-check.json"
    command = [
        "check-vol",
        str(vol),
        "--strict-labels",
        "--required-materials",
        "iron,air,kelvin" if kelvin_domain else "iron",
        "--required-boundaries",
        "iron_air_interface,kelvin_int,kelvin_ext" if kelvin_domain else "iron_boundary",
    ]
    if kelvin_domain:
        command.extend(
            [
                "--air-materials",
                "air,kelvin",
            ]
        )
    command.extend(["--report-json", str(report), "--format", "json"])
    run = _run(command, output_dir, output_dir / f"{vol.stem}.check-vol.log")
    if not report.is_file():
        raise RuntimeError(f"check-vol did not create {report}")
    payload = json.loads(report.read_text(encoding="utf-8"))
    if not payload.get("passed", False):
        raise RuntimeError(f"check-vol rejected {vol}; see {report}")
    return {"run": run, "report": str(report), "payload": payload}


def _ngsolve_inventory(path: Path) -> dict[str, object]:
    import ngsolve as ng

    mesh = ng.Mesh(str(path))
    vertex_counts = sorted({len(element.vertices) for element in mesh.Elements(ng.VOL)})
    return {
        "elements": int(mesh.ne),
        "vertices": int(mesh.nv),
        "materials": list(mesh.GetMaterials()),
        "boundaries": sorted(set(mesh.GetBoundaries())),
        "volume_element_vertex_counts": vertex_counts,
    }


def _reflection_inventory(path: Path) -> dict[str, object]:
    import ngsolve as ng
    import numpy as np

    mesh = ng.Mesh(str(path))
    coordinates = np.asarray([vertex.point for vertex in mesh.vertices], dtype=float)
    vertex_keys = {tuple(np.round(point, 14)) for point in coordinates}
    missing_vertices = 0
    maximum_vertex_error = 0.0
    for point in coordinates:
        reflected = np.asarray((point[0], point[1], -point[2]), dtype=float)
        key = tuple(np.round(reflected, 14))
        if key not in vertex_keys:
            missing_vertices += 1
            continue
        candidates = coordinates[
            np.all(np.isclose(coordinates, reflected, rtol=0.0, atol=1e-13), axis=1)
        ]
        if candidates.size:
            maximum_vertex_error = max(
                maximum_vertex_error,
                float(np.min(np.linalg.norm(candidates - reflected, axis=1))),
            )

    def element_signature(element, *, reflect: bool) -> tuple[object, ...]:
        points = []
        for vertex in element.vertices:
            point = np.asarray(mesh.vertices[vertex.nr].point, dtype=float)
            if reflect:
                point[2] = -point[2]
            points.append(tuple(np.round(point, 14)))
        return str(element.mat), tuple(sorted(points))

    elements = tuple(mesh.Elements(ng.VOL))
    signatures = {element_signature(element, reflect=False) for element in elements}
    missing_elements = sum(
        element_signature(element, reflect=True) not in signatures
        for element in elements
    )
    return {
        "vertex_count": int(len(coordinates)),
        "missing_reflected_vertices": int(missing_vertices),
        "maximum_reflected_vertex_error_m": maximum_vertex_error,
        "volume_element_count": int(len(elements)),
        "missing_reflected_volume_elements": int(missing_elements),
    }


def _kelvin_identification_inventory(path: Path) -> dict[str, object]:
    import ngsolve as ng
    import numpy as np

    mesh = ng.Mesh(str(path))
    pairs = mesh.ngmesh.GetIdentifications()
    points = mesh.ngmesh.Points()
    displacements = np.asarray(
        [
            np.asarray(points[slave].p, dtype=float)
            - np.asarray(points[master].p, dtype=float)
            for master, slave in pairs
        ],
        dtype=float,
    )
    if displacements.size == 0:
        return {
            "pair_count": 0,
            "translation_m": [0.0, 0.0, 0.0],
            "maximum_pair_translation_error_m": None,
        }
    translation = np.mean(displacements, axis=0)
    return {
        "pair_count": int(len(pairs)),
        "translation_m": translation.tolist(),
        "maximum_pair_translation_error_m": float(
            np.max(np.linalg.norm(displacements - translation, axis=1))
        ),
    }


def _kelvin_fes_inventory(path: Path) -> dict[str, object]:
    """Verify that the persisted identification actually couples H1 traces."""
    import ngsolve as ng

    mesh = ng.Mesh(str(path))
    base = ng.H1(mesh, order=1, dirichlet_bbbnd="GND")
    periodic = ng.Periodic(base)
    slaved_dofs = int(sum(base.FreeDofs()) - sum(periodic.FreeDofs()))

    trace = ng.GridFunction(periodic)
    trace.vec[:] = 0.0
    trace.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
    inner_norm = float(
        ng.Integrate(
            trace * trace,
            mesh,
            definedon=mesh.Boundaries("kelvin_int"),
        )
    )
    outer_norm = float(
        ng.Integrate(
            trace * trace,
            mesh,
            definedon=mesh.Boundaries("kelvin_ext"),
        )
    )
    ratio = outer_norm / inner_norm if inner_norm > 0.0 else None
    return {
        "base_h1_ndof": int(base.ndof),
        "periodic_h1_ndof": int(periodic.ndof),
        "slaved_free_dofs": slaved_dofs,
        "kelvin_int_trace_norm": inner_norm,
        "kelvin_ext_trace_norm": outer_norm,
        "trace_norm_ratio": ratio,
    }


def _gap_inventory(path: Path) -> dict[str, object]:
    import ngsolve as ng
    import numpy as np

    mesh = ng.Mesh(str(path))
    elements = 0
    maximum_edge = 0.0
    maximum_z_span = 0.0
    for element in mesh.Elements(ng.VOL):
        material = str(element.mat) if hasattr(element, "mat") else ""
        if material != "air":
            continue
        points = np.asarray(
            [mesh.vertices[vertex.nr].point for vertex in element.vertices],
            dtype=float,
        )
        centroid = np.mean(points, axis=0)
        if not (
            abs(float(centroid[0])) <= 0.0171
            and abs(float(centroid[1])) <= 0.0121
            and abs(float(centroid[2])) <= 0.0051
        ):
            continue
        elements += 1
        maximum_z_span = max(
            maximum_z_span, float(np.max(points[:, 2]) - np.min(points[:, 2]))
        )
        for left in range(len(points)):
            for right in range(left + 1, len(points)):
                maximum_edge = max(
                    maximum_edge,
                    float(np.linalg.norm(points[left] - points[right])),
                )
    return {
        "elements": elements,
        "maximum_edge_m": maximum_edge,
        "maximum_z_span_m": maximum_z_span,
        "gap_height_m": 0.010,
    }


def build(options: argparse.Namespace) -> dict[str, object]:
    import ngsolve as ng

    output_dir = options.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cubit = options.cubit.resolve()
    if not cubit.is_file():
        raise FileNotFoundError(cubit)
    if not CAD_JOURNAL.is_file():
        raise FileNotFoundError(CAD_JOURNAL)
    command_plugin_dir = (
        options.command_plugin_dir.resolve()
        if options.command_plugin_dir is not None
        else None
    )
    if command_plugin_dir is not None:
        command_plugin = command_plugin_dir / "cubit_mesh_export.ccm"
        if not command_plugin.is_file():
            raise FileNotFoundError(command_plugin)

    iron_vol = output_dir / "iron.vol"
    kelvin_vol = output_dir / "kelvin_domain.vol"
    iron_journal = output_dir / "build_iron.jou"
    kelvin_journal = output_dir / "build_kelvin_domain.jou"

    _write_journal(
        iron_journal,
        iron_vol,
        iron_size_m=options.iron_size,
        air_size_m=None,
        gap_size_m=options.gap_size,
        kelvin_radius_m=options.kelvin_radius,
        kelvin_mesh_size_m=options.kelvin_mesh_size,
        curve_order=options.curve_order,
    )
    _write_journal(
        kelvin_journal,
        kelvin_vol,
        iron_size_m=options.iron_size,
        air_size_m=options.air_size,
        gap_size_m=options.gap_size,
        kelvin_radius_m=options.kelvin_radius,
        kelvin_mesh_size_m=options.kelvin_mesh_size,
        curve_order=options.curve_order,
    )

    runs = {}
    for name, journal, vol in (
        ("iron", iron_journal, iron_vol),
        ("kelvin_domain", kelvin_journal, kelvin_vol),
    ):
        command = [str(cubit), "-nographics", "-batch", "-nojournal"]
        if command_plugin_dir is not None:
            command.extend(["-commandplugindir", str(command_plugin_dir)])
        command.extend(["-input", str(journal)])
        run = _run(
            command,
            output_dir,
            output_dir / f"{name}.cubit.log",
        )
        runs[name] = run
        if not vol.is_file():
            raise RuntimeError(f"Cubit did not create {vol}; see {run['log']}")

    checks = {
        "iron": _check_vol(iron_vol, output_dir, kelvin_domain=False),
        "kelvin_domain": _check_vol(kelvin_vol, output_dir, kelvin_domain=True),
    }
    with ng.TaskManager():
        inventories = {
            "iron": _ngsolve_inventory(iron_vol),
            "kelvin_domain": _ngsolve_inventory(kelvin_vol),
        }
        reflection = {
            "iron": _reflection_inventory(iron_vol),
            "kelvin_domain": _reflection_inventory(kelvin_vol),
        }
        kelvin_identification = _kelvin_identification_inventory(kelvin_vol)
        kelvin_fes = _kelvin_fes_inventory(kelvin_vol)
        gap_inventory = _gap_inventory(kelvin_vol)
    for name, inventory in inventories.items():
        if inventory["volume_element_vertex_counts"] != [4]:
            raise RuntimeError(f"{name} is not a pure TET mesh: {inventory}")

    reflection_is_exact = all(
        row["missing_reflected_vertices"] == 0
        and row["maximum_reflected_vertex_error_m"] <= 1e-13
        and row["missing_reflected_volume_elements"] == 0
        for row in reflection.values()
    )
    kelvin_is_periodic = bool(
        kelvin_identification["pair_count"] > 0
        and kelvin_identification["maximum_pair_translation_error_m"] is not None
        and kelvin_identification["maximum_pair_translation_error_m"] <= 1e-12
        and kelvin_fes["slaved_free_dofs"] > 0
        and kelvin_fes["trace_norm_ratio"] is not None
        and abs(kelvin_fes["trace_norm_ratio"] - 1.0) <= 1e-12
    )

    iron_rows = checks["iron"]["payload"].get("materials", [])
    iron_volume = next(
        (float(row["ng_volume"]) for row in iron_rows if row.get("name") == "iron"),
        None,
    )
    volume_error = None if iron_volume is None else (
        (iron_volume - EXACT_IRON_VOLUME_M3) / EXACT_IRON_VOLUME_M3
    )
    gap_is_resolved = bool(
        gap_inventory["elements"] > 0
        and gap_inventory["maximum_z_span_m"] <= 0.5 * gap_inventory["gap_height_m"]
    )
    result = {
        "schema": "radia.validation.c-type-cubit-meshes.v1",
        "passed": bool(
            volume_error is not None
            and abs(volume_error) <= 1e-8
            and gap_is_resolved
            and reflection_is_exact
            and kelvin_is_periodic
        ),
        "machine": platform.node(),
        "cad_authority": str(CAD_JOURNAL),
        "cad_sha256": sha256(CAD_JOURNAL),
        "cubit": str(cubit),
        "command_plugin_dir": (
            str(command_plugin_dir) if command_plugin_dir is not None else None
        ),
        "mesh_policy": {
            "hdiv_mmm": "exact iron-only Cubit/ACIS TET .vol; Coulomb open boundary",
            "reduced_a": "exact iron + locally refined physical air + Kelvin exterior",
            "mixed_total_reduced_omega": (
                "same periodic Kelvin .vol as reduced-A; H1 TOSCA "
                "mixed total/reduced Omega"
            ),
        },
        "kelvin_radius_m": float(options.kelvin_radius),
        "kelvin_physical_center_m": [0.0, 0.0, 0.0],
        "iron_size_m": float(options.iron_size),
        "air_size_m": float(options.air_size),
        "gap_size_m": float(options.gap_size),
        "kelvin_mesh_size_m": float(options.kelvin_mesh_size),
        "curve_order": int(options.curve_order),
        "gap_inventory": gap_inventory,
        "reflection_inventory": reflection,
        "kelvin_identification": kelvin_identification,
        "kelvin_fes": kelvin_fes,
        "exact_iron_volume_m3": EXACT_IRON_VOLUME_M3,
        "mesh_iron_volume_m3": iron_volume,
        "relative_iron_volume_error": volume_error,
        "artifacts": {
            "iron_vol": str(iron_vol),
            "iron_vol_sha256": sha256(iron_vol),
            "kelvin_domain_vol": str(kelvin_vol),
            "kelvin_domain_vol_sha256": sha256(kelvin_vol),
        },
        "inventory": inventories,
        "checks": checks,
        "runs": runs,
    }
    result_path = output_dir / "mesh_result.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not result["passed"]:
        raise RuntimeError(f"Cubit C-type mesh contract failed; see {result_path}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cubit", type=Path, default=CUBIT_DEFAULT)
    parser.add_argument(
        "--command-plugin-dir",
        type=Path,
        help="Cubit command-plugin directory containing cubit_mesh_export.ccm",
    )
    parser.add_argument("--iron-size", type=float, default=0.010)
    parser.add_argument("--air-size", type=float, default=0.015)
    parser.add_argument("--gap-size", type=float, default=0.002)
    parser.add_argument("--kelvin-radius", type=float, default=0.22)
    parser.add_argument("--kelvin-mesh-size", type=float, default=0.025)
    parser.add_argument("--curve-order", type=int, choices=range(2, 6), default=2)
    options = parser.parse_args()
    if any(value <= 0.0 for value in (
        options.iron_size,
        options.air_size,
        options.gap_size,
        options.kelvin_radius,
        options.kelvin_mesh_size,
    )):
        raise ValueError("mesh sizes must be positive")
    build(options)


if __name__ == "__main__":
    main()
