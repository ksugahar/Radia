"""Cubit 2025.12 builder for ESRF Example #3's three-engine FEM mesh.

The HDiv-MMM response uses ``model.vol`` and the independent permanent-magnet
source uses ``fixed_magnetization_sources/magnet_source.vol``.  The HCurl and
mixed-H1 routes instead require a conforming physical air/iron mesh plus the
Kelvin exterior.  PM volumes remain in that mesh as ``air``: their rigid given
magnetization is supplied by the same C++ ``MagnetizationSource`` field, not
reintroduced as a material response unknown.

This is intentionally a Cubit path.  Do not substitute Netgen/OCC geometry in
an acceptance calculation; it would change the CAD and invalidate a three-way
comparison with the HDiv response mesh.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPOSITORY_ROOT = HERE.parents[1]
DEFAULT_CUBIT = Path(
    r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.com"
)
DEFAULT_KELVIN_HELPER = (
    REPOSITORY_ROOT
    / "packages"
    / "cubit-mesh-export"
    / "src"
    / "cubit_mesh_export"
    / "cubit_helpers"
    / "add_kelvin.py"
)


def three_engine_material_contract() -> dict[str, object]:
    """Return the non-negotiable material/source partition for Example #3."""
    return {
        "response_materials": ("iron", "air", "kelvin"),
        "hdiv_response_mesh": "model.vol",
        "fixed_magnetization_source_mesh": (
            "fixed_magnetization_sources/magnet_source.vol"
        ),
        "fem_mesh": "hybrid_undulator_kelvin.vol",
        "pm_fem_material": "air",
        "pm_source": "fixed-given MagnetizationSource",
        "mixed_h1_source_potential": "global_physical",
        "required_boundaries": (
            "iron_air_interface",
            "kelvin_int",
            "kelvin_ext",
        ),
        "required_bbboundary": "GND",
    }


def _volume_ids(cubit) -> set[int]:
    return set(cubit.parse_cubit_list("volume", "all"))


def _new_volume_ids(cubit, before: set[int], label: str) -> set[int]:
    created = _volume_ids(cubit) - before
    if not created:
        raise RuntimeError(f"Cubit did not create {label}")
    return created


def _reflect_meshed_volume(cubit, volume: int, *, axis: str = "z") -> int:
    before = _volume_ids(cubit)
    cubit.cmd(f"volume {volume} copy reflect {axis}")
    reflected = _new_volume_ids(cubit, before, f"reflected volume {volume}")
    if len(reflected) != 1:
        raise RuntimeError(
            f"Expected one reflected copy of volume {volume}, got {sorted(reflected)}"
        )
    return next(iter(reflected))


def _positive_z_volumes(cubit, volumes: set[int], label: str) -> set[int]:
    tolerance = 1.0e-10
    positive = {
        volume for volume in volumes
        if float(cubit.volume(volume).centroid()[2]) > tolerance
    }
    if len(positive) * 2 != len(volumes):
        raise RuntimeError(
            f"{label} must contain z-reflection volume pairs without a "
            f"mid-plane body; all={sorted(volumes)}, positive={sorted(positive)}"
        )
    return positive


def _block_id(cubit, name: str) -> int:
    for block_id in cubit.parse_cubit_list("block", "all"):
        try:
            block_name = cubit.get_block_name(block_id) or ""
        except Exception:
            block_name = cubit.get_exodus_entity_name("block", block_id) or ""
        if block_name.lower() == name.lower():
            return int(block_id)
    raise RuntimeError(f"Cubit block {name!r} was not found")


def _add_block(cubit, block_id: int, name: str, volumes: set[int]) -> None:
    if not volumes:
        raise RuntimeError(f"Cannot create empty Cubit block {name!r}")
    cubit.cmd(f"block {block_id} add volume {' '.join(map(str, sorted(volumes)))}")
    cubit.cmd(f'block {block_id} name "{name}"')


def _add_sideset(cubit, sideset_id: int, name: str, surfaces: set[int]) -> None:
    if not surfaces:
        raise RuntimeError(f"Cannot create empty Cubit sideset {name!r}")
    cubit.cmd(
        f"sideset {sideset_id} add surface {' '.join(map(str, sorted(surfaces)))}"
    )
    cubit.cmd(f'sideset {sideset_id} name "{name}"')


def build_inside_cubit(
    *,
    iron_step: str,
    magnet_step: str,
    kelvin_radius_m: float,
    iron_size_m: float,
    air_size_m: float,
    kelvin_size_m: float,
    kelvin_helper: str | None = None,
) -> dict[str, object]:
    """Create a conforming full-domain FEM mesh in the active Cubit session."""
    import cubit

    if min(kelvin_radius_m, iron_size_m, air_size_m, kelvin_size_m) <= 0.0:
        raise ValueError("all mesh sizes and the Kelvin radius must be positive")

    cubit.cmd("reset")
    before_iron = _volume_ids(cubit)
    cubit.cmd(f'import step "{Path(iron_step).as_posix()}" heal')
    iron_volumes = _new_volume_ids(cubit, before_iron, "iron STEP volumes")
    before_magnet = _volume_ids(cubit)
    cubit.cmd(f'import step "{Path(magnet_step).as_posix()}" heal')
    magnet_volumes = _new_volume_ids(
        cubit, before_magnet, "permanent-magnet STEP volumes"
    )

    # Meshing two CAD-reflected halves independently gives slightly different
    # nodes at z=0.  Kelvin identification then has no one-to-one map.  Keep
    # the positive half, mesh it once, and copy the *meshed* volumes exactly
    # through the legacy z-symmetry plane.
    iron_top = _positive_z_volumes(cubit, iron_volumes, "iron")
    magnet_top = _positive_z_volumes(cubit, magnet_volumes, "permanent magnet")
    deleted = (iron_volumes | magnet_volumes) - iron_top - magnet_top
    cubit.cmd("delete volume " + " ".join(map(str, sorted(deleted))))

    before_sphere = _volume_ids(cubit)
    cubit.cmd(f"create sphere radius {kelvin_radius_m:.17g}")
    air_seed = _new_volume_ids(cubit, before_sphere, "physical air sphere")
    if len(air_seed) != 1:
        raise RuntimeError(f"Expected one physical air sphere, got {sorted(air_seed)}")
    air_seed_id = next(iter(air_seed))
    # A full sphere has no boundary curve from which add_kelvin_cubit can
    # anchor a one-to-one copy mesh.  Retain a z=0 *mesh seam* (not a reduced
    # domain): Example #3 has the legacy z reflection, and the helper then
    # copies/reflects the Kelvin hemisphere without perturbing the full CAD.
    before_cut = _volume_ids(cubit)
    cubit.cmd(f"webcut volume {air_seed_id} with plane zplane")
    air_seeds = {air_seed_id} | (_volume_ids(cubit) - before_cut)
    if len(air_seeds) != 2:
        raise RuntimeError(
            "Expected two physical air hemispheres after the Kelvin mesh seam, "
            f"got {sorted(air_seeds)}"
        )
    air_top_seed = max(
        air_seeds, key=lambda volume: float(cubit.volume(volume).centroid()[2])
    )
    air_bottom_seed = min(
        air_seeds, key=lambda volume: float(cubit.volume(volume).centroid()[2])
    )
    cubit.cmd(f"delete volume {air_bottom_seed}")
    cubit.cmd(
        "subtract volume " + " ".join(map(str, sorted(iron_top | magnet_top)))
        + f" from volume {air_top_seed} keep_tool"
    )
    candidates = _volume_ids(cubit) - iron_top - magnet_top
    if len(candidates) != 1:
        raise RuntimeError(
            "Expected one positive-z physical air hemisphere after subtracting "
            f"the iron/PM solids, got {sorted(candidates)}"
        )
    air_top = next(iter(candidates))

    physical_top = iron_top | magnet_top | {air_top}
    physical_text = " ".join(map(str, sorted(physical_top)))
    cubit.cmd(f"imprint volume {physical_text}")
    cubit.cmd(f"merge volume {physical_text}")

    cubit.cmd("volume " + " ".join(map(str, sorted(physical_top))) + " scheme tetmesh")
    cubit.cmd("volume " + " ".join(map(str, sorted(iron_top)))
              + f" size {iron_size_m:.17g}")
    cubit.cmd("volume " + " ".join(map(str, sorted(magnet_top | {air_top})))
              + f" size {air_size_m:.17g}")
    cubit.cmd("mesh volume " + " ".join(map(str, sorted(physical_top))))

    iron_bottom = {
        _reflect_meshed_volume(cubit, volume) for volume in sorted(iron_top)
    }
    magnet_bottom = {
        _reflect_meshed_volume(cubit, volume) for volume in sorted(magnet_top)
    }
    air_bottom = _reflect_meshed_volume(cubit, air_top)
    iron_volumes = iron_top | iron_bottom
    magnet_volumes = magnet_top | magnet_bottom
    air_volumes = {air_top, air_bottom}
    physical = iron_volumes | magnet_volumes | air_volumes
    physical_text = " ".join(map(str, sorted(physical)))
    cubit.cmd(f"imprint volume {physical_text}")
    cubit.cmd(f"merge volume {physical_text}")

    # add_kelvin_cubit must see only the exterior host in ``air`` while it
    # locates the sphere.  PM volumes join that material only afterwards.
    _add_block(cubit, 1, "air", air_volumes)
    _add_block(cubit, 2, "iron", iron_volumes)
    iron_surfaces = {
        int(surface)
        for volume in iron_volumes
        for surface in cubit.get_relatives("volume", volume, "surface")
    }
    _add_sideset(cubit, 1, "iron_air_interface", iron_surfaces)

    helper = DEFAULT_KELVIN_HELPER if kelvin_helper is None else Path(kelvin_helper)
    if not helper.is_file():
        raise FileNotFoundError(helper)
    spec = importlib.util.spec_from_file_location("_radia_add_kelvin", helper)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {helper}")
    kelvin_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kelvin_module)
    kelvin_module.add_kelvin_cubit(
        R=float(kelvin_radius_m),
        air_block="air",
        symmetry=["z"],
        mesh_size=float(kelvin_size_m),
        kelvin_block="kelvin",
    )
    air_block = _block_id(cubit, "air")
    cubit.cmd(
        f"block {air_block} add volume {' '.join(map(str, sorted(magnet_volumes)))}"
    )
    return {
        "iron_volumes": sorted(iron_volumes),
        "magnet_volumes": sorted(magnet_volumes),
        "physical_air_volumes": sorted(air_volumes),
        "contract": three_engine_material_contract(),
    }


def write_journal(
    output_dir: Path,
    *,
    assets_dir: Path,
    kelvin_radius_m: float,
    iron_size_m: float,
    air_size_m: float,
    kelvin_size_m: float,
    curve_order: int,
    kelvin_helper: Path = DEFAULT_KELVIN_HELPER,
) -> tuple[Path, Path]:
    """Write a self-contained Cubit journal and its in-session Python helper."""
    output_dir.mkdir(parents=True, exist_ok=True)
    iron_step = assets_dir / "iron.step"
    magnet_step = assets_dir / "magnet.step"
    for path in (iron_step, magnet_step):
        if not path.is_file():
            raise FileNotFoundError(path)
    inside = output_dir / "build_hybrid_undulator_inside_cubit.py"
    inside.write_text(
        "import importlib.util\n"
        f'_spec = importlib.util.spec_from_file_location("_hybrid_builder", r"{Path(__file__).resolve()}")\n'
        "_module = importlib.util.module_from_spec(_spec)\n"
        "_spec.loader.exec_module(_module)\n"
        "_result = _module.build_inside_cubit("
        f'iron_step=r"{iron_step.resolve()}", '
        f'magnet_step=r"{magnet_step.resolve()}", '
        f"kelvin_radius_m={float(kelvin_radius_m):.17g}, "
        f"iron_size_m={float(iron_size_m):.17g}, "
        f"air_size_m={float(air_size_m):.17g}, "
        f"kelvin_size_m={float(kelvin_size_m):.17g}, "
        f'kelvin_helper=r"{kelvin_helper.resolve()}")\n'
        "assert _result['contract']['pm_fem_material'] == 'air'\n",
        encoding="utf-8",
    )
    vol = output_dir / "hybrid_undulator_kelvin.vol"
    journal = output_dir / "build_hybrid_undulator_kelvin.jou"
    journal.write_text(
        f'play "{inside.as_posix()}"\n'
        f'export netgen "{vol.as_posix()}" order {int(curve_order)} overwrite\n'
        "exit\n",
        encoding="utf-8",
    )
    return journal, vol


def validate_fem_mesh(path: Path) -> dict[str, object]:
    """Validate labels and Kelvin identification before a three-way solve."""
    import ngsolve as ng
    import numpy as np

    mesh = ng.Mesh(str(path))
    contract = three_engine_material_contract()
    materials = set(map(str, mesh.GetMaterials()))
    boundaries = set(map(str, mesh.GetBoundaries()))
    missing_materials = set(contract["response_materials"]) - materials
    missing_boundaries = set(contract["required_boundaries"]) - boundaries
    if missing_materials or missing_boundaries:
        raise RuntimeError(
            "invalid hybrid-undulator FEM mesh; "
            f"missing_materials={sorted(missing_materials)}, "
            f"missing_boundaries={sorted(missing_boundaries)}"
        )
    if contract["required_bbboundary"] not in set(map(str, mesh.GetBBBoundaries())):
        raise RuntimeError("hybrid-undulator FEM mesh has no Kelvin GND BBND")
    if not mesh.ngmesh.GetIdentifications():
        raise RuntimeError("hybrid-undulator FEM mesh has no Kelvin identification")
    coordinates = np.asarray([vertex.point for vertex in mesh.vertices], dtype=float)
    reflected_vertex = {
        tuple(np.round(point, 14)): point for point in coordinates
    }
    missing_vertices = 0
    maximum_vertex_error = 0.0
    for point in coordinates:
        reflected = point.copy()
        reflected[2] *= -1.0
        key = tuple(np.round(reflected, 14))
        matching = reflected_vertex.get(key)
        if matching is None:
            missing_vertices += 1
            continue
        maximum_vertex_error = max(
            maximum_vertex_error, float(np.linalg.norm(matching - reflected))
        )
    signatures = set()
    element_count = 0
    for element in mesh.Elements(ng.VOL):
        element_count += 1
        points = [
            tuple(np.round(mesh.vertices[vertex.nr].point, 14))
            for vertex in element.vertices
        ]
        signatures.add((str(element.mat), tuple(sorted(points))))
    missing_elements = 0
    for element in mesh.Elements(ng.VOL):
        reflected_points = []
        for vertex in element.vertices:
            point = np.asarray(mesh.vertices[vertex.nr].point, dtype=float)
            point[2] *= -1.0
            reflected_points.append(tuple(np.round(point, 14)))
        signature = (str(element.mat), tuple(sorted(reflected_points)))
        missing_elements += int(signature not in signatures)
    reflection = {
        "missing_reflected_vertices": int(missing_vertices),
        "maximum_reflected_vertex_error_m": float(maximum_vertex_error),
        "volume_element_count": int(element_count),
        "missing_reflected_volume_elements": int(missing_elements),
    }
    if (
        reflection["missing_reflected_vertices"]
        or reflection["maximum_reflected_vertex_error_m"] > 1.0e-13
        or reflection["missing_reflected_volume_elements"]
    ):
        raise RuntimeError(
            "hybrid-undulator FEM mesh lost its reflection-invariant physical "
            f"mesh contract: {reflection}"
        )
    return {
        "path": str(path),
        "elements": int(mesh.ne),
        "materials": list(mesh.GetMaterials()),
        "boundaries": sorted(boundaries),
        "bbboundaries": list(mesh.GetBBBoundaries()),
        "identification_count": int(len(mesh.ngmesh.GetIdentifications())),
        "reflection": reflection,
        "contract": contract,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cubit", type=Path, default=DEFAULT_CUBIT)
    parser.add_argument("--kelvin-radius", type=float, default=0.18)
    parser.add_argument("--iron-size", type=float, default=0.004)
    parser.add_argument("--air-size", type=float, default=0.006)
    parser.add_argument("--kelvin-size", type=float, default=0.012)
    parser.add_argument("--curve-order", type=int, default=2)
    parser.add_argument("--kelvin-helper", type=Path, default=DEFAULT_KELVIN_HELPER)
    options = parser.parse_args(argv)
    if not options.cubit.is_file():
        raise FileNotFoundError(options.cubit)
    journal, vol = write_journal(
        options.output_dir.resolve(),
        assets_dir=options.assets_dir.resolve(),
        kelvin_radius_m=options.kelvin_radius,
        iron_size_m=options.iron_size,
        air_size_m=options.air_size,
        kelvin_size_m=options.kelvin_size,
        curve_order=options.curve_order,
        kelvin_helper=options.kelvin_helper.resolve(),
    )
    completed = subprocess.run(
        [
            str(options.cubit),
            "-nographics",
            "-batch",
            "-nojournal",
            "-input",
            str(journal),
        ],
        cwd=options.output_dir,
        check=False,
    )
    if not vol.is_file():
        raise RuntimeError(f"Cubit did not create {vol}")
    payload = validate_fem_mesh(vol)
    payload["cubit_returncode"] = int(completed.returncode)
    (options.output_dir / "hybrid_undulator_kelvin.mesh.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
