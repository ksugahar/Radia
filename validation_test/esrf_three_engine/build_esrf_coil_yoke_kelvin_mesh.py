"""Build the FEM Kelvin mesh for coil-driven ESRF Examples 6 and 7.

The corresponding HDiv-MMM response mesh is the existing iron-only Q2 Cubit
asset ``model.vol``.  This builder creates its independent FEM partner from
the very same iron STEP: iron, a physical air sphere, and the Kelvin exterior.
Coils remain mesh-free ``CoilBuilder`` sources in every formulation.
"""

from __future__ import annotations

import argparse
import hashlib
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def three_engine_material_contract() -> dict[str, object]:
    """Return the shared material and open-boundary contract."""
    return {
        "response_materials": ("iron", "air", "kelvin"),
        "hdiv_response_mesh": "assets/model.vol",
        "fem_mesh_role": "conforming iron + physical air + Kelvin exterior",
        "coil_source": "shared mesh-free CoilBuilder solid current",
        "mixed_h1_source_potential": "surface_trace_with_cut_gate",
        "required_boundaries": (
            "iron_air_interface",
            "kelvin_int",
            "kelvin_ext",
        ),
        "required_bbboundary": "GND",
        "finite_outer_air_box_forbidden": True,
    }


def _volume_ids(cubit) -> set[int]:
    return set(cubit.parse_cubit_list("volume", "all"))


def _new_volumes(cubit, before: set[int], label: str) -> set[int]:
    result = _volume_ids(cubit) - before
    if not result:
        raise RuntimeError(f"Cubit did not create {label}")
    return result


def _reflect_meshed_volume(cubit, volume: int, *, axis: str = "z") -> int:
    before = _volume_ids(cubit)
    cubit.cmd(f"volume {volume} copy reflect {axis}")
    reflected = _new_volumes(cubit, before, f"reflected volume {volume}")
    if len(reflected) != 1:
        raise RuntimeError(
            f"Expected one reflected copy of volume {volume}, got {sorted(reflected)}"
        )
    return next(iter(reflected))


def _positive_half_volumes(cubit, volumes: set[int], label: str) -> set[int]:
    tolerance = 1.0e-10
    positive = {
        volume
        for volume in volumes
        if float(cubit.volume(volume).centroid()[2]) > tolerance
    }
    if len(positive) * 2 != len(volumes):
        raise RuntimeError(
            f"{label} must split into z-reflection volume pairs; "
            f"all={sorted(volumes)}, positive={sorted(positive)}"
        )
    return positive


def _crosses_z_mesh_seam(cubit, volume: int) -> bool:
    """Return whether a CAD body has interior on both sides of z=0."""
    coordinates = [
        float(cubit.vertex(vertex).coordinates()[2])
        for vertex in cubit.get_relatives("volume", volume, "vertex")
    ]
    if not coordinates:
        raise RuntimeError(f"Cubit volume {volume} has no vertices")
    tolerance = 1.0e-10
    return min(coordinates) < -tolerance and max(coordinates) > tolerance


def _add_block(cubit, block_id: int, name: str, volumes: set[int]) -> None:
    if not volumes:
        raise RuntimeError(f"Cannot create empty block {name!r}")
    cubit.cmd(f"block {block_id} add volume {' '.join(map(str, sorted(volumes)))}")
    cubit.cmd(f'block {block_id} name "{name}"')


def _add_sideset(cubit, sideset_id: int, name: str, surfaces: set[int]) -> None:
    if not surfaces:
        raise RuntimeError(f"Cannot create empty sideset {name!r}")
    cubit.cmd(
        f"sideset {sideset_id} add surface {' '.join(map(str, sorted(surfaces)))}"
    )
    cubit.cmd(f'sideset {sideset_id} name "{name}"')


def _iron_air_surfaces(cubit, iron: set[int], air: set[int]) -> set[int]:
    """Return only physical iron/air faces, excluding iron partitions."""
    surfaces: set[int] = set()
    for volume in iron:
        for surface in cubit.get_relatives("volume", volume, "surface"):
            adjacent = set(cubit.get_relatives("surface", surface, "volume"))
            if adjacent & iron and adjacent & air:
                surfaces.add(int(surface))
    return surfaces


def _gap_interface_surfaces(
    cubit,
    interface: set[int],
    *,
    beam_axis: int,
    radius_m: float,
    half_length_m: float,
) -> set[int]:
    """Select bore-facing interface patches for local physical-air refinement."""
    if beam_axis not in (0, 1, 2):
        raise ValueError("beam_axis must be 0, 1, or 2")
    if radius_m <= 0.0 or half_length_m <= 0.0:
        raise ValueError("gap refinement radius and half length must be positive")
    transverse_axes = tuple(axis for axis in range(3) if axis != beam_axis)
    selected: set[int] = set()
    for surface in interface:
        center = tuple(float(value) for value in cubit.get_center_point("surface", surface))
        transverse_radius = sum(center[axis] ** 2 for axis in transverse_axes) ** 0.5
        if (
            transverse_radius <= radius_m
            and abs(center[beam_axis]) <= half_length_m
        ):
            selected.add(surface)
    if not selected:
        raise RuntimeError(
            "No bore-facing iron/air surfaces matched the requested local "
            "gap-refinement region"
        )
    return selected


def build_inside_cubit(
    *,
    iron_step: str,
    kelvin_radius_m: float,
    iron_size_m: float,
    air_size_m: float,
    kelvin_size_m: float,
    gap_size_m: float,
    gap_refinement_radius_m: float,
    gap_refinement_half_length_m: float,
    beam_axis: int,
    kelvin_helper: str,
) -> dict[str, object]:
    """Construct and mesh the full physical domain in an active Cubit run."""
    import cubit

    if min(
        kelvin_radius_m,
        iron_size_m,
        air_size_m,
        kelvin_size_m,
        gap_size_m,
        gap_refinement_radius_m,
        gap_refinement_half_length_m,
    ) <= 0.0:
        raise ValueError("Kelvin radius and all mesh sizes must be positive")
    cubit.cmd("reset")
    before_iron = _volume_ids(cubit)
    cubit.cmd(f'import step "{Path(iron_step).as_posix()}" heal')
    iron = _new_volumes(cubit, before_iron, "iron STEP volumes")
    # Keep a z=0 *mesh seam* in the full model.  A closed sphere has no
    # boundary curve usable by Cubit's one-to-one Kelvin copy.  Meshing the
    # positive half once and reflecting the completed volumes preserves a
    # full-domain calculation while also making both Kelvin interfaces exact
    # copies of their physical-air partners.
    crossing_iron = {
        volume for volume in iron if _crosses_z_mesh_seam(cubit, volume)
    }
    if crossing_iron:
        cubit.cmd(
            "webcut volume "
            + " ".join(map(str, sorted(crossing_iron)))
            + " with plane zplane"
        )
    iron_halves = _volume_ids(cubit)
    iron_top = _positive_half_volumes(cubit, iron_halves, "iron")
    iron_bottom_seed = iron_halves - iron_top
    cubit.cmd("delete volume " + " ".join(map(str, sorted(iron_bottom_seed))))

    before_air = _volume_ids(cubit)
    cubit.cmd(f"create sphere radius {float(kelvin_radius_m):.17g}")
    air_seed = _new_volumes(cubit, before_air, "physical air sphere")
    if len(air_seed) != 1:
        raise RuntimeError(f"Expected one physical air sphere, got {sorted(air_seed)}")
    air_seed_id = next(iter(air_seed))
    cubit.cmd(f"webcut volume {air_seed_id} with plane zplane")
    air_halves = _volume_ids(cubit) - iron_top
    if air_halves == {air_seed_id} or len(air_halves) != 2:
        raise RuntimeError(
            "Cubit did not split the physical air sphere at the z=0 mesh seam; "
            f"got {sorted(air_halves)}"
        )
    air_top_seed = _positive_half_volumes(cubit, air_halves, "physical air")
    if len(air_top_seed) != 1:
        raise RuntimeError(f"Expected one positive-z air hemisphere, got {sorted(air_top_seed)}")
    cubit.cmd("delete volume " + " ".join(map(str, sorted(air_halves - air_top_seed))))
    air_top_id = next(iter(air_top_seed))
    cubit.cmd(
        "subtract volume "
        + " ".join(map(str, sorted(iron_top)))
        + f" from volume {air_top_id} keep_tool"
    )
    air_top = _volume_ids(cubit) - iron_top
    if len(air_top) != 1:
        raise RuntimeError(
            "Expected one positive-z physical air domain after subtracting iron; "
            f"got {sorted(air_top)}"
        )
    physical_top = iron_top | air_top
    physical_top_text = " ".join(map(str, sorted(physical_top)))
    cubit.cmd(f"imprint volume {physical_top_text}")
    cubit.cmd(f"merge volume {physical_top_text}")
    interface = _iron_air_surfaces(cubit, iron_top, air_top)
    if not interface:
        raise RuntimeError("No physical iron/air interface was found")

    cubit.cmd(f"volume {physical_top_text} scheme tetmesh")
    cubit.cmd("volume " + " ".join(map(str, sorted(iron_top)))
              + f" size {float(iron_size_m):.17g}")
    cubit.cmd("volume " + " ".join(map(str, sorted(air_top)))
              + f" size {float(air_size_m):.17g}")
    gap_surfaces = _gap_interface_surfaces(
        cubit,
        interface,
        beam_axis=int(beam_axis),
        radius_m=float(gap_refinement_radius_m),
        half_length_m=float(gap_refinement_half_length_m),
    )
    cubit.cmd(
        "surface "
        + " ".join(map(str, sorted(gap_surfaces)))
        + f" size {float(gap_size_m):.17g}"
    )
    cubit.cmd(f"mesh volume {physical_top_text}")
    iron_bottom = {
        _reflect_meshed_volume(cubit, volume) for volume in sorted(iron_top)
    }
    air_bottom = {_reflect_meshed_volume(cubit, next(iter(air_top)))}
    iron = iron_top | iron_bottom
    air = air_top | air_bottom
    physical = iron | air
    cubit.cmd("imprint volume " + " ".join(map(str, sorted(physical))))
    cubit.cmd("merge volume " + " ".join(map(str, sorted(physical))))
    _add_block(cubit, 1, "iron", iron)
    _add_block(cubit, 2, "air", air)
    interface = _iron_air_surfaces(cubit, iron, air)
    _add_sideset(cubit, 1, "iron_air_interface", interface)

    helper = Path(kelvin_helper)
    if not helper.is_file():
        raise FileNotFoundError(helper)
    spec = importlib.util.spec_from_file_location("_radia_add_kelvin", helper)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load Kelvin helper {helper}")
    kelvin_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kelvin_module)
    kelvin = kelvin_module.add_kelvin_cubit(
        R=float(kelvin_radius_m),
        air_block="air",
        symmetry=["z"],
        mesh_size=float(kelvin_size_m),
        kelvin_block="kelvin",
    )
    return {
        "iron_volumes": sorted(iron),
        "physical_air_volumes": sorted(air),
        "iron_air_surface_count": int(len(interface)),
        "gap_refinement_surface_count": int(len(gap_surfaces)),
        "kelvin": dict(kelvin),
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
    gap_size_m: float,
    gap_refinement_radius_m: float,
    gap_refinement_half_length_m: float,
    beam_axis: int,
    curve_order: int,
    kelvin_helper: Path = DEFAULT_KELVIN_HELPER,
) -> tuple[Path, Path]:
    """Write a reproducible Cubit journal and active-session helper."""
    output_dir.mkdir(parents=True, exist_ok=True)
    iron_step = assets_dir / "iron.step"
    if not iron_step.is_file():
        raise FileNotFoundError(iron_step)
    if int(curve_order) < 1:
        raise ValueError("curve_order must be positive")
    inside = output_dir / "build_esrf_coil_yoke_inside_cubit.py"
    inside.write_text(
        "import importlib.util\n"
        f'_spec = importlib.util.spec_from_file_location("_esrf_coil_yoke_builder", r"{Path(__file__).resolve()}")\n'
        "_module = importlib.util.module_from_spec(_spec)\n"
        "_spec.loader.exec_module(_module)\n"
        "_result = _module.build_inside_cubit("
        f'iron_step=r"{iron_step.resolve()}", '
        f"kelvin_radius_m={float(kelvin_radius_m):.17g}, "
        f"iron_size_m={float(iron_size_m):.17g}, "
        f"air_size_m={float(air_size_m):.17g}, "
        f"kelvin_size_m={float(kelvin_size_m):.17g}, "
        f"gap_size_m={float(gap_size_m):.17g}, "
        f"gap_refinement_radius_m={float(gap_refinement_radius_m):.17g}, "
        f"gap_refinement_half_length_m={float(gap_refinement_half_length_m):.17g}, "
        f"beam_axis={int(beam_axis)}, "
        f'kelvin_helper=r"{kelvin_helper.resolve()}")\n'
        "assert _result['contract']['finite_outer_air_box_forbidden']\n",
        encoding="utf-8",
    )
    vol = output_dir / "coil_yoke_kelvin.vol"
    journal = output_dir / "build_coil_yoke_kelvin.jou"
    journal.write_text(
        f'play "{inside.as_posix()}"\n'
        f'export netgen "{vol.as_posix()}" order {int(curve_order)} overwrite\n'
        "exit\n",
        encoding="utf-8",
    )
    return journal, vol


def validate_fem_mesh(path: Path) -> dict[str, object]:
    """Reject any FEM mesh that lacks the required Kelvin topology."""
    import ngsolve as ng

    from radia.kelvin_identify_ngsolve import has_kelvin_identification

    mesh = ng.Mesh(str(path))
    contract = three_engine_material_contract()
    materials = set(map(str, mesh.GetMaterials()))
    boundaries = set(map(str, mesh.GetBoundaries()))
    missing_materials = set(contract["response_materials"]) - materials
    missing_boundaries = set(contract["required_boundaries"]) - boundaries
    if missing_materials or missing_boundaries:
        raise RuntimeError(
            "invalid coil-yoke FEM mesh; "
            f"missing_materials={sorted(missing_materials)}, "
            f"missing_boundaries={sorted(missing_boundaries)}"
        )
    if contract["required_bbboundary"] not in set(map(str, mesh.GetBBBoundaries())):
        raise RuntimeError("coil-yoke FEM mesh has no Kelvin GND BBND")
    if not has_kelvin_identification(mesh):
        raise RuntimeError("coil-yoke FEM mesh has no Kelvin identification")
    return {
        "path": str(path),
        "elements": int(mesh.ne),
        "vertices": int(mesh.nv),
        "curve_order": int(mesh.GetCurveOrder()),
        "materials": list(mesh.GetMaterials()),
        "boundaries": sorted(boundaries),
        "bbboundaries": list(mesh.GetBBBoundaries()),
        "identification_count": int(len(mesh.ngmesh.GetIdentifications())),
        "contract": contract,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--case", choices=(6, 7), type=int, required=True)
    parser.add_argument("--cubit", type=Path, default=DEFAULT_CUBIT)
    parser.add_argument("--kelvin-radius", type=float, default=None)
    parser.add_argument("--iron-size", type=float, default=None)
    parser.add_argument("--air-size", type=float, default=None)
    parser.add_argument("--kelvin-size", type=float, default=None)
    parser.add_argument("--gap-size", type=float, default=0.003)
    parser.add_argument("--curve-order", type=int, default=2)
    parser.add_argument("--kelvin-helper", type=Path, default=DEFAULT_KELVIN_HELPER)
    options = parser.parse_args(argv)
    if not options.cubit.is_file():
        raise FileNotFoundError(options.cubit)
    from esrf_coil_yoke import get_case

    case = get_case(options.case)
    journal, vol = write_journal(
        options.output_dir.resolve(),
        assets_dir=options.assets_dir.resolve(),
        kelvin_radius_m=(
            case.kelvin_radius_m
            if options.kelvin_radius is None
            else options.kelvin_radius
        ),
        iron_size_m=(case.iron_size_m if options.iron_size is None else options.iron_size),
        air_size_m=(
            case.outer_air_size_m if options.air_size is None else options.air_size
        ),
        kelvin_size_m=(case.kelvin_size_m if options.kelvin_size is None else options.kelvin_size),
        gap_size_m=options.gap_size,
        gap_refinement_radius_m=case.gap_refinement_radius_m,
        gap_refinement_half_length_m=case.gap_refinement_half_length_m,
        beam_axis=case.beam_axis,
        curve_order=options.curve_order,
        kelvin_helper=options.kelvin_helper.resolve(),
    )
    completed = subprocess.run(
        [str(options.cubit), "-nographics", "-batch", "-nojournal", "-input", str(journal)],
        cwd=options.output_dir,
        check=False,
    )
    if completed.returncode != 0 or not vol.is_file():
        raise RuntimeError(f"Cubit did not create {vol}; returncode={completed.returncode}")
    report = validate_fem_mesh(vol)
    asset_manifest = options.assets_dir.resolve() / "manifest.json"
    if asset_manifest.is_file():
        report["asset_manifest_sha256"] = _sha256(asset_manifest)
    report["cubit_returncode"] = int(completed.returncode)
    output = options.output_dir.resolve() / "coil_yoke_kelvin.mesh.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
