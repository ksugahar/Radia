"""Build a reflection-invariant full Cubit mesh from positive-z CAD.

The C-type ACIS journal authors only the positive-z iron.  This helper meshes
that physical half-domain once and lets Cubit copy the *meshed* volumes through
z=0.  Meshing two pre-reflected geometric halves independently is forbidden:
it breaks the roundoff-level reflection contract used by the three-engine
validation.
"""

from __future__ import annotations


def _volume_ids(cubit) -> set[int]:
    return set(cubit.parse_cubit_list("volume", "all"))


def _single(values: set[int], label: str) -> int:
    if len(values) != 1:
        raise RuntimeError(f"Expected one {label} volume, found {sorted(values)}")
    return next(iter(values))


def _new_volume(cubit, before: set[int], label: str) -> int:
    return _single(_volume_ids(cubit) - before, label)


def _reflect_meshed_volume(cubit, source: int, axis: str = "z") -> int:
    before = _volume_ids(cubit)
    cubit.cmd(f"volume {source} copy reflect {axis}")
    return _new_volume(cubit, before, f"reflected copy of volume {source}")


def _add_block(cubit, block_id: int, name: str, volumes: list[int]) -> None:
    joined = " ".join(str(value) for value in volumes)
    cubit.cmd(f"block {block_id} add volume {joined}")
    cubit.cmd(f'block {block_id} name "{name}"')


def _add_sideset(cubit, sideset_id: int, name: str, surfaces: set[int]) -> None:
    if not surfaces:
        raise RuntimeError(f"Cannot create empty sideset {name!r}")
    joined = " ".join(str(value) for value in sorted(surfaces))
    cubit.cmd(f"sideset {sideset_id} add surface {joined}")
    cubit.cmd(f'sideset {sideset_id} name "{name}"')


def build_reflection_invariant_physical_mesh(
    *,
    iron_size: float,
    air_size: float | None = None,
    gap_size: float = 0.002,
    kelvin_radius: float = 0.22,
) -> dict[str, object]:
    """Mesh positive-z physical volumes and copy their mesh through z=0."""
    import cubit

    if iron_size <= 0.0 or gap_size <= 0.0 or kelvin_radius <= 0.0:
        raise ValueError("mesh sizes and Kelvin radius must be positive")
    if air_size is not None and air_size <= 0.0:
        raise ValueError("air_size must be positive")

    cubit.cmd("unite volume all")
    iron_up = _single(_volume_ids(cubit), "positive-z iron")

    if air_size is None:
        cubit.cmd(f"volume {iron_up} scheme tetmesh")
        cubit.cmd(f"volume {iron_up} size {iron_size:.17g}")
        cubit.cmd(f"mesh volume {iron_up}")
        iron_down = _reflect_meshed_volume(cubit, iron_up)
        cubit.cmd(f"imprint volume {iron_up} {iron_down}")
        cubit.cmd(f"merge volume {iron_up} {iron_down}")
        _add_block(cubit, 1, "iron", [iron_up, iron_down])
        iron_surfaces = set(cubit.get_relatives("volume", iron_up, "surface"))
        iron_surfaces.update(cubit.get_relatives("volume", iron_down, "surface"))
        _add_sideset(cubit, 1, "iron_boundary", iron_surfaces)
        return {
            "iron_up": iron_up,
            "iron_down": iron_down,
            "air_up": None,
            "air_down": None,
        }

    before_sphere = _volume_ids(cubit)
    cubit.cmd(f"create sphere radius {kelvin_radius:.17g}")
    sphere = _new_volume(cubit, before_sphere, "physical air sphere")
    before_cut = _volume_ids(cubit)
    cubit.cmd(f"webcut volume {sphere} with plane zplane")
    sphere_halves = {sphere} | (_volume_ids(cubit) - before_cut)
    if len(sphere_halves) != 2:
        raise RuntimeError(
            f"Expected two physical sphere halves, found {sorted(sphere_halves)}"
        )
    air_up_seed = max(sphere_halves, key=lambda value: cubit.volume(value).centroid()[2])
    air_down_seed = min(sphere_halves, key=lambda value: cubit.volume(value).centroid()[2])
    if cubit.volume(air_up_seed).centroid()[2] <= 0.0:
        raise RuntimeError("Cubit did not produce a positive-z physical hemisphere")
    cubit.cmd(f"delete volume {air_down_seed}")

    cubit.cmd(f"subtract volume {iron_up} from volume {air_up_seed} keep_tool")
    remaining = _volume_ids(cubit)
    if iron_up not in remaining:
        raise RuntimeError("Physical-air subtraction destroyed the iron tool")
    air_up = _single(remaining - {iron_up}, "positive-z physical air")
    cubit.cmd(f"imprint volume {iron_up} {air_up}")
    cubit.cmd(f"merge volume {iron_up} {air_up}")

    gap_half = 0.5 * 0.010
    gap_band = gap_half + 1.0e-7
    gap_surfaces = {
        int(surface)
        for surface in cubit.get_relatives("volume", iron_up, "surface")
        if cubit.surface(surface).area() < 0.002
        and -gap_band < cubit.get_center_point("surface", surface)[2] < gap_band
    }
    if not gap_surfaces:
        raise RuntimeError("No C-type pole surface was selected for 2 mm gap refinement")

    cubit.cmd(f"volume {iron_up} {air_up} scheme tetmesh")
    cubit.cmd(f"volume {iron_up} size {iron_size:.17g}")
    cubit.cmd(f"volume {air_up} size {air_size:.17g}")
    cubit.cmd(
        "surface "
        + " ".join(str(value) for value in sorted(gap_surfaces))
        + f" size {gap_size:.17g}"
    )
    cubit.cmd(f"mesh volume {iron_up} {air_up}")

    iron_down = _reflect_meshed_volume(cubit, iron_up)
    air_down = _reflect_meshed_volume(cubit, air_up)
    physical_volumes = [iron_up, air_up, iron_down, air_down]
    joined = " ".join(str(value) for value in physical_volumes)
    cubit.cmd(f"imprint volume {joined}")
    cubit.cmd(f"merge volume {joined}")

    _add_block(cubit, 1, "iron", [iron_up, iron_down])
    _add_block(cubit, 2, "air", [air_up, air_down])
    iron_surfaces = set(cubit.get_relatives("volume", iron_up, "surface"))
    iron_surfaces.update(cubit.get_relatives("volume", iron_down, "surface"))
    _add_sideset(cubit, 1, "iron_air_interface", iron_surfaces)

    return {
        "iron_up": iron_up,
        "iron_down": iron_down,
        "air_up": air_up,
        "air_down": air_down,
    }
