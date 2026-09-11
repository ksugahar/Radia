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


GAP_HALF = 0.5 * 0.010
# The column is the CAD's `thick` x `width` footprint: the chamfer base around
# the 34 x 24 mm pole tip, i.e. the tip plus 8 mm on every side.
GAP_COLUMN_FOOTPRINT = (0.050, 0.040)
_BOX_TOLERANCE = 1.0e-9


def _box(cubit, volume: int) -> tuple[float, ...]:
    return tuple(float(value) for value in cubit.volume(volume).bounding_box())


def _is_column_slab(cubit, volume: int, bottom: float, top: float) -> bool:
    sx, sy = GAP_COLUMN_FOOTPRINT
    expected = (-0.5 * sx, -0.5 * sy, bottom, 0.5 * sx, 0.5 * sy, top)
    return all(abs(a - b) <= _BOX_TOLERANCE for a, b in zip(_box(cubit, volume), expected))


def _cut_gap_layers(cubit, air: int, half_layers: int) -> tuple[int, list[int]]:
    """Split the positive-z gap air under the pole into ``half_layers`` slabs.

    Sizing the pole faces refines the gap IN PLANE only; a tetrahedron may
    still span the whole half-gap in z.  Parallel slabs force element faces
    at every layer interface, so any line through the gap crosses at least
    one element per slab and never runs further than one slab thickness
    inside a single element.  Every produced volume is identified by its
    bounding box and checked against the intended slab, never by id.
    """
    if half_layers < 1:
        raise ValueError("half_layers must be positive")
    sx, sy = GAP_COLUMN_FOOTPRINT
    before = _volume_ids(cubit)
    cubit.cmd(f"create brick x {sx:.17g} y {sy:.17g} z {GAP_HALF:.17g}")
    tool = _new_volume(cubit, before, "gap column tool")
    cubit.cmd(f"move volume {tool} z {0.5 * GAP_HALF:.17g} include_merged")
    before_cut = _volume_ids(cubit)
    cubit.cmd(f"webcut volume {air} with tool volume {tool}")
    pieces = ({air} | (_volume_ids(cubit) - before_cut)) - {tool}
    if tool in _volume_ids(cubit):
        cubit.cmd(f"delete volume {tool}")
    columns = [v for v in pieces if _is_column_slab(cubit, v, 0.0, GAP_HALF)]
    if len(columns) != 1 or len(pieces) != 2:
        boxes = {v: _box(cubit, v) for v in sorted(pieces)}
        raise RuntimeError(f"gap column webcut produced {boxes}")
    column = columns[0]
    remainder = _single(pieces - {column}, "air outside the gap column")

    thickness = GAP_HALF / half_layers
    layers = [column]
    for index in range(1, half_layers):
        target = max(layers, key=lambda v: _box(cubit, v)[5])
        before_layer = _volume_ids(cubit)
        cubit.cmd(f"webcut volume {target} with plane zplane offset "
                  f"{index * thickness:.17g}")
        layers.append(_new_volume(cubit, before_layer, f"gap layer {index}"))
    layers.sort(key=lambda v: _box(cubit, v)[2])
    for index, volume in enumerate(layers):
        bottom = index * thickness
        top = GAP_HALF if index == half_layers - 1 else bottom + thickness
        if not _is_column_slab(cubit, volume, bottom, top):
            raise RuntimeError(
                f"gap layer {index} is {_box(cubit, volume)}; expected "
                f"z=[{bottom:.9e}, {top:.9e}] under the pole footprint")
        expected_volume = sx * sy * (top - bottom)
        actual_volume = float(cubit.volume(volume).volume())
        if abs(actual_volume - expected_volume) > 1.0e-9 * expected_volume:
            raise RuntimeError(f"gap layer {index} volume {actual_volume:.12e} "
                               f"!= {expected_volume:.12e} m^3")
    return remainder, layers


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
    gap_layers: int = 0,
) -> dict[str, object]:
    """Mesh positive-z physical volumes and copy their mesh through z=0.

    ``gap_layers`` counts slabs through the FULL gap; the positive half is cut
    into ``gap_layers // 2`` of them and the reflection supplies the rest, so
    it must be even.  Zero keeps the historical unlayered gap.
    """
    import cubit

    if iron_size <= 0.0 or gap_size <= 0.0 or kelvin_radius <= 0.0:
        raise ValueError("mesh sizes and Kelvin radius must be positive")
    if air_size is not None and air_size <= 0.0:
        raise ValueError("air_size must be positive")
    if isinstance(gap_layers, bool) or int(gap_layers) != gap_layers \
            or gap_layers < 0 or gap_layers % 2:
        raise ValueError("gap_layers must be a non-negative even integer")
    if gap_layers and air_size is None:
        raise ValueError("gap_layers needs the physical air domain")

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
    layers_up: list[int] = []
    if gap_layers:
        air_up, layers_up = _cut_gap_layers(cubit, air_up, gap_layers // 2)
    positive = " ".join(str(value) for value in [iron_up, air_up, *layers_up])
    cubit.cmd(f"imprint volume {positive}")
    cubit.cmd(f"merge volume {positive}")

    gap_half = GAP_HALF
    gap_band = gap_half + 1.0e-7
    gap_surfaces = {
        int(surface)
        for surface in cubit.get_relatives("volume", iron_up, "surface")
        if cubit.surface(surface).area() < 0.002
        and -gap_band < cubit.get_center_point("surface", surface)[2] < gap_band
    }
    if not gap_surfaces:
        raise RuntimeError("No C-type pole surface was selected for 2 mm gap refinement")

    cubit.cmd(f"volume {positive} scheme tetmesh")
    cubit.cmd(f"volume {iron_up} size {iron_size:.17g}")
    cubit.cmd(f"volume {air_up} size {air_size:.17g}")
    if layers_up:
        cubit.cmd("volume " + " ".join(str(value) for value in layers_up)
                  + f" size {gap_size:.17g}")
    cubit.cmd(
        "surface "
        + " ".join(str(value) for value in sorted(gap_surfaces))
        + f" size {gap_size:.17g}"
    )
    cubit.cmd(f"mesh volume {positive}")

    iron_down = _reflect_meshed_volume(cubit, iron_up)
    air_down = _reflect_meshed_volume(cubit, air_up)
    layers_down = [_reflect_meshed_volume(cubit, value) for value in layers_up]
    physical_volumes = [iron_up, air_up, *layers_up, iron_down, air_down, *layers_down]
    joined = " ".join(str(value) for value in physical_volumes)
    cubit.cmd(f"imprint volume {joined}")
    cubit.cmd(f"merge volume {joined}")

    _add_block(cubit, 1, "iron", [iron_up, iron_down])
    _add_block(cubit, 2, "air", [air_up, *layers_up, air_down, *layers_down])
    iron_surfaces = set(cubit.get_relatives("volume", iron_up, "surface"))
    iron_surfaces.update(cubit.get_relatives("volume", iron_down, "surface"))
    _add_sideset(cubit, 1, "iron_air_interface", iron_surfaces)

    return {
        "iron_up": iron_up,
        "iron_down": iron_down,
        "air_up": air_up,
        "air_down": air_down,
        "gap_layers_up": layers_up,
        "gap_layers_down": layers_down,
    }
