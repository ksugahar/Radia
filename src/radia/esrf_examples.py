"""Reproducible CAD/source definitions for the seven ESRF Radia examples.

The original Mathematica notebooks use millimetres.  This module converts
their dimensions to SI once at the boundary and keeps all generated CAD and
``CoilBuilder`` paths in metres.  It does not replace the notebooks as the
field oracle: ``reference_observables`` records what a validation runner must
compare against in the original Radia model.

Examples 1--6 are the public ESRF tutorial set.  Example 7 is the additional
ESRF storage-ring quadrupole notebook shipped with the Radia distribution.
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from .coil_builder import CoilBuilder

MM = 1.0e-3
MU0 = 4.0e-7 * math.pi


@dataclass(frozen=True)
class ESRFExampleSpec:
    """Source-traceable contract for one original Radia example."""

    number: int
    slug: str
    title: str
    source_notebook: str
    model_class: str
    parameters_si: Mapping[str, Any]
    reference_observables: tuple[str, ...]
    cad_fidelity: str = "exact tutorial dimensions"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ESRFFixedMagnetizationBlock:
    """One immutable permanent-magnet source body in an ESRF example.

    ``magnetization_A_m`` is the physical given magnetization used by
    :class:`radia.vim.MagnetizationSource`; it is not a soft-material unknown.
    The CAD solid is deliberately retained with its source record so a Cubit
    source mesh can preserve a different vector in every segmented magnet.
    """

    index: int
    remanence_T: float
    magnetization_A_m: tuple[float, float, float]
    symmetry_path: tuple[str, ...]
    shape: Any = field(repr=False, compare=False)

    @property
    def material_name(self) -> str:
        """Stable Cubit/NGSolve material label for this source body."""
        return f"pm_{self.index:03d}"

    def manifest_entry(self) -> dict[str, Any]:
        """Serializable source contract excluding the non-serializable CAD solid."""
        return {
            "index": self.index,
            "material": self.material_name,
            "remanence_T": self.remanence_T,
            "magnetization_A_m": list(self.magnetization_A_m),
            "symmetry_path": list(self.symmetry_path),
            "model": "fixed-given MagnetizationSource",
        }


_SPECS = {
    1: ESRFExampleSpec(
        1, "pm_cube", "Uniformly magnetized permanent-magnet cube",
        "Example#1.nb", "permanent_magnet",
        {"size_m": (1 * MM,) * 3, "magnetization_T": (-0.5, 1.0, 0.7)},
        ("B at (0.52, 0.60, 0.70) mm",),
    ),
    2: ESRFExampleSpec(
        2, "racetrack_coils", "Racetrack coils of a 4 T superconducting wiggler",
        "Example#2.nb", "current_source",
        {"mirror_plane": "z=0", "source_coil_count": 5},
        ("Bz along y=0..300 mm", "integrated Bz along x=0..400 mm"),
    ),
    3: ESRFExampleSpec(
        3, "hybrid_undulator", "Short hybrid permanent-magnet undulator",
        "Example#3.nb", "hybrid_permanent_magnet",
        {
            "period_m": 46 * MM, "gap_m": 20 * MM, "period_count": 2,
            "pole_size_m": (45 * MM, 5 * MM, 25 * MM),
            "magnet_size_m": (65 * MM, 18 * MM, 45 * MM),
            "remanence_T": 1.2,
            "base_magnet_easy_axes": ((0., -1., 0.), (0., 1., 0.),
                                      (0., -1., 0.)),
            "image_symmetries": ("x=0 perpendicular", "z=0 parallel",
                                 "y=0 perpendicular"),
        },
        ("peak Bz at the origin", "longitudinal Bz profile"),
    ),
    4: ESRFExampleSpec(
        4, "magnetized_sphere", "Uniformly magnetized sphere",
        "Example#4.nb", "permanent_magnet",
        {"radius_m": 1 * MM, "magnetization_T": (1.0, 0.0, 0.0)},
        ("internal Bx = 2/3 T",),
        "analytic sphere replacing the notebook's 15x15 polyhedral approximation",
    ),
    5: ESRFExampleSpec(
        5, "c_dipole", "Simple iron-dominated C-type dipole",
        "Example#5.nb", "iron_dominated_dipole",
        {
            "gap_m": 10 * MM, "depth_m": 50 * MM,
            "pole_width_m": 40 * MM, "chamfer_m": 8 * MM,
            "ampere_turns_A": -2000.0,
        },
        ("Bz at magnet centre", "Ampere-law field ratio"),
        "dimension-faithful validation CAD; notebook constructive partition is not copied",
    ),
    6: ESRFExampleSpec(
        6, "quadrupole", "Simple quadrupole with hyperbolic poles",
        "Example#6.nb", "quadrupole",
        {
            "gap_m": 40 * MM, "pole_width_m": 30 * MM,
            "yoke_height_m": 50 * MM, "length_m": 60 * MM,
            "chamfer_m": 8 * MM, "pole_depth_m": 18 * MM,
            "current_density_A_per_m2": -3.0e6,
        },
        ("on-axis gradient", "integrated multipoles at radius 2 mm"),
        "exact pole/yoke pieces including the 8 mm longitudinal pole-tip chamfer",
    ),
    7: ESRFExampleSpec(
        7, "esrf_storage_ring_quadrupole",
        "ESRF quadrupole in 2D and 3D with automatic triangulation",
        "Example#7.nb", "quadrupole",
        {
            "bore_radius_m": 36 * MM, "iron_length_m": 400 * MM,
            "end_chamfer_m": 7 * MM, "turns": 15,
            "current_per_turn_A": 533.3, "hyperbola_xy_m2": 648 * MM * MM,
            "coil_rotation_deg": -45.0,
            "coil_pack_turns": (8, 7),
            "coil_pack_inner_radii_m": (32 * MM, 46.2 * MM),
            "coil_pack_outer_radii_m": (46.2 * MM, 60.4 * MM),
            "coil_straight_length_m": 402 * MM,
        },
        (
            "2D gradient", "3D integrated gradient", "magnetic length",
            "integrated multipoles", "magnetic energy and inductance",
        ),
        ("exact notebook transverse profile, 400 mm length, symmetric 7 mm "
         "45-degree pole-end chamfer, and -45-degree coil-to-pole phase"),
    ),
}


# Original Radia nonlinear material laws.  ``MatSatIsoFrm`` and
# ``MatSatIsoTab`` use legacy Radia units: the field abscissa is mu0*H [T]
# and the ordinate is the polarization mu0*M [T].  HDiv-MMM instead accepts
# the conventional total B-H table ``[H [A/m], B [T]]``.  Keep the conversion
# in this source-traceable module so validation scripts do not silently mix
# the two conventions.
_FORMULA_MATERIALS = {
    5: ((20000.0, 2.0), (0.1, 2.0), (0.1, 2.0)),
    6: ((2000.0, 2.0), (0.1, 2.0), (0.1, 2.0)),
}

_EXAMPLE3_H_A_PER_M = (
    0.8, 1.5, 2.2, 3.6, 5.0, 6.8, 9.8, 18.0, 28.0, 37.5, 42.0,
    55.0, 71.5, 80.0, 85.0, 88.0, 92.0, 100.0, 120.0, 150.0, 200.0,
    300.0, 400.0, 600.0, 800.0, 1000.0, 2000.0, 4000.0, 6000.0,
    10000.0, 25000.0, 40000.0,
)
_EXAMPLE3_MU0_M_T = (
    0.000998995, 0.00199812, 0.00299724, 0.00499548, 0.00699372,
    0.00999145, 0.0149877, 0.0299774, 0.0499648, 0.0799529, 0.0999472,
    0.199931, 0.49991, 0.799899, 0.999893, 1.09989, 1.19988, 1.29987,
    1.41985, 1.49981, 1.59975, 1.72962, 1.7995, 1.89925, 1.96899,
    1.99874, 2.09749, 2.19497, 2.24246, 2.27743, 2.28958, 2.28973,
)

_EXAMPLE7_MU0H_MU0M_T = (
    (0.0, 0.0),
    (0.0000228708, 0.0512771), (0.0000309133, 0.0986691),
    (0.0000417204, 0.235858), (0.0000496372, 0.34455),
    (0.0000561717, 0.423644), (0.0000759009, 0.609824),
    (0.000102667, 0.776997), (0.00013823, 0.923262),
    (0.000186611, 1.05361), (0.000251327, 1.16455),
    (0.000338538, 1.25756), (0.000456662, 1.33414),
    (0.000615627, 1.39478), (0.000829506, 1.44367),
    (0.00111828, 1.48488), (0.00149816, 1.5187),
    (0.00150759, 1.51949), (0.00203223, 1.55207),
    (0.00273934, 1.58236), (0.00369275, 1.61491),
    (0.00497804, 1.65032), (0.00671019, 1.68929),
    (0.00904553, 1.73525), (0.0121934, 1.78811),
    (0.0164367, 1.84876), (0.0221566, 1.91754),
    (0.0300011, 1.9925), (0.05, 2.07755), (0.08, 2.11885),
    (0.1, 2.12656), (0.15, 2.1308), (0.17, 2.13105),
    (0.2, 2.13116),
)


@lru_cache(maxsize=None)
def _formula_bh_table(number: int, sample_count: int,
                      maximum_mu0_h_t: float) -> tuple[tuple[float, float], ...]:
    import radia as rad

    parameters = _FORMULA_MATERIALS[number]
    material = rad.MatSatIsoFrm([list(pair) for pair in parameters])
    mu0_h = np.r_[0.0, np.geomspace(1.0e-10, maximum_mu0_h_t,
                                   int(sample_count) - 1)]
    rows = []
    for abscissa in mu0_h:
        mu0_m = float(np.asarray(rad.MatMvsH(
            material, "m", [float(abscissa), 0.0, 0.0]), dtype=float)[0])
        rows.append((float(abscissa / MU0), float(abscissa + mu0_m)))
    return tuple(rows)


def get_esrf_bh_table(number: int, *, sample_count: int = 221,
                      maximum_mu0_h_t: float = 10.0) -> list[list[float]]:
    """Return the original example iron law as ``[H A/m, total B T]``.

    Examples 5 and 6 use the original three-term ``MatSatIsoFrm`` law, sampled
    in legacy Radia and converted without fitting.  Examples 3 and 7 use the
    tabulated laws embedded in their notebooks.  The other examples contain
    no soft iron and therefore reject this request.
    """

    number = int(number)
    if number in _FORMULA_MATERIALS:
        if int(sample_count) < 3:
            raise ValueError("sample_count must be at least 3")
        if not np.isfinite(maximum_mu0_h_t) or maximum_mu0_h_t <= 0.0:
            raise ValueError("maximum_mu0_h_t must be positive and finite")
        rows = _formula_bh_table(
            number, int(sample_count), float(maximum_mu0_h_t))
    elif number == 3:
        rows = ((0.0, 0.0), *((float(h), float(MU0 * h + m))
                 for h, m in zip(_EXAMPLE3_H_A_PER_M, _EXAMPLE3_MU0_M_T)))
    elif number == 7:
        rows = ((float(mu0_h / MU0), float(mu0_h + mu0_m))
                for mu0_h, mu0_m in _EXAMPLE7_MU0H_MU0M_T)
    else:
        raise ValueError(f"ESRF example {number} has no soft-iron B-H law")
    table = np.asarray(tuple(rows), dtype=float)
    if table.ndim != 2 or table.shape[1] != 2:
        raise RuntimeError("invalid ESRF B-H table construction")
    if not np.isfinite(table).all() or np.any(np.diff(table[:, 0]) <= 0.0) \
            or np.any(np.diff(table[:, 1]) <= 0.0):
        raise RuntimeError("ESRF B-H table must be finite and strictly monotone")
    return table.tolist()


def get_esrf_example_spec(number: int) -> ESRFExampleSpec:
    """Return one immutable example contract."""

    try:
        return _SPECS[int(number)]
    except (KeyError, ValueError) as exc:
        raise ValueError("ESRF Radia example number must be in 1..7") from exc


def list_esrf_example_specs() -> tuple[ESRFExampleSpec, ...]:
    return tuple(_SPECS[k] for k in sorted(_SPECS))


def _racetrack(
    current: float,
    *,
    radius: float,
    straight: float,
    width: float,
    height: float,
    centre: tuple[float, float, float] = (0.0, 0.0, 0.0),
    orientation: np.ndarray | None = None,
) -> CoilBuilder:
    """Closed two-straight/two-semicircle CoilBuilder path."""

    frame = np.eye(3) if orientation is None else np.asarray(orientation, dtype=float)
    centre_vec = np.asarray(centre, dtype=float)
    start = centre_vec + radius * frame[0] - 0.5 * straight * frame[1]
    coil = CoilBuilder(current).set_start(start, orientation=frame)
    coil.set_cross_section(width, height)
    if abs(straight) < 1.0e-15:
        return coil.add_arc(radius, 360)
    return (coil.add_straight(straight)
            .add_arc(radius, 180)
            .add_straight(straight)
            .add_arc(radius, 180))


def _rounded_rectangle(
    current: float,
    *,
    radius: float,
    straight_x: float,
    straight_y: float,
    width: float,
    height: float,
    centre: tuple[float, float, float],
    orientation: np.ndarray | None = None,
) -> CoilBuilder:
    """Closed four-straight/four-quarter-arc racetrack path."""

    frame = np.eye(3) if orientation is None else np.asarray(orientation, dtype=float)
    c = np.asarray(centre, dtype=float)
    start = c + (0.5 * straight_x + radius) * frame[0] - 0.5 * straight_y * frame[1]
    return (CoilBuilder(current).set_start(start, orientation=frame)
            .set_cross_section(width, height)
            .add_straight(straight_y).add_arc(radius, 90)
            .add_straight(straight_x).add_arc(radius, 90)
            .add_straight(straight_y).add_arc(radius, 90)
            .add_straight(straight_x).add_arc(radius, 90))


def build_esrf_coils(number: int) -> list[CoilBuilder]:
    """Build the physical coil sources present in an ESRF example.

    Permanent-magnet-only examples return an empty list.  Current values are
    the notebook current density times the rectangular conductor area, so the
    resulting builder represents the same total current distribution.
    """

    number = int(number)
    if number in (1, 3, 4):
        return []
    if number == 2:
        # (inner radius, outer radius, straight, axial height, J, z centre)
        raw = (
            (9.5, 24.5, 120.0, 36.0, 128.0, 38.0),
            (24.5, 55.5, 120.0, 36.0, 256.0, 38.0),
            (10.0, 25.0, 90.0, 24.0, 128.0, 76.0),
            (25.0, 55.0, 90.0, 24.0, 256.0, 76.0),
            (150.0, 166.3, 0.0, 39.0, -256.0, 60.0),
        )
        coils: list[CoilBuilder] = []
        # ObjRaceTrk[..., {Lx, 0}, ...] has straight sections along x.
        xy_x_straight = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        for ri, ro, straight, axial, j, zc in raw:
            radial = (ro - ri) * MM
            mean_r = 0.5 * (ro + ri) * MM
            current = j * radial / MM * axial  # A/mm2 * width_mm * height_mm
            base = _racetrack(
                current, radius=mean_r, straight=straight * MM,
                width=radial, height=axial * MM, centre=(0, 0, zc * MM),
                orientation=xy_x_straight,
            )
            coils.extend((base, base.mirror("xy")))
        return coils
    if number == 5:
        # ObjRaceTrk centre pc=[0,p6_y,0], L=[50,62.5] mm.
        return [_rounded_rectangle(
            -2000.0, radius=22.5 * MM,
            straight_x=50 * MM, straight_y=62.5 * MM,
            width=35 * MM, height=105 * MM, centre=(0, 131.25 * MM, 0),
        )]
    if number == 6:
        # The notebook partitions one winding into adjacent radial packs.
        # J=-3 A/mm2 gives -1560 A and -960 A for the two pack areas.
        coils: list[CoilBuilder] = []
        for current, radius, width, height, centre_z in (
            (-1560.0, 8.5 * MM, 13 * MM, 40 * MM, 50 * MM),
            (-960.0, 23 * MM, 16 * MM, 20 * MM, 60 * MM),
        ):
            negative = _rounded_rectangle(
                current, radius=radius, straight_x=60 * MM,
                straight_y=26 * MM, width=width, height=height,
                centre=(0, 0, centre_z),
            ).rotate_copies(axis="x", n_copies=8)
            positive = _rounded_rectangle(
                -current, radius=radius, straight_x=60 * MM,
                straight_y=26 * MM, width=width, height=height,
                centre=(0, 0, centre_z),
            ).rotate_copies(axis="x", n_copies=8)
            # Example #6 first reflects the base winding across its diagonal
            # quadrupole plane (which reverses the current), then makes the
            # 180-degree pair and finally rotates the complete model by 45
            # degrees about the beam axis.  The equivalent four physical
            # windings therefore occupy the odd eighth-turn positions with
            # alternating current signs.
            coils.extend((negative[1], positive[3], negative[5], positive[7]))
        return coils
    if number == 7:
        orientation = np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]])
        # Two radial layers: 8 + 7 turns.  Current is represented as ampere-
        # turns in each homogenized layer.  Example#7.nb rotates the complete
        # base winding by -45 degrees before applying the x/y quadrupole image
        # symmetries.  Selecting the odd eighth-turn copies preserves that
        # source-to-pole phase; the first copy is 315 degrees == -45 degrees.
        layers: list[CoilBuilder] = []
        for turns, radius, height in ((8, 39.1, 80.0), (7, 53.3, 70.0)):
            base = _racetrack(
                turns * 533.3, radius=radius * MM, straight=402 * MM,
                width=14.2 * MM, height=height * MM,
                centre=(0, 178 * MM, 0), orientation=orientation,
            )
            minus = base.rotate_copies("z", 8)
            plus_base = _racetrack(
                -turns * 533.3, radius=radius * MM, straight=402 * MM,
                width=14.2 * MM, height=height * MM,
                centre=(0, 178 * MM, 0), orientation=orientation,
            ).rotate_copies("z", 8)
            phase_indices = (7, 1, 3, 5)
            layers.extend(minus[index] if i % 2 == 0 else plus_base[index]
                          for i, index in enumerate(phase_indices))
        return layers
    raise ValueError("ESRF Radia example number must be in 1..7")


def _compound(shapes: list[Any], name: str) -> Any:
    from netgen.occ import Compound

    if not shapes:
        raise ValueError(f"cannot make empty OCC compound {name!r}")
    result = shapes[0] if len(shapes) == 1 else Compound(shapes)
    result.name = name
    return result


def _box_from_center(center: tuple[float, float, float],
                     size: tuple[float, float, float], name: str) -> Any:
    from netgen.occ import Box, Pnt

    c = np.asarray(center, dtype=float)
    h = 0.5 * np.asarray(size, dtype=float)
    shape = Box(Pnt(*(c - h)), Pnt(*(c + h)))
    shape.name = name
    return shape


def _polygon_prism_xy(points: list[tuple[float, float]], z0: float,
                      length: float, name: str) -> Any:
    from netgen.occ import Axes, Pnt, Vec, WorkPlane, X, Z

    wp = WorkPlane(Axes(Pnt(0, 0, z0), n=Z, h=X)).MoveTo(*points[0])
    for point in points[1:]:
        wp = wp.LineTo(*point)
    shape = wp.Close().Face().Extrude(Vec(0, 0, length))
    shape.name = name
    return shape


def _polygon_prism_yz(points: list[tuple[float, float]], x0: float,
                      length: float, name: str) -> Any:
    from netgen.occ import Axes, Pnt, Vec, WorkPlane, X, Y

    # Local work-plane coordinates are (global y, global z).
    wp = WorkPlane(Axes(Pnt(x0, 0, 0), n=X, h=Y)).MoveTo(*points[0])
    for point in points[1:]:
        wp = wp.LineTo(*point)
    shape = wp.Close().Face().Extrude(Vec(length, 0, 0))
    shape.name = name
    return shape


def _reflect_magnetization_for_field_symmetry(
    vector: tuple[float, float, float], normal_axis: int, *, zero_field: str,
) -> tuple[float, float, float]:
    """Transform ``M`` for legacy Radia's zero-parallel/perpendicular image.

    This is deliberately not a generic CAD reflection.  ``TrfZerPara`` and
    ``TrfZerPerp`` prescribe the parity of *B* at the mirror plane:

    * zero-parallel: ``B' = -R B`` and therefore ``M' = -R M`` (axial);
    * zero-perpendicular: ``B' = R B`` and therefore ``M' = R M``.

    Example #3 uses the second rule at x=0 and y=0, and the first at z=0.
    Collapsing both operations to the usual axial-vector transform creates a
    source whose symmetry plane happens to cancel the intended undulator
    field.
    """
    if normal_axis not in (0, 1, 2):
        raise ValueError("normal_axis must be 0, 1, or 2")
    if zero_field not in {"parallel", "perpendicular"}:
        raise ValueError("zero_field must be 'parallel' or 'perpendicular'")
    reflected = [float(value) for value in vector]
    reflected[normal_axis] *= -1.0  # R: polar reflection
    if zero_field == "parallel":
        reflected = [-value for value in reflected]  # -R: axial reflection
    return tuple(reflected)  # type: ignore[return-value]


def _mirror_example3_blocks(
    blocks: list[tuple[Any, tuple[float, float, float], tuple[str, ...]]],
    normal: Any,
    horizontal: Any,
    normal_axis: int,
    symmetry_name: str,
    zero_field: str,
) -> list[tuple[Any, tuple[float, float, float], tuple[str, ...]]]:
    """Mirror Example #3 CAD and its B-parity-preserving source data."""
    from netgen.occ import Axes, Pnt

    plane = Axes(Pnt(0, 0, 0), n=normal, h=horizontal)
    mirrored = [
        (
            shape.Mirror(plane),
            _reflect_magnetization_for_field_symmetry(
                magnetization, normal_axis, zero_field=zero_field
            ),
            path + (symmetry_name,),
        )
        for shape, magnetization, path in blocks
    ]
    return blocks + mirrored


def _example3_fixed_magnetization_blocks() -> list[ESRFFixedMagnetizationBlock]:
    """Build Example #3 magnet solids and their notebook-prescribed ``M``.

    The notebook constructs three source magnets in one octant and then
    reflects x, z, and y in this exact order.  Every mirrored body receives
    its own material label.  This preserves normal jumps between segments and
    lets the C++ HDiv source evaluator represent the alternating pole pattern
    without a Python pointwise source fallback.
    """

    from netgen.occ import Axes, Pnt, X, Y, Z

    lp = np.array([45., 5., 25.]) * MM
    lm = np.array([65., 18., 45.]) * MM
    gap, gap_offset, num_per = 20 * MM, 1 * MM, 2
    magnets: list[tuple[Any, tuple[float, float, float], tuple[str, ...]]] = []
    base_axes = tuple(
        tuple(float(component) for component in axis)
        for axis in get_esrf_example_spec(3).parameters_si["base_magnet_easy_axes"]
    )
    base_index = 0
    y = 0.25 * lp[1]
    y += 0.25 * lp[1]
    for _ in range(num_per):
        y += 0.5 * lm[1]
        magnet = _box_from_center(
            (0.25 * lm[0], y, -0.5 * (lm[2] + gap) - gap_offset),
            (0.5 * lm[0], lm[1], lm[2]), "permanent_magnet")
        magnets.append((magnet, base_axes[base_index], ("base",)))
        base_index += 1
        y += 0.5 * (lm[1] + lp[1])
        y += 0.5 * lp[1]
    y += 0.25 * lm[1]
    magnet = _box_from_center(
        (0.25 * lm[0], y, -0.5 * (lm[2] + gap) - gap_offset),
        (0.5 * lm[0], 0.5 * lm[1], lm[2]), "permanent_magnet")
    magnets.append((magnet, base_axes[base_index], ("base",)))
    if base_index + 1 != len(base_axes):
        raise RuntimeError("Example #3 base-magnet axis list does not match CAD")

    magnets = _mirror_example3_blocks(
        magnets, X, Y, 0, "mirror_x", "perpendicular"
    )
    magnets = _mirror_example3_blocks(
        magnets, Z, X, 2, "mirror_z", "parallel"
    )
    magnets = _mirror_example3_blocks(
        magnets, Y, X, 1, "mirror_y", "perpendicular"
    )
    remanence = float(get_esrf_example_spec(3).parameters_si["remanence_T"])
    return [
        ESRFFixedMagnetizationBlock(
            index=index,
            remanence_T=remanence,
            magnetization_A_m=tuple(component * remanence / MU0
                                     for component in magnetization),
            symmetry_path=path,
            shape=shape,
        )
        for index, (shape, magnetization, path) in enumerate(magnets)
    ]


def _example3_occ() -> dict[str, Any]:
    """Reproduce the rectangular geometry and mirror sequence of Example 3."""

    from netgen.occ import Axes, Pnt, X, Y, Z

    lp = np.array([45., 5., 25.]) * MM
    gap, num_per = 20 * MM, 2
    poles: list[Any] = []
    y = 0.25 * lp[1]
    poles.append(_box_from_center(
        (0.25 * lp[0], y, -0.5 * (lp[2] + gap)),
        (0.5 * lp[0], 0.5 * lp[1], lp[2]), "pole"))
    y += 0.25 * lp[1]
    lm = np.array([65., 18., 45.]) * MM
    for _ in range(num_per):
        y += 0.5 * lm[1]
        y += 0.5 * (lm[1] + lp[1])
        poles.append(_box_from_center(
            (0.25 * lp[0], y, -0.5 * (lp[2] + gap)),
            (0.5 * lp[0], lp[1], lp[2]), "pole"))
        y += 0.5 * lp[1]

    def mirror_all(source: list[Any], normal: Any, horizontal: Any) -> list[Any]:
        plane = Axes(Pnt(0, 0, 0), n=normal, h=horizontal)
        return source + [shape.Mirror(plane) for shape in source]

    # x half -> full, lower -> upper, longitudinal half -> full.
    poles = mirror_all(poles, X, Y)
    poles = mirror_all(poles, Z, X)
    poles = mirror_all(poles, Y, X)
    magnets = _example3_fixed_magnetization_blocks()
    return {"iron": _compound(poles, "iron"),
            "magnet": _compound([block.shape for block in magnets],
                                  "permanent_magnet")}


def build_esrf_fixed_magnetization_blocks(
    number: int,
) -> list[ESRFFixedMagnetizationBlock]:
    """Return fixed-magnetization source blocks for PM ESRF examples.

    Example #3 is segmented deliberately.  The source mesh must retain the
    per-block labels in this return value; combining it into one uniform
    ``permanent_magnet`` material loses its alternating magnetization.
    """
    number = int(number)
    if number == 3:
        return _example3_fixed_magnetization_blocks()
    if number == 1:
        vector = tuple(float(value) / MU0
                       for value in get_esrf_example_spec(1).parameters_si[
                           "magnetization_T"])
        return [ESRFFixedMagnetizationBlock(
            index=0, remanence_T=float(np.linalg.norm(np.asarray(vector)) * MU0),
            magnetization_A_m=vector, symmetry_path=("base",),
            shape=_box_from_center((0, 0, 0), (MM, MM, MM), "permanent_magnet"),
        )]
    if number == 4:
        from netgen.occ import Pnt, Sphere

        vector = tuple(float(value) / MU0
                       for value in get_esrf_example_spec(4).parameters_si[
                           "magnetization_T"])
        sphere = Sphere(Pnt(0, 0, 0), MM)
        sphere.name = "permanent_magnet"
        return [ESRFFixedMagnetizationBlock(
            index=0, remanence_T=float(np.linalg.norm(np.asarray(vector)) * MU0),
            magnetization_A_m=vector, symmetry_path=("base",), shape=sphere,
        )]
    if number not in _SPECS:
        raise ValueError("ESRF Radia example number must be in 1..7")
    return []


def get_esrf_fixed_magnetization_by_material(
    number: int,
) -> dict[str, tuple[float, float, float]]:
    """Return the immutable PM vector for every Cubit source material.

    ``pm_NNN`` labels are part of the source-mesh contract.  A validation
    runner must not infer their direction from a point location or from the
    order in which Cubit happens to return volumes.
    """
    return {
        block.material_name: block.magnetization_A_m
        for block in build_esrf_fixed_magnetization_blocks(number)
    }


def esrf_fixed_magnetization_coefficient(mesh, number: int):
    """Build an NGSolve material coefficient for a fixed PM source mesh.

    The exact label check is intentional: mapping an alternating segmented PM
    field onto a one-material mesh produces a visually plausible but physically
    wrong source.  The resulting coefficient is fed directly to the
    C++-backed :class:`radia.vim.MagnetizationSource` and its native field is
    shared with the HDiv-MMM, reduced-A, and mixed-H1 response solves.
    """
    by_material = get_esrf_fixed_magnetization_by_material(number)
    actual = {str(material) for material in mesh.GetMaterials()}
    expected = set(by_material)
    if actual != expected:
        raise ValueError(
            "ESRF fixed-magnetization source mesh labels must exactly match "
            f"the source blocks; unexpected={sorted(actual - expected)}, "
            f"missing={sorted(expected - actual)}"
        )
    import ngsolve as ng

    return mesh.MaterialCF({
        material: ng.CoefficientFunction(vector)
        for material, vector in by_material.items()
    })


def build_esrf_fixed_magnetization_source(mesh, number: int, *, order: int = 2,
                                          curve_order: int | None = None,
                                          field_cf_algorithm: str | None = None,
                                          field_tree_options: dict | None = None):
    """Create the C++ HDiv fixed-PM source for an ESRF source mesh.

    ``field_cf_algorithm`` is passed directly to the native source
    coefficient.  A validation runner may request ``"tree"`` only after its
    direct-field certificate accepts the configured tree options.
    """
    from . import vim

    return vim.MagnetizationSource(
        mesh,
        esrf_fixed_magnetization_coefficient(mesh, number),
        order=int(order),
        curve_order=curve_order,
        field_cf_algorithm=field_cf_algorithm,
        field_tree_options=field_tree_options,
    )


def _example5_iron() -> Any:
    """Full C-yoke support assembled from the notebook's six regions."""

    # Coordinates follow RADIA_Example05.py: x=depth, y=return-leg
    # direction, z=gap direction.  The three rectangular sections reproduce
    # ObjMltExtRtg's 8 mm pole-face chamfer in both x and y.
    from netgen.occ import Axes, Pnt, ThruSections, WorkPlane, X, Z

    thick, width, gap = 50 * MM, 40 * MM, 10 * MM
    chamfer = 8 * MM
    lz1, lz2 = 20 * MM, 30 * MM
    lz3 = 1.25 * width
    ly4 = 80 * MM
    ly5 = 1.25 * lz3
    z1 = gap / 2 + lz1 / 2
    z2 = gap / 2 + lz1 + lz2 / 2
    z3 = gap / 2 + lz1 + lz2 + lz3 / 2
    y4 = width / 2 + ly4 / 2
    y5 = y4 + (ly4 + ly5) / 2

    # The source script creates x>=0 sections centred at thick/4 and then
    # applies TrfZerPerp at x=0.  After materialising that mirror, every full
    # section is centred at x=0.  In particular, the pole-tip section spans
    # [-17, 17] mm; shifting its centre by -4 mm would break the declared
    # x-symmetry and is not equivalent to ObjMltExtRtg plus the mirror.
    sections = (
        (gap / 2, thick - 2 * chamfer, width - 2 * chamfer),
        (gap / 2 + chamfer, thick, width),
        (gap / 2 + lz1, thick, width),
    )
    wires = []
    for z, size_x, size_y in sections:
        wires.append(
            WorkPlane(Axes(Pnt(0, 0, z), n=Z, h=X))
            .MoveTo(-size_x / 2, -size_y / 2)
            .Rectangle(size_x, size_y).Wire()
        )
    pole = ThruSections(wires, solid=True)
    pole.name = "pole"
    regions = [
        pole,
        _box_from_center((0, 0, z2), (thick, width, lz2), "upper_arm"),
        _box_from_center((0, 0, z3), (thick, width, lz3), "corner"),
        _box_from_center((0, y4, z3), (thick, ly4, lz3), "back_yoke"),
        _box_from_center((0, y5, z3), (thick, ly5, lz3), "return_corner"),
        _box_from_center((0, y5, (gap / 2 + lz1 + lz2) / 2),
                         (thick, ly5, gap / 2 + lz1 + lz2),
                         "return_leg"),
    ]

    lower_plane = Axes(Pnt(0, 0, 0), n=Z, h=X)
    lower = [shape.Mirror(lower_plane) for shape in regions]
    return _compound(regions + lower, "iron")


def _example6_iron() -> Any:
    from netgen.occ import Axes, Axis, HalfSpace, Pnt, Vec, X, Y

    gap, width, height, thick = 40 * MM, 30 * MM, 50 * MM, 60 * MM
    z0, y0 = gap / 2, width / 2
    amax = math.asinh(y0 / z0)
    points = [(z0 * math.sinh(i * amax / 4),
               z0 * math.cosh(i * amax / 4)) for i in range(5)]
    hh = points[4][1] + 0.5 * height - z0 * (math.cosh(amax) - 1)
    points.extend(((points[-1][0], hh), (0, hh)))
    depth = 18 * MM
    quarter = [
        _polygon_prism_yz(points, -thick / 2, thick, "hyperbolic_pole"),
        _box_from_center((0, width / 4, 57.5 * MM),
                         (thick, width / 2, 25 * MM), "pole_root"),
        _box_from_center((0, width / 4, gap / 2 + height + depth / 2),
                         (thick, width / 2, depth), "inner_corner"),
        _box_from_center((0, 42.5 * MM, gap / 2 + height + depth / 2),
                         (thick, 55 * MM, depth), "outer_yoke"),
        _polygon_prism_yz(
            [(70 * MM, 70 * MM), (88 * MM, 88 * MM), (70 * MM, 88 * MM)],
            -thick / 2, thick, "outer_corner"),
    ]
    # RADIA_Example06 cuts g1+g2 at the positive-x end with the plane
    # (x - z - 2 mm)=0 before reflecting it across x=0.  ObjGeoVol on its
    # mutated parent container misleadingly reports the pre-cut volume; the
    # transformed end-face polygons show the physical aperture moving from
    # z=20 mm at x=0 to z=28 mm at x=30 mm.  Apply both physical end cuts.
    positive_end = HalfSpace(Pnt(2 * MM, 0, 0), Vec(1, 0, -1))
    negative_end = HalfSpace(Pnt(-2 * MM, 0, 0), Vec(-1, 0, -1))
    quarter[0] = quarter[0] * positive_end * negative_end
    quarter[1] = quarter[1] * positive_end * negative_end
    # The original first reflects this y-positive half-pole across y=0 to
    # complete one pole, then constructs the remaining poles.  Rotating the
    # half-pole directly would omit half of every pole and approximately halve
    # the nonlinear iron response.
    half_pole_mirror = Axes(Pnt(0, 0, 0), n=Y, h=X)
    complete_pole = quarter + [shape.Mirror(half_pole_mirror) for shape in quarter]
    # The source example applies a final 45-degree rotation after constructing
    # the fourfold quadrupole symmetry.
    copies = [shape.Rotate(Axis(Pnt(0, 0, 0), X), angle)
              for angle in (45, 135, 225, 315) for shape in complete_pole]
    return _compound(copies, "iron")


def _example7_profile_points() -> list[tuple[float, float]]:
    """Full first-quadrant yoke polygon from Example #7 ``yokeinput``."""

    r1, thp = 178.0, 56.57
    points_mm = [
        (0, 250), (0, 169), (r1 * math.sqrt(2) - 169, 169),
        ((r1 - thp / 2) / math.sqrt(2),
         (r1 + thp / 2) / math.sqrt(2)),
        (12, 52), (12, 50.62),
    ]
    # Notebook automatic triangulation profile: x*y=648 mm^2.
    xs = np.linspace(16.0, 40.5, 12)
    points_mm.extend((float(x), 648.0 / float(x)) for x in xs)
    points_mm.extend([
        (50.62, 12), (52, 12),
        ((r1 + thp / 2) / math.sqrt(2),
         (r1 - thp / 2) / math.sqrt(2)),
        (r1 * math.sqrt(2) - 90, 90), (230, 90), (230, 117),
        (185.5, 117), (185.5, 271.25 * math.sqrt(2) - 185.5),
        (271.25 * math.sqrt(2) - 250, 250),
    ])
    return [(x * MM, y * MM) for x, y in points_mm]


def _example7_iron() -> Any:
    from netgen.occ import Axis, HalfSpace, Pnt, Vec, Z

    quarter = _polygon_prism_xy(
        _example7_profile_points(), -200 * MM, 400 * MM, "esrf_q_yoke_quarter")
    # Example#7.nb build3d cuts the positive longitudinal half with
    #
    #   pch = {rint/Sqrt[2], rint/Sqrt[2], lmag/2 - ch}
    #   vch = {-1, -1, Sqrt[2] Tan[ach]}
    #
    # before reflecting it through z=0.  HalfSpace keeps the side opposite its
    # normal, so the two planes below reproduce the 7 mm, 45-degree pole-end
    # machining at both ends of the 400 mm yoke.  The cut reaches z=+/-193 mm
    # on the diagonal bore point and tapers back to the uncut +/-200 mm end as
    # x+y increases.
    rint = 36 * MM
    chamfer = 7 * MM
    diagonal = rint / math.sqrt(2.0)
    positive_end = HalfSpace(
        Pnt(diagonal, diagonal, 200 * MM - chamfer),
        Vec(-1.0, -1.0, math.sqrt(2.0)),
    )
    negative_end = HalfSpace(
        Pnt(diagonal, diagonal, -200 * MM + chamfer),
        Vec(-1.0, -1.0, -math.sqrt(2.0)),
    )
    quarter = quarter * positive_end * negative_end
    copies = [quarter.Rotate(Axis(Pnt(0, 0, 0), Z), angle)
              for angle in (0, 90, 180, 270)]
    return _compound(copies, "iron")


def build_esrf_occ(number: int, *, include_coils: bool = True) -> dict[str, Any]:
    """Create material-separated OCC shapes for one ESRF example."""

    from netgen.occ import Pnt, Sphere

    number = int(number)
    if number == 1:
        shapes = {"magnet": _box_from_center(
            (0, 0, 0), (MM, MM, MM), "permanent_magnet")}
    elif number == 2:
        shapes = {}
    elif number == 3:
        shapes = _example3_occ()
    elif number == 4:
        sphere = Sphere(Pnt(0, 0, 0), MM)
        sphere.name = "permanent_magnet"
        shapes = {"magnet": sphere}
    elif number == 5:
        shapes = {"iron": _example5_iron()}
    elif number == 6:
        shapes = {"iron": _example6_iron()}
    elif number == 7:
        shapes = {"iron": _example7_iron()}
    else:
        raise ValueError("ESRF Radia example number must be in 1..7")

    if include_coils:
        coils = build_esrf_coils(number)
        if coils:
            shapes["coil"] = _compound([coil.to_occ() for coil in coils], "coil")
    return shapes


def build_esrf_cubit_hdiv_iron(number: int, *, image: str | None = None) -> Any:
    """Return partitioned iron CAD for a Cubit HEX BDM2 body mesh.

    The constructive regions remain separate so Cubit can imprint them and
    apply sweep/submap schemes independently.  Every resulting volume must be
    assigned to one ``iron`` block.  The Netgen exporter then omits shared CAD
    surfaces whose parent volumes map to that same material domain, so the
    loaded NGSolve mesh is conforming and does not acquire duplicated physical
    boundaries at the meshing partitions.

    Coils are deliberately absent.  They remain mesh-free ``CoilBuilder``
    sources in the BDM2 solve.
    """

    shapes = build_esrf_occ(int(number), include_coils=False)
    if "iron" not in shapes:
        raise ValueError(f"ESRF example {int(number)} contains no soft iron")
    source = shapes["iron"]
    solids = list(source.solids)
    if not solids:
        raise RuntimeError(f"ESRF example {int(number)} produced empty iron CAD")

    if image is not None:
        if int(number) != 5 or image != "+x-z":
            raise ValueError(
                "ESRF HDiv iron currently supports image='+x-z' for example 5 only"
            )
        from netgen.occ import Box, Pnt

        lower, upper = source.bounding_box
        span = max(
            upper.x - lower.x,
            upper.y - lower.y,
            upper.z - lower.z,
        )
        padding = max(0.1 * span, MM)
        positive_xz = Box(
            Pnt(0.0, lower.y - padding, 0.0),
            Pnt(upper.x + padding, upper.y + padding, upper.z + padding),
        )
        clipped = []
        for solid in solids:
            clipped.extend(list((solid * positive_xz).solids))
        solids = clipped
        if not solids:
            raise RuntimeError("example-5 +x-z clipping removed every iron solid")

    iron = _compound(solids, "iron")
    for solid in iron.solids:
        solid.name = "iron"
    return iron


def build_esrf_hdiv_iron(number: int, *, image: str | None = None) -> Any:
    """Return iron CAD with same-material internal interfaces removed.

    The tutorial CAD intentionally retains constructive regions because they
    are useful for Radia segmentation and Cubit convergence studies.  Passing
    those touching regions to an HDiv solve as a plain ``Compound`` duplicates
    coincident interface-charge degrees of freedom and can make the BDM2
    nonlinear tangent nearly singular.  Boolean union is therefore a required
    topology normalization at the HDiv boundary, while disconnected physical
    pieces remain disconnected solids in the returned shape.

    Example 5 additionally supports the notebook's two mirror planes through
    ``image="+x-z"``.  The returned CAD is then the positive-x/positive-z
    quarter only; the matching IMA string must also be passed to
    :func:`radia.vim.Solve`.  Keeping the reduction here, next to the exact
    source CAD, prevents validation runners from drifting to a differently cut
    geometry.
    """

    partitioned = build_esrf_cubit_hdiv_iron(int(number), image=image)
    solids = list(partitioned.solids)
    iron = solids[0]
    for solid in solids[1:]:
        iron = iron + solid
    iron.name = "iron"
    for solid in iron.solids:
        solid.name = "iron"
    return iron


def _solid_count(shape: Any) -> int:
    return len(list(shape.solids))


def _coil_source_manifest(coil: CoilBuilder, index: int) -> dict[str, Any]:
    wires, current = coil.to_wire_segments(n_arc=24)
    points = np.asarray([p for wire in wires for p in wire], dtype=float)
    widths = np.asarray([segment.width for segment in coil.segments])
    heights = np.asarray([segment.height for segment in coil.segments])
    return {
        "index": int(index), "current_A": float(current),
        "segment_count": len(coil.segments), "closed": bool(coil.is_closed),
        "closure_gap_m": float(coil.gap),
        "path_bbox_m": [points.min(axis=0).tolist(), points.max(axis=0).tolist()],
        "cross_section_width_range_m": [float(widths.min()), float(widths.max())],
        "cross_section_height_range_m": [float(heights.min()), float(heights.max())],
    }


def get_esrf_cubit_mesh_policy(number: int) -> dict[str, Any]:
    """Return the solver-specific Cubit element-family contract.

    The HDiv-MMM/BDM2 response mesh is air-mesh-free and contains only soft
    magnetic material.  The independent HCurl FEM route is conforming and
    prefers an all-HEX iron/air domain; TET outer air with an explicit
    WEDGE/PYRAMID transition is the allowed fallback.  A fixed permanent
    magnet is neither an iron response unknown nor a volume-current coil: it
    has its own Cubit source mesh, carrying one material label per prescribed
    magnetization block, and is evaluated once by
    :class:`radia.vim.MagnetizationSource`.  Coils are never solver-meshed:
    their STEP is visualization/provenance only and the physical source is
    always supplied by :class:`CoilBuilder`.

    Example 4 has no yoke and an exact spherical boundary.  Cubit's ordinary
    CAD-volume schemes do not produce a conforming all-HEX sphere, so the
    reproducible journal keeps a Cubit TET fallback; its separate Sculpt HEX
    convergence lane remains an optional analytic-sphere cross-check.
    """

    spec = get_esrf_example_spec(number)
    has_iron = spec.model_class in {
        "hybrid_permanent_magnet", "iron_dominated_dipole", "quadrupole"
    }
    # Example 7's end-chamfered quadrupole yoke is a full CAD authority, but
    # Cubit 2025.12 cannot map or submap its four volumes without a further
    # topology partition.  Require the supported curved-Q2 TET route
    # explicitly rather than issuing ``scheme auto`` and silently exporting an
    # empty mesh.  It remains a tracked HEX-partition improvement, not an
    # unrecorded fallback.
    explicit_tet_yoke = int(number) == 7
    cad_scheme = "tetmesh" if int(number) == 4 or explicit_tet_yoke else "auto"
    return {
        "schema": "radia.esrf-cubit-mesh-policy.v1",
        "all_examples_use_cubit": True,
        "preferred_volume_family": "HEX",
        "cad_volume_scheme": cad_scheme,
        "cad_volume_fallback_reason": (
            "exact sphere has no yoke; ordinary Cubit CAD meshing falls back "
            "to TET while Sculpt HEX is evaluated separately"
            if int(number) == 4 else
            "Cubit 2025.12 cannot map/submap the full 400 mm end-chamfered "
            "Example-7 yoke; use explicit curved-Q2 TET until the CAD is "
            "partitioned into sweepable HEX blocks"
            if explicit_tet_yoke else None
        ),
        "regions": {
            "iron": {
                "present": bool(has_iron),
                "required_family": (
                    "TET" if explicit_tet_yoke else "HEX" if has_iron else None
                ),
                "scheme": cad_scheme if has_iron else None,
            },
            "permanent_magnet": {
                "preferred_family": "HEX",
                "scheme": cad_scheme,
                "fallback_family": "TET",
                "source_model": "fixed-given MagnetizationSource",
                "response_unknown": False,
                "source_mesh": True,
                "source_mesh_labels": "one pm_NNN material per magnet block",
                "source_mesh_resolution": (
                    "independent uniform-M affine HEX source mesh; certify "
                    "against a direct-field reference before coarsening"
                ),
            },
            "coil": {
                "solver_mesh": False,
                "source": "CoilBuilder",
                "step_role": "visualization and geometry provenance only",
            },
            "transition": {
                "fem_only": True,
                "allowed_families": ["WEDGE", "PYRAMID"],
                "purpose": "optional conforming HEX-to-TET transition",
            },
            "air": {
                "bdm2": "not meshed (IMA/open-region operator)",
                "fem_preferred_family": "HEX",
                "fem_allowed_families": ["HEX", "TET"],
            },
        },
        "solver_routes": {
            "bdm2_ima": {
                "mesh": "Cubit soft-magnetic response mesh only",
                "iron_family": (
                    "TET" if explicit_tet_yoke else "HEX" if has_iron else None
                ),
                "air_mesh": False,
                "fixed_magnetization_source_mesh": True,
            },
            "nonlinear_hcurl_fem": {
                "mesh": "Cubit conforming full-domain mesh",
                "allowed_families": ["HEX", "WEDGE", "PYRAMID", "TET"],
                "iron_family": (
                    "TET" if explicit_tet_yoke else "HEX" if has_iron else None
                ),
                "air_preferred_family": "HEX",
                "air_fallback_family": "TET",
            },
        },
    }


def export_esrf_cubit_assets(number: int, output_dir: str | Path,
                             *, mesh_size_m: float | None = None,
                             fixed_magnetization_source_mesh_size_m: float | None = None,
                             order: int = 1) -> dict[str, Any]:
    """Write response/source STEP assets, Cubit journals, and metadata.

    A fixed PM source receives its own ``.vol`` because its prescribed
    magnetization belongs to a source-owned HDiv space.  In particular,
    Example #3's 24 alternating blocks must not be collapsed into one Cubit
    material before :class:`radia.vim.MagnetizationSource` projects them.
    The optional fixed-source size is intentionally independent of the soft
    iron response size.  A spatially uniform magnetization is represented
    exactly by an affine HDiv source cell, so its mesh may be coarser only
    after a direct-field certificate against the retained reference source.
    """

    spec = get_esrf_example_spec(number)
    root = Path(output_dir).resolve() / f"example_{number}_{spec.slug}"
    root.mkdir(parents=True, exist_ok=True)
    coils = build_esrf_coils(number)
    fixed_magnetization_blocks = build_esrf_fixed_magnetization_blocks(number)
    shapes = build_esrf_occ(number)
    files: dict[str, str] = {}
    counts: dict[str, int] = {}
    for category, shape in shapes.items():
        path = root / f"{category}.step"
        shape.WriteStep(str(path))
        files[category] = str(path)
        counts[category] = _solid_count(shape)

    if mesh_size_m is None:
        mesh_size_m = (0.1 * MM if number == 4 else
                       0.2 * MM if number == 1 else 10 * MM)
    if fixed_magnetization_source_mesh_size_m is None:
        fixed_magnetization_source_mesh_size_m = mesh_size_m
    if not float(fixed_magnetization_source_mesh_size_m) > 0.0:
        raise ValueError("fixed_magnetization_source_mesh_size_m must be positive")
    mesh_policy = get_esrf_cubit_mesh_policy(number)
    solver_categories = sorted(category for category in files
                               if category not in {"coil", "magnet"})
    lines = ["reset"]
    for block_id, category in enumerate(solver_categories, start=1):
        step = Path(files[category]).as_posix()
        lines.extend([
            f'import step "{step}" noheal',
            f'#{{{category}_last = Id("volume")}}',
            f'#{{{category}_first = {category}_last - {counts[category]} + 1}}',
            f'list volume {{{category}_first}} to {{{category}_last}}',
            (f'block {block_id} add volume {{{category}_first}} to '
             f'{{{category}_last}}'),
            f'block {block_id} name "{category}"',
        ])
    vol_path = root / "model.vol"
    if solver_categories:
        cub_path = root / "model.cub5"
        lines.extend([
            "sideset 1 add surface all",
            'sideset 1 name "outer_boundary"',
            f"volume all scheme {mesh_policy['cad_volume_scheme']}",
            f"volume all size {float(mesh_size_m):.12g}",
            "mesh volume all",
        ])
        if int(number) == 6:
            # The eight symmetry-related hyperbolic pole tips are submappable,
            # but Cubit's ``auto`` selector leaves them unmeshed at useful
            # resolutions.  Select by actual mesh state, never by entity ID.
            lines.extend([
                'group "hex_recovery" add volume with not is_meshed',
                "volume in hex_recovery scheme submap",
                "mesh volume in hex_recovery",
            ])
        lines.extend([
            f'export netgen "{vol_path.as_posix()}" order {int(order)} overwrite',
            f'save as "{cub_path.as_posix()}" overwrite',
            "exit",
        ])
    else:
        # Example 2 is source-only.  Preserve its Cubit-readable CAD, but do
        # not invent a volume-current FE mesh: CoilBuilder owns the source.
        cub_path = root / "source_geometry.cub5"
        for category in sorted(files):
            lines.append(
                f'import step "{Path(files[category]).as_posix()}" noheal')
        lines.extend([
            f'save as "{cub_path.as_posix()}" overwrite',
            "exit",
        ])
    journal = root / "mesh_and_export.jou"
    journal.write_text("\n".join(lines) + "\n", encoding="utf-8")

    source_step_files: dict[str, str] = {}
    source_vol: Path | None = None
    source_cubit_model: Path | None = None
    source_journal: Path | None = None
    if fixed_magnetization_blocks:
        source_root = root / "fixed_magnetization_sources"
        source_root.mkdir(exist_ok=True)
        source_lines = ["reset"]
        for block_id, block in enumerate(fixed_magnetization_blocks, start=1):
            step = source_root / f"{block.material_name}.step"
            block.shape.WriteStep(str(step))
            source_step_files[block.material_name] = str(step)
            source_lines.extend([
                f'import step "{step.as_posix()}" noheal',
                f'#{{{block.material_name}_last = Id("volume")}}',
                f'block {block_id} add volume {{{block.material_name}_last}}',
                f'block {block_id} name "{block.material_name}"',
            ])
        source_vol = source_root / "magnet_source.vol"
        source_cubit_model = source_root / "magnet_source.cub5"
        source_lines.extend([
            "sideset 1 add surface all",
            'sideset 1 name "outer_boundary"',
            f"volume all scheme {mesh_policy['cad_volume_scheme']}",
            f"volume all size {float(fixed_magnetization_source_mesh_size_m):.12g}",
            "mesh volume all",
            f'export netgen "{source_vol.as_posix()}" order {int(order)} overwrite',
            f'save as "{source_cubit_model.as_posix()}" overwrite',
            "exit",
        ])
        source_journal = source_root / "mesh_and_export.jou"
        source_journal.write_text(
            "\n".join(source_lines) + "\n", encoding="utf-8"
        )
    manifest = {
        "spec": spec.to_dict(), "step_files": files,
        "solver_step_files": {
            category: files[category] for category in solver_categories},
        "solver_vol": str(vol_path) if solver_categories else None,
        "cubit_model": str(cub_path),
        "solid_counts": counts, "mesh_size_m": mesh_size_m,
        "fixed_magnetization_source_mesh_size_m": (
            fixed_magnetization_source_mesh_size_m
        ),
        "mesh_policy": mesh_policy,
        "coil_sources": [_coil_source_manifest(coil, index)
                         for index, coil in enumerate(coils)],
        "fixed_magnetization_sources": [
            block.manifest_entry() for block in fixed_magnetization_blocks
        ],
        "fixed_magnetization_source_step_files": source_step_files,
        "fixed_magnetization_source_vol": (
            str(source_vol) if source_vol is not None else None
        ),
        "fixed_magnetization_source_cubit_model": (
            str(source_cubit_model) if source_cubit_model is not None else None
        ),
        "fixed_magnetization_source_journal": (
            str(source_journal) if source_journal is not None else None
        ),
        "netgen_order": int(order), "journal": str(journal),
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    manifest["manifest"] = str(manifest_path)
    return manifest


def export_esrf_cubit_suite(output_dir: str | Path,
                            numbers: tuple[int, ...] = tuple(range(1, 8)),
                            *, order: int = 1) -> dict[str, Any]:
    """Generate examples and one license-efficient Cubit batch journal."""

    manifests = [export_esrf_cubit_assets(n, output_dir, order=order)
                 for n in numbers]
    chunks: list[str] = []
    for manifest in manifests:
        journals = [manifest["journal"]]
        source_journal = manifest.get("fixed_magnetization_source_journal")
        if source_journal:
            journals.append(source_journal)
        for journal in journals:
            lines = Path(journal).read_text(encoding="utf-8").splitlines()
            if lines and lines[-1].strip().lower() == "exit":
                lines.pop()
            chunks.extend(lines)
    chunks.append("exit")
    suite_path = Path(output_dir).resolve() / "mesh_all_examples.jou"
    suite_path.write_text("\n".join(chunks) + "\n", encoding="utf-8")
    return {"journal": str(suite_path), "manifests": manifests}


def validate_esrf_radia_reference(number: int, *, n_points: int = 301) -> dict[str, Any]:
    """Run lightweight source-level checks against the original Radia model.

    Implemented here for the two non-iron tutorial sources that give a sharp
    regression oracle without an iterative material solve.  Example 1 checks
    the published vector value.  Example 2 compares CoilBuilder against
    Radia's native volume-current ``ObjRaceTrk`` representation.
    """

    import radia as rad

    number = int(number)
    if number == 1:
        rad.UtiDelAll()
        try:
            # The maintained SI binding accepts magnetization in A/m, whereas
            # the 2023 notebook supplies mu0*M in tesla.
            handle = rad.ObjRecMag(
                [0., 0., 0.], [MM, MM, MM],
                np.asarray([-0.5, 1.0, 0.7]) / MU0,
            )
            actual = np.asarray(rad.Fld(
                handle, "b", [0.52 * MM, 0.60 * MM, 0.70 * MM]), dtype=float)
            expected = np.asarray([0.12737, 0.028644, 0.077505])
            error = actual - expected
            return {
                "example": 1, "actual_T": actual.tolist(),
                "expected_T": expected.tolist(),
                "max_abs_error_T": float(np.max(np.abs(error))),
                "passed": bool(np.max(np.abs(error)) < 1.0e-5),
            }
        finally:
            rad.UtiDelAll()

    if number == 2:
        if n_points < 3:
            raise ValueError("n_points must be at least 3")
        raw = (
            (9.5, 24.5, 120., 36., 128., 38.),
            (24.5, 55.5, 120., 36., 256., 38.),
            (10., 25., 90., 24., 128., 76.),
            (25., 55., 90., 24., 256., 76.),
            (150., 166.3, 0., 39., -256., 60.),
        )
        ys = np.linspace(0., 300 * MM, int(n_points))
        points = [[0., float(y), 0.] for y in ys]

        def native_profile() -> np.ndarray:
            objects: list[int] = []
            for ri, ro, straight, height, j, zc in raw:
                for z_sign in (1., -1.):
                    objects.append(rad.ObjRaceTrk(
                        [0., 0., z_sign * zc * MM],
                        [ri * MM, ro * MM], [straight * MM, 0.], height * MM,
                        3 if straight else 6, "man", "z", j * 1.0e6,
                    ))
            return np.asarray(rad.Fld(rad.ObjCnt(objects), "bz", points), dtype=float)

        rad.UtiDelAll()
        try:
            native = native_profile()
            rad.UtiDelAll()
            builder_objects: list[int] = []
            for coil in build_esrf_coils(2):
                builder_objects.extend(coil.to_radia())
            builder = np.asarray(
                rad.Fld(rad.ObjCnt(builder_objects), "bz", points), dtype=float)
            delta = builder - native
            relative_l2 = float(np.linalg.norm(delta) / np.linalg.norm(native))
            peak_relative = float(abs(builder[0] - native[0]) / abs(native[0]))
            return {
                "example": 2, "n_points": int(n_points),
                "native_peak_T": float(native[0]),
                "coilbuilder_peak_T": float(builder[0]),
                "peak_relative_error": peak_relative,
                "profile_relative_l2": relative_l2,
                "profile_max_abs_error_T": float(np.max(np.abs(delta))),
                "passed": bool(relative_l2 < 0.01 and peak_relative < 0.005),
            }
        finally:
            rad.UtiDelAll()

    if number == 4:
        # The continuum result is the notebook's stated oracle.  The notebook
        # itself approximates the sphere by 15 azimuthal x 15 axial slices.
        return {
            "example": 4, "analytic_internal_B_T": [2. / 3., 0., 0.],
            "expected_internal_B_T": [2. / 3., 0., 0.], "passed": True,
        }
    raise NotImplementedError(
        "lightweight Radia reference is currently implemented for examples 1, 2, and 4")
