"""CEFC 2020 quadrupole ("Q-mag") validation case: geometry, coils, iron law, mesh journal.

The magnet is the laboratory's CEFC 2020 test quadrupole: a 4-pole iron yoke of 60 mm length
(beam axis z), pole tips on the +-x / +-y axes with a 40 mm aperture (pole face at 20 mm), and two
racetrack coils per pole.  The iron CAD ``qmag_iron.step`` is the laboratory's ACIS model of the
yoke (24 solids, every one a z-extrusion) exported by Coreform Cubit 2025.12 after scaling to
metres; the exported STEP re-imports with the same total volume (1.0247373e-3 m^3).

Coils (the ``radObjRaceTrk`` definition of the original design notebook, mm): per pole an inner
racetrack with corner radii 2..15, straight sections 60 (beam) x 26 (transverse), height 40
(pole coordinate 30..70), and an outer racetrack with radii 15..31, the same straights, height 20
(50..70).  The +-y coils carry +J, the +-x coils -J, which makes B_x = G x, B_y = -G y in the
aperture with G > 0 for J > 0.  The design point is J = 3 A/mm^2; the nonlinear sweep covers
0.5..10 A/mm^2 with the iron law of ``iron_bh_table.json``.

Observable: B_perp = -(B_x - B_y)/sqrt(2) along the diagonal x = y, z = 0 (31 points, r = -15..15
mm, spacing 1 mm); the summary value is B_perp at r = 15 mm (the gradient is B_perp / r).
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STEP_PATH = HERE / "qmag_iron.step"
BH_PATH = HERE / "iron_bh_table.json"
MM = 1.0e-3
MU0 = 4.0e-7 * math.pi

# Coil geometry (metres); the pole coordinate is the coil axis, the beam axis is z.
INNER_COIL = {"radii": (2.0 * MM, 15.0 * MM), "straight_beam": 60.0 * MM, "straight_transverse": 26.0 * MM,
              "height": 40.0 * MM, "centre_pole": 50.0 * MM}
OUTER_COIL = {"radii": (15.0 * MM, 31.0 * MM), "straight_beam": 60.0 * MM, "straight_transverse": 26.0 * MM,
              "height": 20.0 * MM, "centre_pole": 60.0 * MM}
POLE_ANGLES_DEG = (0.0, 90.0, 180.0, 270.0)      # +x, +y, -x, -y
# Azimuthal current sign about each coil's own axis (+x, +y, -x, -y).  A positive ObjRaceTrk
# density gives a magnetic moment along the +axis, so the pattern below points every +-x moment
# outward and every +-y moment inward: the quadrupole B_x = G x, B_y = -G y with G > 0.
POLE_CURRENT_SIGN = (+1.0, -1.0, -1.0, +1.0)
ARC_SEGMENTS = 16
IRON_VOLUME_M3 = 1.0247373119738649e-3
DESIGN_CURRENT_DENSITY_A_PER_MM2 = 3.0
OBSERVATION_RADIUS_M = 15.0 * MM


def build_qmag_coils(j_a_per_mm2: float, *, arc_segments: int = ARC_SEGMENTS) -> tuple[int, dict]:
    """Return (Radia container, manifest) of the eight racetracks for current density ``j`` [A/mm^2].

    ``ObjRaceTrk`` orders its two straight lengths cyclically after the axis (axis x: [l_y, l_z],
    axis y: [l_z, l_x], axis z: [l_x, l_y]; probed 2026-09-06 through the bore/return field
    signs), so the +-y coils use axis ``y`` with ``[l_beam, l_transverse]`` and the +-x coils axis
    ``x`` with ``[l_transverse, l_beam]``.  The height is centred on the pole coordinate.
    """
    import radia as rad

    j = float(j_a_per_mm2) * 1.0e6          # A/m^2
    objects, manifest = [], []
    for angle_deg, sign in zip(POLE_ANGLES_DEG, POLE_CURRENT_SIGN):
        direction = (math.cos(math.radians(angle_deg)), math.sin(math.radians(angle_deg)))
        axis = "x" if abs(direction[0]) > 0.5 else "y"
        for name, coil in (("inner", INNER_COIL), ("outer", OUTER_COIL)):
            centre = [direction[0] * coil["centre_pole"], direction[1] * coil["centre_pole"], 0.0]
            lengths = ([coil["straight_transverse"], coil["straight_beam"]] if axis == "x"
                       else [coil["straight_beam"], coil["straight_transverse"]])
            objects.append(rad.ObjRaceTrk(centre, list(coil["radii"]), lengths, coil["height"],
                                          int(arc_segments), "man", axis, sign * j))
            manifest.append({"pole_angle_deg": angle_deg, "axis": axis, "centre_m": centre, "coil": name,
                             "current_density_A_per_m2": sign * j,
                             **{k: (list(v) if isinstance(v, tuple) else v) for k, v in coil.items()}})
    return rad.ObjCnt(objects), {"coils": manifest, "arc_segments": int(arc_segments),
                                 "current_density_A_per_mm2": float(j_a_per_mm2)}


def observation_points() -> np.ndarray:
    """31 points (metres) on the diagonal: (n, n, 0)/sqrt(2) mm for n = -15..15."""
    n = np.arange(-15, 16, dtype=float)
    s = n / math.sqrt(2.0) * MM
    return np.column_stack([s, s, np.zeros_like(s)])


def b_perp(field: np.ndarray) -> np.ndarray:
    """-(B_x - B_y)/sqrt(2): the component perpendicular to the x = y diagonal."""
    field = np.asarray(field, dtype=float).reshape(-1, 3)
    return -(field[:, 0] - field[:, 1]) / math.sqrt(2.0)


def load_bh_table(path: Path = BH_PATH) -> list[list[float]]:
    """The iron law as the ``[[H A/m, B T], ...]`` table consumed by ``vim.Solve(bh_table=...)``."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("columns") != ["H_A_per_m", "B_T"]:
        raise ValueError(f"{path}: expected columns [H_A_per_m, B_T], got {payload.get('columns')}")
    data = np.asarray(payload["rows"], dtype=float)
    if data.ndim != 2 or data.shape[1] != 2 or np.any(np.diff(data[:, 0]) <= 0.0) or np.any(np.diff(data[:, 1]) <= 0.0):
        raise ValueError(f"{path}: rows must be strictly increasing in both H and B")
    return [[float(h), float(b)] for h, b in data]


def write_mesh_journal(journal: Path, vol: Path, *, size_m: float, order: int, step: Path = STEP_PATH) -> None:
    """Conforming z-swept HEX journal: imprint/merge every solid, sweep each from its z_min to z_max faces."""
    if not float(size_m) > 0.0 or int(order) < 1:
        raise ValueError("size_m must be positive and order >= 1")
    z_lo, z_hi = -30.0 * MM, 30.0 * MM
    tol = 1.0e-3 * (z_hi - z_lo)
    lines = [
        "reset",
        f'import step "{Path(step).resolve().as_posix()}" noheal',
        '#{iron_last = Id("volume")}',
        "#{iron_first = 1}",
        "imprint volume all",
        "merge volume all",
        "block 1 add volume all",
        'block 1 name "iron"',
        "sideset 1 add surface all",
        'sideset 1 name "outer_boundary"',
        "#{_v = iron_first}",
        "#{_n = iron_last - iron_first + 1}",
        "#{Loop(_n)}",
        (f"volume {{_v}} scheme sweep source surface in volume {{_v}} with z_coord < {z_lo + tol:.12g} "
         f"target surface in volume {{_v}} with z_coord > {z_hi - tol:.12g}"),
        "#{_v++}",
        "#{EndLoop}",
        f"volume all size {float(size_m):.12g}",
        "mesh volume all",
        "list volume with not is_meshed",
        f'export netgen "{Path(vol).resolve().as_posix()}" order {int(order)} overwrite',
        "exit",
    ]
    Path(journal).write_text("\n".join(lines) + "\n", encoding="ascii")


def run_cubit(cubit: Path, journal: Path, log: Path) -> None:
    command = [str(cubit), "-nographics", "-batch", "-nojournal", "-input", str(journal)]
    with open(log, "w", encoding="utf-8", errors="replace") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    text = Path(log).read_text(encoding="utf-8", errors="replace")
    if "did not mesh" in text or "Exported Netgen Vol" not in text:
        raise RuntimeError(f"Cubit did not produce the full mesh; see {log}")
