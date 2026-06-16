"""Electromagnetic material presets and unified property management.

Provides a single source of truth for material properties (sigma, mu_r,
BH curves) used across all Radia calc scripts, PEEC topology, and GUI
panels.  Replaces the scattered ``sigma_for_material()``,
``get_bh_curve()``, and per-script argparse defaults that previously
made it easy for ``--material copper`` to silently keep the steel
default ``sigma=2e6``.

Usage (calc script)::

    from em_material import EMMaterial, add_material_args

    add_material_args(parser)                # consistent --material / --sigma / --mu-r
    args = parser.parse_args()
    mat = EMMaterial.from_args(args)          # auto-resolves sigma, mu_r, BH
    Z_s = mat.dowell_Zs(frequency, R)        # surface impedance
    esim = mat.create_esim_solver(freq, R)   # ESIM cell-problem solver

Usage (PEEC / Radia C++ bridge)::

    mat = EMMaterial.from_name("copper")
    builder.add_segment(n1, n2, w, h, sigma=mat.sigma)
    radia_mat = mat.create_radia_material()  # rad.MatLin or rad.MatSatIsoTab
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

MU_0 = 4e-7 * np.pi

# ============================================================
# BH curve data (100 points, from CEFC 2020 nonlinear model)
# Source data file: examples/cubit_panels/accel_magnet/BH.txt
# Units: H [A/m], B [T]
# ============================================================
STEEL_BH: List[List[float]] = [
    [0.0, 0.0],
    [13.898, 0.22296], [15.397, 0.25304], [17.058, 0.28380],
    [18.898, 0.31552], [20.936, 0.34852], [23.194, 0.38323],
    [25.696, 0.42011], [28.467, 0.45974], [31.538, 0.50272],
    [34.939, 0.54965], [38.708, 0.60110], [42.883, 0.65744],
    [47.508, 0.71868], [52.632, 0.78437], [58.309, 0.85340],
    [64.598, 0.92403], [71.565, 0.99400], [79.284, 1.06090],
    [87.836, 1.12254], [97.310, 1.17738], [107.806, 1.22465],
    [119.433, 1.26440], [132.315, 1.29727], [146.587, 1.32427],
    [162.397, 1.34654], [179.913, 1.36518], [199.319, 1.38116],
    [220.817, 1.39530], [244.634, 1.40821], [271.020, 1.42039],
    [300.252, 1.43217], [332.636, 1.44381], [368.514, 1.45547],
    [408.262, 1.46728], [452.296, 1.47930], [501.081, 1.49157],
    [555.127, 1.50410], [615.002, 1.51691], [681.335, 1.52999],
    [754.823, 1.54332], [836.238, 1.55689], [926.433, 1.57068],
    [1026.357, 1.58467], [1137.059, 1.59883], [1259.701, 1.61315],
    [1395.571, 1.62761], [1546.096, 1.64220], [1712.856, 1.65688],
    [1897.603, 1.67166], [2102.276, 1.68651], [2329.025, 1.70142],
    [2580.231, 1.71638], [2858.532, 1.73137], [3166.850, 1.74640],
    [3508.423, 1.76144], [3886.837, 1.77649], [4306.067, 1.79154],
    [4770.514, 1.80658], [5285.057, 1.82162], [5855.097, 1.83664],
    [6486.621, 1.85165], [7186.261, 1.86663], [7961.362, 1.88158],
    [8820.066, 1.89651], [9771.388, 1.91142], [10825.319, 1.92629],
    [11992.926, 1.94114], [13286.469, 1.95596], [14719.532, 1.97077],
    [16307.164, 1.98555], [18066.037, 2.00033], [20014.619, 2.01510],
    [22173.373, 2.02987], [24564.968, 2.04466], [27214.517, 2.05948],
    [30149.844, 2.07434], [33401.772, 2.08927], [37004.449, 2.10428],
    [40995.707, 2.11941], [45417.458, 2.13468], [50316.133, 2.15014],
    [55743.174, 2.16583], [61755.569, 2.18181], [68416.455, 2.19815],
    [75795.776, 2.21493], [83971.022, 2.23225], [93028.041, 2.25023],
    [103061.940, 2.26901], [114178.084, 2.28877], [126493.203, 2.30972],
    [140136.616, 2.33211], [155251.592, 2.35623], [171996.852, 2.38247],
    [190548.237, 2.41122], [211100.554, 2.44300], [233869.620, 2.47834],
    [259094.531, 2.51786], [287040.173, 2.56214], [318000.0, 2.61173],
]


def load_bh_file(path: str) -> List[List[float]]:
    """Load a 2-column BH file (H [A/m], B [T]).

    Public helper used by `calc_inductance` / `calc_fem_coilmesh` to
    parse `--bh-file` for ESIM Karl iteration.
    """
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Invalid BH file format (need 2+ columns): {path}")
    return data[:, :2].tolist()


# Backward-compatibility alias: existing callers may still import
# `_load_bh_file`.  Prefer the public `load_bh_file` in new code.
_load_bh_file = load_bh_file


@dataclass
class EMMaterial:
    """Electromagnetic material properties.

    Attributes:
        name:      Human-readable name (e.g. "steel", "copper").
        sigma:     Electrical conductivity [S/m].
        mu_r:      Relative permeability (used when *bh_curve* is None).
        bh_curve:  Nonlinear BH data [[H, B], ...] or None for linear.
        hys_file:  Path to .hys hysteresis file (accelerator use), or "".
    """

    name: str
    sigma: float
    mu_r: float
    bh_curve: Optional[List[List[float]]] = None
    hys_file: str = ""

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def rho(self) -> float:
        """Resistivity [Ohm*m]."""
        return 1.0 / self.sigma

    @property
    def is_magnetic(self) -> bool:
        """True if mu_r > 1 or a BH curve is present."""
        return self.bh_curve is not None or self.mu_r > 1.0

    @property
    def is_nonlinear(self) -> bool:
        """True if a BH curve governs permeability."""
        return self.bh_curve is not None

    @property
    def is_conducting(self) -> bool:
        """True if sigma > 0."""
        return self.sigma > 0

    def skin_depth(self, frequency: float) -> float:
        """Classical skin depth [m] at *frequency* [Hz].

        Uses *mu_r* (ignoring BH nonlinearity).  Returns inf for DC or
        non-conducting materials.
        """
        if frequency <= 0 or self.sigma <= 0:
            return float("inf")
        omega = 2 * math.pi * frequency
        return math.sqrt(2.0 * self.rho / (omega * MU_0 * self.mu_r))

    # ------------------------------------------------------------------
    # ESIM / SIBC helpers
    # ------------------------------------------------------------------

    def dowell_Zs(self, frequency: float, R: float) -> complex:
        """Linear surface impedance via Dowell tanh formula.

        Args:
            frequency: Operating frequency [Hz].
            R: Characteristic radius / half-thickness [m].

        Returns:
            Complex surface impedance Z_s.
        """
        omega = 2 * math.pi * frequency
        if omega <= 0 or R <= 0 or self.sigma <= 0:
            return complex(0, 0)
        mu_eff = MU_0 * self.mu_r
        delta = math.sqrt(2 * self.rho / (omega * mu_eff))
        xi = R / delta
        gamma_a = complex(1, 1) * xi
        if xi > 30.0:
            # Thick-conductor / high-frequency limit: tanh(gamma_a) -> 1 to far
            # below machine precision (|tanh - 1| < 1e-26 at xi = 30).  Use it
            # explicitly: np.tanh OVERFLOWS to nan past xi ~ 710 emitting only a
            # RuntimeWarning (NOT an OverflowError), so the old `except
            # OverflowError` never fired and dowell_Zs SILENTLY returned nan for
            # thick conductors (e.g. an induction-heating steel workpiece at
            # high frequency).
            return (self.rho / R) * gamma_a          # = (1+1j) * rho / delta
        return (self.rho / R) * gamma_a * np.tanh(gamma_a)

    def create_esim_solver(self, frequency: float, half_thickness: float,
                           geometry: str = "cylinder"):
        """Create an ESIMFiniteSlabSolver for this material.

        Lazy-imports esim_cell_problem to avoid import overhead when not
        needed.
        """
        import sys as _sys
        radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)))
        if radia_src not in _sys.path:
            _sys.path.insert(0, radia_src)
        from esim_cell_problem import ESIMFiniteSlabSolver

        return ESIMFiniteSlabSolver(
            half_thickness=half_thickness,
            bh_curve=self.bh_curve,
            sigma=self.sigma,
            frequency=frequency,
            mu_r=self.mu_r if self.bh_curve is None else None,
            n_nodes=200,
            geometry=geometry,
        )

    # ------------------------------------------------------------------
    # Radia C++ material bridge
    # ------------------------------------------------------------------

    def create_radia_material(self):
        """Create a Radia material handle (MatLin or MatSatIsoTab).

        Requires ``import radia as rad`` to be available.
        """
        import radia as rad

        if self.bh_curve is not None:
            return rad.MatSatIsoTab(self.bh_curve)
        return rad.MatLin(self.mu_r)

    # ------------------------------------------------------------------
    # BH curve access (for callers that need raw data)
    # ------------------------------------------------------------------

    def get_bh_curve(self) -> Tuple[Optional[List[List[float]]], Optional[float]]:
        """Return (bh_curve, mu_r) in the legacy format.

        Returns:
            (bh_curve, None) for nonlinear materials, or
            (None, mu_r) for linear materials.
        """
        if self.bh_curve is not None:
            return self.bh_curve, None
        return None, self.mu_r

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_name(cls, name: str, *, bh_file: str = "",
                  hys_file: str = "") -> "EMMaterial":
        """Look up a material by preset name.

        Supported names: "steel", "elf_steel", "copper", "aluminum".
        Raises ValueError for unknown names.
        """
        key = name.lower().strip()
        if key in ("steel", "elf_steel"):
            bh = STEEL_BH
            if bh_file and bh_file != "(built-in Steel)" and os.path.exists(bh_file):
                bh = _load_bh_file(bh_file)
            return cls(name=key, sigma=2e6, mu_r=100.0,
                       bh_curve=bh, hys_file=hys_file)
        if key == "copper":
            return cls(name="copper", sigma=5.8e7, mu_r=1.0)
        if key == "aluminum":
            return cls(name="aluminum", sigma=3.5e7, mu_r=1.0)
        if key == "linear":
            # Linear mode: mu_r set by caller via --mu-r (default 1000
            # for accelerator yokes).  No BH curve.
            return cls(name="linear", sigma=2e6, mu_r=1000.0)
        if key == "hysteresis":
            # Hysteresis mode: uses .hys file, BH from preset steel.
            return cls(name="hysteresis", sigma=2e6, mu_r=100.0,
                       bh_curve=STEEL_BH, hys_file=hys_file)
        raise ValueError(
            f"Unknown material preset: {name!r}.  "
            f"Available: steel, copper, aluminum, linear, "
            f"hysteresis (or use 'custom')."
        )

    @classmethod
    def custom(cls, sigma: float, mu_r: float, *,
               bh_file: str = "", hys_file: str = "") -> "EMMaterial":
        """Create a custom material with explicit properties."""
        bh = None
        if bh_file and os.path.exists(bh_file):
            bh = _load_bh_file(bh_file)
        return cls(name="custom", sigma=sigma, mu_r=mu_r,
                   bh_curve=bh, hys_file=hys_file)

    @classmethod
    def from_args(cls, args) -> "EMMaterial":
        """Build an EMMaterial from an argparse.Namespace.

        Expects the arguments added by :func:`add_material_args`.  When
        ``--material`` is a preset name the corresponding sigma/mu_r are
        used **regardless** of ``--sigma`` / ``--mu-r`` defaults; explicit
        user overrides (non-default values) still win.

        The logic:
        - ``--material steel`` -> sigma=2e6, mu_r=100, BH=STEEL_BH
        - ``--material copper`` -> sigma=5.8e7, mu_r=1.0
        - ``--material custom --sigma 1e7 --mu-r 50`` -> as given
        - ``--material steel --sigma 5e6`` -> sigma=5e6 (override), rest from preset
        """
        material = getattr(args, "material", "steel")
        sigma_arg = getattr(args, "sigma", None)
        mu_r_arg = getattr(args, "mu_r", None)
        bh_file = getattr(args, "bh_file", "") or ""
        hys_file = getattr(args, "hys_file", "") or ""

        if material == "custom":
            sigma = sigma_arg if sigma_arg is not None else 2e6
            mu_r = mu_r_arg if mu_r_arg is not None else 1.0
            return cls.custom(sigma, mu_r, bh_file=bh_file,
                              hys_file=hys_file)

        # Preset: start from canonical values
        mat = cls.from_name(material, bh_file=bh_file, hys_file=hys_file)

        # Allow explicit overrides.  We detect "explicit" by checking
        # whether the parser default was touched.  For robustness we
        # also accept 0 as "not set" for mu_r (legacy convention in
        # calc_inductance.py where --mu-r default=0 means "auto").
        if sigma_arg is not None and sigma_arg != _SENTINEL_SIGMA:
            mat = EMMaterial(name=mat.name, sigma=sigma_arg, mu_r=mat.mu_r,
                             bh_curve=mat.bh_curve, hys_file=mat.hys_file)
        if mu_r_arg is not None and mu_r_arg > 0 and mu_r_arg != _SENTINEL_MU_R:
            mat = EMMaterial(name=mat.name, sigma=mat.sigma, mu_r=mu_r_arg,
                             bh_curve=mat.bh_curve, hys_file=mat.hys_file)

        return mat


# Sentinel defaults so from_args can distinguish "user passed --sigma"
# from "argparse default".  These are intentionally unusual values that
# no real material would have.
_SENTINEL_SIGMA = -1.0
_SENTINEL_MU_R = -1.0


# ============================================================
# argparse integration
# ============================================================

def add_material_args(parser, *,
                      default_material: str = "steel",
                      include_custom: bool = True,
                      include_hys: bool = False,
                      sigma_help: str = "Conductivity [S/m] "
                                        "(overrides preset if given)",
                      mu_r_help: str = "Relative permeability "
                                       "(overrides preset if given)"):
    """Add standard ``--material``, ``--sigma``, ``--mu-r``, ``--bh-file``
    arguments to *parser*.

    The defaults are sentinel values so that :meth:`EMMaterial.from_args`
    can distinguish "user explicitly passed --sigma 2e6" from "argparse
    default".  This way, ``--material copper`` correctly sets sigma to
    5.8e7 without the user having to also pass ``--sigma 5.8e7``.
    """
    choices = ["steel", "copper", "aluminum"]
    if include_custom:
        choices.append("custom")
    if include_hys:
        choices.extend(["elf_steel", "linear", "hysteresis"])

    parser.add_argument(
        "--material", default=default_material, choices=choices,
        help=f"Material preset (default: {default_material}).  "
             f"Each preset sets sigma and mu_r automatically.")
    parser.add_argument(
        "--sigma", type=float, default=_SENTINEL_SIGMA,
        help=sigma_help)
    parser.add_argument(
        "--mu-r", type=float, default=_SENTINEL_MU_R,
        help=mu_r_help)
    parser.add_argument(
        "--bh-file", default="",
        help="BH curve file (2-column: H[A/m] B[T], overrides built-in)")
    if include_hys:
        parser.add_argument(
            "--hys-file", default="",
            help=".hys hysteresis file")
