"""radia_vim — Volume Integral Method on a single hex cell, NGSolve extension.

Implements the Newton-kernel volume Galerkin method for HDiv div-free interior
basis on an axis-aligned hex (cuboid). Used to extract Cauer ladder networks
for eddy currents in isolated 3D conductors with the BC J·n = 0 (correct
"isolated conductor in vacuum" formulation, no air-box pathology).

Templated on numeric precision:
    - double: 6-25 clean Cauer pairs (orders 3-4)
    - cpp_bin_float_quad: 25-40 clean pairs (orders 5-6) -- requires Boost
    - mpfr_float (60+ digits): 50+ pairs (orders 7-8) -- requires MPFR

See `precision_sensitivity_study.py` (in radia-axifemm/scripts) for the
empirical precision/depth tradeoff this design is based on.

Quick usage (Phase F-1 bootstrap):

    from radia_vim import version, precision_info, newton_kernel_double

    print(version())
    print(precision_info())
    print(newton_kernel_double([0, 0, 0], [1, 0, 0]))    # -> 1.0

Reference:
    Sauter & Schwab, "Boundary Element Methods", Springer 2011, Chap. 5.
    Schöberl & Zaglmayr, "High order Nédélec elements with local complete
    sequence properties," IJ-COMP 24(2):374-384, 2005.
    Hiruma & Igarashi, "Eddy-Current Analysis Using Cauer Ladder Network
    Method," IEEE Trans. Magn. 56(3), 2020.
    Nagamine et al., "Verified Numerical Computations of the Cauer Network
    Representation of a Square Prism Conductor," 2026 (preprint).
"""

# IMPORTANT: import ngsolve before the C++ module so that NGSolve's shared
# libraries are loaded into the process before ours.
import ngsolve  # noqa: F401

from .radia_vim import (  # noqa: F401
    version, hello,
    PrecisionInfo, precision_info,
    newton_kernel_double,
)

# Conditionally-compiled high-precision entry points.
try:
    from .radia_vim import newton_kernel_quad  # noqa: F401
except ImportError:
    pass

try:
    from .radia_vim import newton_kernel_mpfr  # noqa: F401
except ImportError:
    pass

# Future phase exports (Phase F-2, F-3, F-4) — wrap in try/except so the
# bootstrap build still imports cleanly.
try:
    from .radia_vim import (  # noqa: F401
        HDivDivFreeHexBasis,
        HexBasisBlock,
        HexBasisDofKey,
    )
except ImportError:
    pass

try:
    from .radia_vim import (  # noqa: F401
        compute_K_entry_spherical_duffy,
        gauss_legendre_m1_1,
    )
except ImportError:
    pass

try:
    from .radia_vim import (  # noqa: F401
        assemble_K_bare,
    )
except ImportError:
    pass

__version__ = "0.1.0"
