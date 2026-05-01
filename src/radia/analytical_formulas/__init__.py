"""Analytical formulas for electromagnetic field computation.

This package collects closed-form analytical expressions taken from the
review series

    若尾真治, 五十嵐一, 藤原耕二, 野口聡, 松尾哲司, 亀有昭久,
    "Useful Formulas of Analytical Integration in Electromagnetic Field
    Computations (Part 1..5)", IEE Japan Joint Technical Meeting on
    Static Apparatus and Rotating Machinery (SA / RM), 2002-2004.

Each module pins the corresponding equation numbers from the original
PDF series in its docstring so that future contributors can locate the
derivation. The formulas are reference implementations: priority is
correctness and traceability, not raw speed.

Modules
-------
ellipsoid       Rotational ellipsoid demagnetization factor + torque
                (Part 5, eq 39-44)

ac_locus        Major/minor axis of the time-locus ellipse traced by an
                AC vector quantity (B, J, ...) in linear steady-state
                eddy current analysis (Part 5, eq 29-37)

shielding       Static shielding factor of a magnetic cylindrical or
                spherical shell in a uniform external field
                (Part 1, eq 23-24)

rect_magnet_2d  2D uniformly magnetized rectangular bar: vector
                potential A_z and field B_x, B_y (Part 2, eq 2-3)

plate_eddy      Eddy current in a thin rectangular plate under a
                slowly-varying perpendicular B-field, in terms of the
                current vector potential T_z (Part 1, eq 26-27)
"""

from .ellipsoid import (
    demag_factor_prolate,
    demag_factor_oblate,
    demag_factor_rotational,
    ellipsoid_internal_field,
    ellipsoid_torque,
)
from .ac_locus import (
    ac_locus_axes,
    ac_locus_axes_batch,
)
from .shielding import (
    shielding_factor_cylinder,
    shielding_factor_sphere,
)
from .rect_magnet_2d import (
    rect_magnet_2d_A,
    rect_magnet_2d_B,
)
from .plate_eddy import (
    plate_eddy_T,
    plate_eddy_J,
)

__all__ = [
    # ellipsoid
    "demag_factor_prolate",
    "demag_factor_oblate",
    "demag_factor_rotational",
    "ellipsoid_internal_field",
    "ellipsoid_torque",
    # ac_locus
    "ac_locus_axes",
    "ac_locus_axes_batch",
    # shielding
    "shielding_factor_cylinder",
    "shielding_factor_sphere",
    # rect_magnet_2d
    "rect_magnet_2d_A",
    "rect_magnet_2d_B",
    # plate_eddy
    "plate_eddy_T",
    "plate_eddy_J",
]
