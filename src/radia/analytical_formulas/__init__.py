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

Closed-form fields and global quantities:

ellipsoid          Rotational ellipsoid demag factor + torque
                   (Part 5, eq 39-44)
ac_locus           Major / minor axis of an AC vector locus
                   (Part 5, eq 29-37)
shielding          Static shielding factor of a magnetic shell
                   (Part 1, eq 23-24)
rect_magnet_2d     2D rectangular bar A_z, B_x, B_y (Part 2, eq 2-3)
plate_eddy         Thin rectangular-plate eddy current
                   (Part 1, eq 26-27)

Coil and line geometries:

solenoid_central   Rectangular cross-section solenoid central /
                   axial field via Fabri's form factor
                   (Part 4, eq 26-27)
three_phase_line   Three-phase straight-line and helical-line field;
                   triangle, planar, hexagon arrangements + far-field
                   asymptotes (Part 4 §5, Part 5 §3)

Numerics and special functions:

elliptic_integrals K(k), E(k) Hastings polynomial approximations
                   (Part 3 §3, Tables 1-2)
gauss_legendre     Gauss-Legendre nodes and weights and convenience
                   integrators on a general interval (Part 3 §4)
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
from .solenoid_central import (
    fabri_F,
    central_field as solenoid_central_field,
    axial_field as solenoid_axial_field,
)
from .three_phase_line import (
    vector_potential_z,
    field_xy,
    triangle_far_field_amplitude,
    planar_far_field_amplitude,
    helical_near_field_amplitude,
    helical_far_field_amplitude,
    triangle_positions,
    planar_positions,
    hexagon_positions,
    balanced_three_phase_currents,
    balanced_six_phase_currents,
)
from .elliptic_integrals import (
    K_hastings_2,
    E_hastings_2,
    K_hastings_4,
    E_hastings_4,
)
from .gauss_legendre import (
    nodes_weights as gauss_legendre_nodes_weights,
    integrate as gauss_legendre_integrate,
    integrate_2d as gauss_legendre_integrate_2d,
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
    # solenoid_central
    "fabri_F",
    "solenoid_central_field",
    "solenoid_axial_field",
    # three_phase_line
    "vector_potential_z",
    "field_xy",
    "triangle_far_field_amplitude",
    "planar_far_field_amplitude",
    "helical_near_field_amplitude",
    "helical_far_field_amplitude",
    "triangle_positions",
    "planar_positions",
    "hexagon_positions",
    "balanced_three_phase_currents",
    "balanced_six_phase_currents",
    # elliptic_integrals
    "K_hastings_2",
    "E_hastings_2",
    "K_hastings_4",
    "E_hastings_4",
    # gauss_legendre
    "gauss_legendre_nodes_weights",
    "gauss_legendre_integrate",
    "gauss_legendre_integrate_2d",
]
