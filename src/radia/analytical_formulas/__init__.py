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

Eddy-current and AC analysis:

plate_eddy           Plate eddy current T_z, J + total dissipation P
                     (Part 1 §6.1 + Part 6 §3)
shielding            Static + thin-shell AC shielding factor; spherical
                     shell internal field (Part 1 §5, Part 6 §2,
                     Part 8 §2)
conductor_impedance  Planar surface impedance and full Bessel cylinder
                     AC impedance (Part 6 §4 + §5)
sphere_eddy          Complex magnetic polarizability of a solid conducting
                     permeable sphere in a uniform AC field (BVP closed
                     form; Wait 1951, Landau-Lifshitz sec. 59)

Numerics and special functions:

elliptic_integrals   K(k), E(k) Hastings polynomial approximations
                     (Part 3 §3, Tables 1-2)
gauss_legendre       Gauss-Legendre nodes and weights and convenience
                     integrators on a general interval (Part 3 §4)
adaptive_quadrature  Adaptive Gauss-Patterson with node reuse
                     (Part 9 §2, n=0..3)

Field averaging:

cuboid_average_field  Average B over a target box from a uniform-M
                      source box (Part 6 §7; closed-form C++ kernel
                      from sympy-derived G1, G2 antiderivatives,
                      v4.22.0+; method="numerical" fallback kept)

induction_heating     Canonical AC eddy-current Joule-loss formulas
                      from textbook references (Smythe, Landau-
                      Lifshitz, Jackson) -- not the Stafl 1967
                      transcription, see commit history (v4.23.0+)
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
from .plate_eddy import plate_eddy_dissipation
from .shielding import (
    shielding_factor_sphere_thin_ac,
    shielding_factor_cylinder_thin_ac,
    spherical_shell_internal_field,
)
from .conductor_impedance import (
    skin_depth,
    planar_surface_impedance,
    cylinder_ac_impedance,
    cylinder_dc_resistance,
    cylinder_internal_inductance,
)
from .sphere_eddy import sphere_complex_polarizability
from .adaptive_quadrature import (
    patterson_nodes_weights,
    adaptive_integrate,
)
from .cuboid_average_field import average_B_in_box, average_demag_tensor
from .induction_heating import (
    cylinder_axial_eddy_loss,
    cylinder_axial_eddy_loss_small_ka,
    cylinder_axial_eddy_loss_thin_skin,
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
    # plate_eddy (extension)
    "plate_eddy_dissipation",
    # shielding (extensions)
    "shielding_factor_sphere_thin_ac",
    "shielding_factor_cylinder_thin_ac",
    "spherical_shell_internal_field",
    # conductor_impedance (Part 6 §4-§5)
    "skin_depth",
    "planar_surface_impedance",
    "cylinder_ac_impedance",
    "cylinder_dc_resistance",
    "cylinder_internal_inductance",
    # sphere_eddy (conducting permeable sphere in a uniform AC field)
    "sphere_complex_polarizability",
    # adaptive_quadrature (Part 9 §2)
    "patterson_nodes_weights",
    "adaptive_integrate",
    # cuboid_average_field (Part 6 §7 -- closed-form C++ kernel + numerical fallback)
    "average_B_in_box",
    "average_demag_tensor",
    # induction_heating (canonical AC cylinder eddy-loss formulas, v4.23.0+)
    "cylinder_axial_eddy_loss",
    "cylinder_axial_eddy_loss_small_ka",
    "cylinder_axial_eddy_loss_thin_skin",
]
