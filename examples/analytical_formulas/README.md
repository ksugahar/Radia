# Analytical Formula Examples

These five short scripts exercise [`radia.analytical_formulas`](../../src/radia/analytical_formulas/),
which collects closed-form expressions taken from the IEE Japan review series

> 若尾真治, 五十嵐一, 藤原耕二, 野口聡, 松尾哲司, 亀有昭久,
> "Useful Formulas of Analytical Integration in Electromagnetic Field
> Computations (Part 1..5)", IEE Japan Joint Technical Meeting,
> 2002–2004.

The PDFs of the review live in [`to_developers/`](../../to_developers/);
each example below pins the relevant Part and equation numbers.

| Script | Topic | PDF reference |
|--------|-------|---------------|
| [`ellipsoid_demag_torque.py`](ellipsoid_demag_torque.py) | Demagnetization factor and torque of a rotational ellipsoid; sweep over aspect ratio. | Part 5 §5, eq 39–44 |
| [`ac_locus_demo.py`](ac_locus_demo.py) | Major / minor axis of the time-locus ellipse traced by an AC phasor (B, J, ...) — compares the closed-form against a brute-force time sweep. | Part 5 §4, eq 29–37 |
| [`cylinder_sphere_shielding.py`](cylinder_sphere_shielding.py) | Static shielding factor of cylindrical and spherical magnetic shells in a uniform external field; thin-shell asymptote check. | Part 1 §5, eq 23–24 |
| [`rect_magnet_2d_field.py`](rect_magnet_2d_field.py) | 2D field of a uniformly magnetised rectangular bar; verifies far-field 2D dipole law. | Part 2 §2, eq 2–3 |
| [`rectangular_plate_eddy.py`](rectangular_plate_eddy.py) | Eddy current in a thin rectangular plate under a slowly-varying perpendicular B-field; plots streamlines and cross-sections. | Part 1 §6.1, eq 26–27 |
| [`cross_validation_3d_vs_2d.py`](cross_validation_3d_vs_2d.py) | Cross-check the 2D bar formula against `radia.analytical_magnet.CuboidMagnet` in the long-bar limit (`L_z -> infty`); convergence is `O(1 / L_z**2)`. | Part 2 §2 + 3D Yang/Camacho |
| [`solenoid_axial_field.py`](solenoid_axial_field.py) | Central- and axial-field profile of a rectangular-section solenoid; verifies the long-coil limit `B_0 -> mu_0 J (a_2 - a_1)`. | Part 4 §4, eq 26–27 |
| [`three_phase_line_field.py`](three_phase_line_field.py) | Three-phase line far-field decay: triangle (1/r²), planar (1/r², √3 prefactor), helical (exp(-2π r/p)). | Part 4 §5, Part 5 §3 |
| [`elliptic_integrals_accuracy.py`](elliptic_integrals_accuracy.py) | Absolute-error plot of degree-2 vs degree-4 Hastings approximations of K(k), E(k) against `scipy.special`. | Part 3 §3, Tables 1–2 |
| [`gauss_legendre_demo.py`](gauss_legendre_demo.py) | Gauss-Legendre nodes / weights table and spectral-convergence plot for two smooth integrands. | Part 3 §4, Table 3 |
| [`cross_validation_solenoid_currentloop.py`](cross_validation_solenoid_currentloop.py) | Compare `solenoid_axial_field` (closed form) against a stacked set of `radia.analytical_magnet.CurrentLoop` rings; midpoint-rule converges as `O(h^2)` (4x error reduction per grid halving). | Part 4 §4 |
| [`eddy_current_complete.py`](eddy_current_complete.py) | Combined Part 6 / 8 / 9 extensions: skin depth + planar Z_s, cylinder AC impedance frequency sweep, thin-shell AC shielding, plate Joule dissipation, magnetic-shell internal-field decomposition, average B from cuboid magnet. | Part 6 / 8 / 9 |

All scripts run in seconds on a laptop, write a PNG next to themselves
when matplotlib is available, and print a tabular numerical summary
either way.
