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

All scripts run in seconds on a laptop, write a PNG next to themselves
when matplotlib is available, and print a tabular numerical summary
either way.
