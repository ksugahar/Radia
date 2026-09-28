# Induction-heating validation records

These JSON files retain the source commits, runtime versions, commands and
native hashes of their original runs. They are historical numerical records,
not acceptance results for the current working tree or the NGSolve 6.2.2607
migration. Command paths identify disposable run directories and are not
required installation locations.

The axisymmetric FEM/SIBC comparison does not establish convergence of every
surface loss component. For example, the fine 150 kHz result differs from the
axisymmetric reference by about 0.57% in total loss but 18% at the ends and
63% in the small bore contribution. Total-loss agreement must not be used to
certify local heat deposition near the sharp edges.

The through-hole rotor records exercise discrete rotor states and thermal
source transfer at their recorded settings. They do not establish convergence
for arbitrary angular spacing, geometry, temperature-dependent material, or a
different solver build. Production acceptance must refine the relevant mesh,
angular sampling and time step and retain the resulting power/energy checks.

Run `../test_fem_ams_parity.py` for the separate small complex HCurl comparison
of automatic AMS/BDDC+AMS selection against SparseCholesky at orders 1–3.
It checks true residuals, loss and inductance; it is not a performance benchmark
or a validation of the compound A–V formulation.
