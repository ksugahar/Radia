# ESIM deep-skin slab adoption and surrogate test

The adopted evidence is recorded in:
- [Impedance convergence](../../validation_records/induction_heating/20261023_esim_slab_impedance/record.json).
- [Complex field convergence and phase](../../validation_records/induction_heating/20261023_esim_slab_phase/record.json).

These records use a constant-real-permeability, constant-conductivity slab with
half-depth 5 mm, 222 kHz, relative permeability 100 and conductivity 5e6 S/m.
The skin depth is about 48 micrometres, so a/delta is about 105: this tests
linear diffusion phase/sign and the deep-skin limit, **not finite-thickness or
midplane/symmetric-slab behaviour**. It does not certify nonlinear BH, thermal
coupling, full-coil models, GPU or independent FEM.

The fixed smooth-grid ladder uses N=25,49,97,193,385. Each adjacent pair must
have strictly decreasing error and observed q>=1.5, using actual h_max divided
by fixed half-depth. The finest relative error must be <=0.003. The field
quantity is complex L2 over the half-depth, integrating the piecewise-linear
interpolant element by element with high-order quadrature; it is not a nodal
error. Unchanged h_max or roundoff plateaus are unresolved. The adopted field
reference is H0*cosh(k*x)/cosh(k*a), x from centre and k=sqrt(j*omega*mu*sigma),
with exp(+j*omega*t). The impedance reference is rho*k*tanh(k*a).

h_max lies in a region where H is negligible. The rate is interpretable here
because the fixed smooth grading refines everywhere; the records also report
h_min/delta. The separate adaptive tolerance sequence is reference-only: its
h_max is unchanged as it preferentially refines the active skin layer, so no
rate is manufactured from those pairs.

The N25 surrogate test compares Re(Z_s) and Im(H1)/H0 with saved historical
values. Its roughly 2% field error against the analytic reference is not a test
failure and is not a reason to enlarge it. It asserts the diffusion-phase sign
and names its dependence on mesher grading in failure messages. Timings are
metadata only. Both a conductance x1.001 perturbation and sign reversal were
killed by the phase validation and the candidate historical phase sensor;
the previous Re(Z)-only guard could not detect sign reversal. Baseline and
surrogate-tolerance changes require reviewed regeneration after validation.
Validation is trigger-driven, not scheduled by these records.
