# `validate_curved_hdiv.py`: what its comparison means

Status: **proposed, not adopted.** No record in `validation_records/` adopts this
gate. The input candidate of the 774-DoF run below is not pinned yet, so that run
cannot support an adoption.

## Two comparisons, not one

| | Original reference | Matched-profile reference |
|---|---|---|
| CPU reference settings | `intorder=4`, `ho_far_factor=1e30` | The device near/far profile (those two settings dropped) |
| What it tests | Independent reference on its own quadrature profile | Operator and implementation parity |
| Result on the checked coarse sphere | Failed the unchanged limits: coefficient error about 7e-8 | Accepted, see below |

The original failure stands at its source. It is not relabelled, and the later
pass does not turn it into one. The 7e-8 is the difference between two quadrature
profiles; it was not shown to be a device defect, and it is also not shown here
to be harmless to physical accuracy.

The matched-profile comparison shares the rule arithmetic with the route under
test. Agreement says that the device implements the same discrete operator as the
CPU route. It says nothing about whether that operator is accurate. Physical
accuracy needs an independent reference, such as the ideal sphere in
`validate_nonlinear_sphere.py` or the analytic cases of the CPU route.

The thresholds did not change, but an unchanged threshold does not authorize a
changed reference on an adopted gate. Adopting this comparison needs the user's
approval and a pinned input.

## Matched-profile run on record

Affine geometry, order 1, checked linear coarse sphere, 102 tetrahedra, 774 DoF:
rhs error 1.3e-16, coefficient error 5.5e-11, physical error 1.6e-13, independent
residual 9.0e-12, CPU residual 6.0e-12. Curved geometry and order 2 have not run
through this entry.

## Correction to the record of `810256e5b`

That commit message says the unchanged limits were measuring the ACA tolerance.
That is wrong. A dense Gram gave identical errors, so ACA was not the cause; the
cause is the quadrature profile difference above. The dense comparison it added
stays. The history is kept as it is; this paragraph corrects the reading of it.
