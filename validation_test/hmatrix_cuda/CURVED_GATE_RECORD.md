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

Every family and order on small generated public meshes: a 102-tetrahedron unit
sphere, linear and P2, and the two-hexahedron fixture. These inputs are job-owned,
not pinned, and the limits stay the ones in the script.

| Family, order | DoF | rhs | Coefficients | Physical | Independent residual |
|---|---|---|---|---|---|
| Affine tet, 0 | 258 | 1.2e-16 | 1.0e-11 | 1.1e-12 | 9.6e-12 |
| Affine tet, 1 | 774 | 1.3e-16 | 5.5e-11 | 1.6e-13 | 9.0e-12 |
| Affine tet, 2 | 2160 | 1.5e-16 | 5.0e-10 | 1.3e-12 | 9.6e-12 |
| Curved tet, 1 | 774 | 1.3e-16 | 5.7e-11 | 4.0e-13 | 9.5e-12 |
| Curved tet, 2 | 2160 | 1.5e-16 | 4.7e-10 | 3.7e-13 | 9.9e-12 |
| Hex, 1 | 68 | 1.7e-16 | 2.6e-11 | 5.2e-13 | 6.1e-12 |
| Hex, 2 | 207 | 2.1e-16 | 1.9e-10 | 7.3e-13 | 5.9e-12 |

RT0 on hexahedra is refused, as declared. Pyramid and mixed meshes do not run
through this entry.

## Curved mass rule

On curved geometry this gate is not a comparison with the CPU route as it stands. Its
reference is a direct solve with the CPU Gram and charge map and the device mass rule. The Gram
and charge map still come from the CPU code, and the rule is written out in the script, so a
device operator error or a drift of the device rule still fails. The gate cannot see a
regression in the CPU curved mass or a growth of the gap below. Each curved run reports that
gap as `curved_cpu_route_difference`; nothing gates it. Before this entry the CPU route was the
reference, so curved results are not comparable across that change.

On curved tetrahedra the matched reference also takes the device material mass
rule, `IntegrationRule(TET, 2*max(6, 2p+2))`. The CPU default mass of the curved
charge basis uses NGSolve's default rule, which under-integrates the rational
Piola mass there. Measured on a 102-tetrahedron P2 unit sphere at order 1, against
rule degree 30: 6.8e-2 for the default rule and 1.3e-3 for the device rule. With
the default mass in the reference, the device differs by 1.5e-3 at order 1 and
6.6e-3 at order 2 in coefficients, and 7.9e-6 and 2.3e-6 in the mean field. With
the device rule in the reference, it differs by 5.7e-11 and 4.7e-10.

This is a finding about the CPU linear route, which is unchanged here. It is
not evidence that either mass is accurate enough for an adopted curved gate.

## Correction to the record of `810256e5b`

That commit message says the unchanged limits were measuring the ACA tolerance.
That is wrong. A dense Gram gave identical errors, so ACA was not the cause; the
cause is the quadrature profile difference above. The dense comparison it added
stays. The history is kept as it is; this paragraph corrects the reading of it.
