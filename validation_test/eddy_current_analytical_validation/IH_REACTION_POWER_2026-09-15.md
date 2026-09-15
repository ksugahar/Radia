# IH reaction-power investigation: repository-owned reproducer

Status: **reproduced and narrowed; production correction and full acceptance remain open**.
No Takahashi artifacts are required. No installed solver was repointed or changed.

## Independent FEM baseline

`validate_axisym_power_balance.py` extends the existing axisymmetric full-conductor
FEM reference with a volume Joule integral and an optional 7.5 mm bore. It checks
the actual conductor volume against the analytic volume so a requested bore
cannot silently be omitted. The source is an impressed one-ampere peak current;
the source reaction and the volume loss are evaluated separately. Conductive
coil renormalization is explicitly outside this energy gate.

At 1 kHz, copper, radius 25 mm, height 25 mm, coil centre radius 30 mm and
coil cross-section radius 0.5 mm, both solid and bored workpieces pass at p=2,3:
maximum relative power imbalance is 2.8e-9. The solid p=3 loss is
1.9350330e-5 W. These are 2D axisymmetric solves, not 3D volume FEM solves.
Geometry is constructed in memory with OCC/Netgen for this explicit validation.
No VOL Save/reload/Curve workaround is used.

## BEM reproducer and missing term

`validate_sibc_reaction_power.py` uses the same solid cylinder and material,
with a 720-segment ideal source ring and planar half-space SIBC. It runs the
production scalar BIE kernel. The source is axisymmetric; the current BEM
implementation requires a 3D **surface** discretization. It does not perform a
3D volume FEM solve.

| BEM maxh | Surface loss / magnetic-only reaction | Complete reaction imbalance |
|---|---:|---:|
| 6 mm | 2.83076 | 1.0762% |
| 3 mm | 2.87224 | 0.4560% |
| 2 mm | 2.88524 | 0.1867% |

At 2 mm, surface loss is 1.9968017e-5 W. The existing magnetic-only reaction
term gives 6.9207394e-6 W. The electric reciprocity term contributes
1.3084556e-5 W; their sum is 2.0005295e-5 W. None of these reaction terms is
defined by reversing the surface loss into a resistance.

For peak phasors with exp(+i omega t), total magnetic potential phi and
incident potential psi, the scalar, single-valued, constant-Zs discretization is

    delta_L_m = integral(phi n.B_inc dS) / I^2
    delta_L_e = Zs/(i omega I^2) * psi^T K phi
    P_reaction = -omega/2 * Im(delta_L_m + delta_L_e) * I^2

There is **no complex conjugation** in this reciprocal bilinear pairing.
The electric term comes from the second term in the surface reciprocity
integral, E_inc cross H minus E cross H_inc. The existing phi.B helper
accounts for the first contribution, not both. The current kernel source hash
on hibino was checked against this worktree before the probe.

The two-term surface impedance formula is also given in equation 21 of
[Luo and Di Rienzo, High-order surface impedance boundary conditions in
three-dimensional boundary element modeling of eddy-current problems](https://www.sciencedirect.com/science/article/pii/S0955799726002420).
The derivation and numerical checks above, not a citation alone, are the basis
for the candidate correction.

## What this does not establish

- The received case's exact sixfold discrepancy is not yet explained in full.
- BEM surface loss is 3.19% above the finite-coil full-conductor FEM value.
  Ideal filament versus finite coil, finite FEM exterior, SIBC approximation,
  curved versus faceted geometry, and mesh convergence are not identical.
  This is an independent scale check, **not** an apples-to-apples acceptance.
- The bored case passes FEM energy balance, but loop-extended BEM has not been
  checked here. A multivalued potential needs its cut contribution; do not
  apply the single-valued phi.B formula unchanged.
- Spatially varying Zs requires a weighted surface form, not mean(Zs)*K.
- Nonlinear ESIM and spatial heat-map generation are not certified here.
- `calc_inductance.py` is not modified by this validation commit. Its loop path
  currently retains plain phi while replacing dissipation; that remains a
  separate consistency problem for reaction and local heat maps.

## Reproduction and acceptance remaining

Run in a foreground process on an idle compute host, inside caller-owned
NGSolve TaskManager (the scripts provide it):

    python validate_axisym_power_balance.py --coil-radius .0005 --output fem.json
    python validate_sibc_reaction_power.py --output bem.json

The stored `ih_fem_power_balance_20260915.json` and
`ih_sibc_reaction_power_20260915.json` include actual runtime and source hashes.
The BEM gate checks the complete reciprocal pairing at the finest mesh, not
agreement with a resistance back-calculated from surface loss.

Next: implement the complete reaction operator in the production source bridge,
verify identical-source FEM-SIBC versus BEM, then validate the full loop field,
variable/nonlinear impedance and local heat extraction before claiming closure.
