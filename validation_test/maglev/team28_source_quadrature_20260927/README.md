# TEAM 28: coil-source quadrature versus gap (2026-09-27)

Diagnostic evidence, not TEAM 28 acceptance.  development host, Radia 5.0.1 release
tree, NGSolve 6.2.2606, 6 threads.  `team28_source_quadrature.json` records
the script hash and the hashes of the two imported drivers
(`../team28_hcurl_vim_force.py`, `../team28_coilbuilder_eddy_bubble.py`).

Model: the 3-D HCurl eddy-bubble + VIM lane of `team28_hcurl_vim_force.py`
(p = 6 parent, `maxh = 0.025`, rank-3 reduction, 50 Hz, 20 A) with the
CoilBuilder solid winding packs evaluated by Radia.  The reduced system
(sampled resistance, tet-volume inductance) is assembled once on the
production current basis (`intorder = 10`).  The drive `int mode . A_ext`
and the force `int J x B_ext` are sampled on bases with `intorder` 14 and 18
built from the same mesh, parent space and ports; their response vectors
match the production ones to 1.3e-12, so only the source quadrature changes.
The disk moves towards the coil through the field-evaluation offset with the
mesh fixed.  Each case solves the reduced harmonic system
`(R + sL) c = -s P i` directly.  At the nominal gap the production sampling reproduces the
stored lane to 3.8e-5 (that lane was produced by release 4.95.48).

| gap | Fz (order 18) | order 10 relative error | order 14 |
|---|---|---|---|
| 10.8 mm | 1.1019 N | 2.0e-8 | 7.6e-10 |
| 5 mm | 2.7728 N | 8.9e-7 | 1.0e-7 |
| 2 mm | 4.3121 N | 1.2e-5 | 1.6e-6 |
| 1 mm | 4.9764 N | 3.1e-5 | 4.8e-6 |
| 0.5 mm | 5.3425 N | 5.0e-5 | 8.8e-6 |

Force sweep over gaps 1-3 mm in 0.05 mm steps: the residual after a quartic
fit is 4.2e-8 (order 10) and 2.1e-8 (order 18) of the peak force, and the
order-10 minus order-18 difference has a detrended oscillation of 1.1e-7.
The stiffness dF/dz differs by at most 2.3e-4 relative, a smooth bias.

Findings:

- The production source quadrature error grows as the gap closes but stays
  at 5e-5 at 0.5 mm; a near-element correction is not needed for this case.
- Axial motion of the disk shows no mesh-periodic force ripple; the
  element-to-quadrature layout does not change under this motion.  Lateral
  motion (a brake or a moving mover) changes that layout and is not measured.
- The lane sets `intorder = 10` explicitly.  The default `intorder` of the
  eddy-hybrid sampling helpers is 2, which would sample the source far more
  coarsely.

Reproduce from a checkout (`--repo` is the checkout holding
`validation_test/maglev`):

```text
python team28_source_quadrature.py --repo <checkout> --output team28_source_quadrature.json
```
