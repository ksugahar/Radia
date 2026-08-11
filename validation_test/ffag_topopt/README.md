# FFAG HDiv-MMM topology-optimization validation

`validation_ffag_cell_targets.py` builds the Bell--Abell non-scaling FFAG
soft-edge one-cell target family at seven proton energies from 31 to 250 MeV.
It checks the three prerequisites for each reduced closed orbit (cell bend,
periodic position, and periodic tangent), then checks the analytic map
transform and symplectic residual.

This lane establishes the momentum-indexed target contract.  The first magnet
PoC treats each supplied entrance-to-exit design orbit as fixed and reproduces
its one-pass transfer matrix; it does not solve a ring closure problem.
Periodic closed-orbit recovery remains a separate, harder FFAG-ring validation
lane.  Vacuum observation points are not an air volume mesh.

## Ferromagnetic DUCAS lineage

The Lego proposal is a direct extension of the ferromagnetic-shape method of
Murata, Abe, Ando, and Nakayama, *IEEE Transactions on Applied
Superconductivity* 19 (2009), DOI 10.1109/TASC.2009.2018104.  Their explicit
loop solves the error-field inverse problem with DUCAS/SVD, interprets the
recovered magnetizing current as an equivalent saturated-ferromagnet volume,
updates the material mesh, and repeats the field solution.

Radia keeps that continuous step as an auditable predictor.  Candidate columns
are complete positive-material HDiv-MMM responses, the band-normalized
ACA--QR--TSVD pseudo-inverse returns signed material fractions, and multiplying
each fraction by its candidate volume gives the equivalent material-volume
change.  Positive feasible fractions propose insertion; negative feasible
fractions propose removal.  A 0--1 LP projects the continuous predictor onto
whole HEX states, and an exact Schur/full active-system solve decides whether
the binary state is accepted.  Thus the DUCAS predictor introduces neither
gray material nor a design finite difference.  Each fixed-orbit validation
iteration records singular values, normalized modal field strengths, modal
material amplitudes, and the dominant equivalent-volume changes.

Transfer-matrix error is not sent directly to those material candidates.  On
each fixed design orbit, a separate small dense TSVD first identifies the
reachable field-to-map subspace.  A Chebyshev LP inside that retained subspace
minimizes the largest engineering-band error and produces a target correction
of the sampled binormal field and normal gradient.  The default relative TSVD
cutoff is 1e-3 so near-null optics modes do not create unbounded field steps.
Only that field target enters the HDiv-MMM material inverse.  Candidate
screening therefore performs no particle tracking or closed-orbit root solve.
After each accepted binary batch, the one-pass transfer matrices are
recomputed on the same frozen paths and the optics inverse is relinearized.
The validation history stores both TSVD layers independently.

The original method required a suitable stochastic initial shape because a
large continuous volume update could diverge.  This lane instead uses a fixed
return-yoke seed, a physical changed-volume trust region, connectivity gates,
and exact re-solves.  Nonlocal graph/beam exploration is a separate discrete
search layer; it does not alter the DUCAS material predictor.

The all-candidate ACA--QR--TSVD screen remains global.  The
`--exact-candidate-limit` option bounds only the costly conditional
block-Schur representative front after that screen.  A positive
`--exact-beam-width` and `--exact-beam-depth` retain fully resolved,
topology-valid states across a shallow objective barrier.  The returned design
is always the best fully solved incumbent, never an intermediate worsening
state.

`validation_ffag_hdiv_mmm_poc.py` is the next, study-scale lane.  It places a
structured BDM1 HEX candidate pole slab beside the unmeshed orbit aperture,
uses native Laplace field rows at several orbit energies, and requires at least
one exact whole-element move to improve the fused bend/map objective.  The
source and one-sided slab are intentionally simple so this remains a method
PoC; a final design needs the two poles, return yoke, and CoilBuilder geometry.

The saved study-scale result uses 192 HEX / 6,912 BDM1 DoFs and three momentum
points.  Two accepted exact re-solves remove 13 whole cells and reduce the
maximum normalized bend/map residual by 5.11 %.  The target bands are not yet
met; the result is evidence for the fused optimization path and its descent,
not a completed FFAG magnet.
