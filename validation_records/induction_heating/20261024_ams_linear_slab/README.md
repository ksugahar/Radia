# Linear volumetric slab: saved F1 validation evidence

This independently authored homogeneous slab uses a/delta = 2 and a boundary
driven analytic tangential vector-potential trace. It is unrelated to a private
engineering geometry. All inputs use SI units and peak phasors exp(+i omega t).

The six cases use refinement multipliers m = 2, 4 and p = 1, 2, 3. Lowest-order
HCurl uses AMS; full order-2/3 HCurl uses BDDC with the AMS wirebasket coarse
solver. Each iterative field is compared with the same-matrix sparsecholesky
field and the independent constant-coefficient analytic diffusion solution.

The executed runner SHA is e0fe40c98580ac87b3c89c65235b649ad6e6aa67a10c3478501b79da82e3d7c9.
The derived result identifies the original runner, manifest, result, gates,
terminal, wheel and native hashes. It preserves scientific values exactly and
excludes runtime paths, process identities, allocation details and operational
logs. `entry.py` audits the recorded gate arithmetic; it is not a FEM replay.

Approvals FROZEN-01, FROZEN-02 and FROZEN-03 are listed with full IDs and hashes
in result.json. FROZEN-03 makes p=1 finest errors recorded-only; its h-rate and
strict p-monotonicity still gate. The 1% finest-error criterion applies to p=2/3.
The actual largest iterative true residual is 1.0465108862687692e-12, below the
approved 1e-6 gate; internal GMRES target flags are retained without relabeling.

## Historical adoption

This record binds the actual radia5.3.3 source8a988152/native4b09e365/run04.
It is explicitly stale with respect to current-main qualification. It does not
backfill source hashes to current main. The original6-case
campaign is not replayed here. Re-certification belongs to the release wheel.

The surrogate covers m2p1 complex AMS and m2p2 full-HCurl BDDC+AMS.
Its rtol=1e-8/atol=0 is user-approved by APPROVAL-20261024-IH-AMS-F1-SURROGATE-TOL-01. Both production zero-output
mutations failed the unchanged residual contract and clean sources/native
were restored. Python tracing does not establish native internal coverage;
strict A4 rejects historical/current source and artifact identity differences.
The unchanged-method source comparison is informative, not recertification.

The reference verification is manual symbolic algebra, not a CAS receipt.
The official series checker covers per-component h/p rates; entry.py also
audits all102 saved gates, including algebraic parity and plateau rules.
Case runtimes were not recorded by run04 and are explicitly null. Baseline
runtime metadata comes from the separately identified surrogate measurement.

Scope excludes nonlinear BH, GPU, thermal/PEEC-BEM/fullcoil and other a/delta.

The current-main native build has no empty-vertex guard. The previously measured
m2 cases never triggered that guard; its removal changes neither case's numerical
path. Recorded repeatability and production-native mutations remain measurement
evidence, not a replay or re-certification of the historical campaign.
