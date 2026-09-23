# Comparator follow-up plan

Offline comparison only. Fine mixed Omega differs from numerical BDM2 finer
by 0.2374773% core RMS; reduced-A by 0.2597175%. BDM2 finer is not an exact
solution and the cross-run diagnostic is not a newly passed three-engine lane.
No runtime ratio or superiority claim follows from this calculation.

Next full run: checked finer mesh family, all three engines, HDiv BDM2,
FEM order2, mu_r1000, eight threads, Gram eps1e-14, common existing driver
build/solve/direct field timing. Reuse exact-point cache wrapper with its
preparation time included. Capture native HDiv residual, full observation
contract, auxiliary spaces, sampled memory and complete provenance.

Duration estimate before launch: prior fine mixed Omega 1421.54s with14774
FEM iron elements; finer23484 gives a linear estimate2259.6s (~37.7min).
Allow40-55min for all3, setup and nonlinear scaling of direct work. This is
an estimate, not measured performance. HDiv finer already measured57.28s.
Increase the diagnostic watchdog from2400 to5400s for this run; timeout remains
failure, never acceptance. Do not launch duplicates or interfere with CI.
No new remote computation or remote staging was performed during this audit.

After finer result, compare FEM increments and field spread before selecting
matched-error timing points; independent-host reproduction remains outstanding.
Recover/hash/commit all evidence, then remove dedicated remote staging.
