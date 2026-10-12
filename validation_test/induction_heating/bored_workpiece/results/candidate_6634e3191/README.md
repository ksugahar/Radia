# Radia 5.3.1 independent one-hole IH acceptance

The original synthetic stepped-bore case passed all declared gates using
the exact tagged Python source and native payloads from the accepted wheel.
[Numerical record](results.json), [artifact identity](release-identity.json),
[focused tests](focused.xml), [source verification](source-verification.json).

- Three BEM meshes: 2470, 3778, 6232 nodes.
- Independent FEM meshes: 601029, 486921, 517256 DOFs.
- Finest BEM versus finest independent FEM power disagreement: 1.396357% (gate 2%).
- Projected-carrier power width: 0.118411%.
- Point-trace comparison width: 1.340039%.
- Cohomology, loop-work and HACApK focused regressions: 47 passed, no skips.

Both numerical routes assume linear surface impedance. This does not validate
a skin-resolved conducting-volume model, nonlinear material or general
off-axis/multiple-hole geometry. Timings are diagnostic, not performance acceptance.

The numerical JSON is unchanged from the executed record. Its `source_dirty`
flag includes an untracked native-payload manifest; tracked solver source is
clean at the release tag and all executed source hashes were checked.
The focused XML removes only identifying runtime metadata; its original is
retained privately with its hash. Physical values and acceptance criteria are unchanged.
