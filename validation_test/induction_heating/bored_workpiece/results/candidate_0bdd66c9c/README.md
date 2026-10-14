# Radia 5.3.2 independent one-hole IH acceptance

The original synthetic stepped-bore case passed all declared gates using the exact tagged Python source and native payloads from the accepted wheel.
[Numerical record](results.json), [artifact identity](release-identity.json), [focused tests](focused.xml), [source verification](source-verification.json).

- Three BEM meshes: 2470, 3778, 6232 nodes.
- Independent FEM meshes: 601029, 486921, 517256 DOFs.
- Finest BEM versus independently finest FEM power disagreement: 1.396357% (gate 2%).
- Projected-carrier power width: 0.118411%.
- Point-trace comparison width: 1.340039%.
- Cohomology, loop-work and HACApK regressions: 47 passed, no skips.

Both routes assume linear surface impedance. This does not validate skin-resolved conducting-volume FEM, nonlinear material or off-axis/multiple-hole geometry. Timings are diagnostic, not performance acceptance.

The numerical JSON is unchanged from the executed record. Its `source_dirty` flag includes an untracked native-payload manifest; tracked solver source is clean at the release tag and executed source hashes were checked. Focused XML removes only identifying runtime metadata; originals and hashes are retained privately. Physical values and acceptance criteria are unchanged.
