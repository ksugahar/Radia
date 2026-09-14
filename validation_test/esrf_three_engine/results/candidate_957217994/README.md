# ESRF6 candidate input evidence

Status: **HOLD**. These files identify the candidate and selected meshes; they
do not certify field agreement, quadrature convergence, or performance.

## Candidate

- Native build source: `95721799420952a446d44793fc1baf83fa6c92f9`.
- CI run: `34703027075`, successful build/import/native smoke/wheel production.
- Wheel SHA256: `a4c64938ac814fc12b6631d7ae1e2bddf5da01efd25238ddda052ea080ccafe7`.
- Native SHA256: `f99ec26e15d880431ed560c3f5cac7a537056e66f1bea0872186cfe32c2158a0`.
- The isolated Hibino venv passed `pip check`. `preflight.json` records the
  installed Python/native files checked against the actual wheel bytes.
- No native binary was copied into an existing installation.

## Which mesh was used?

The old nonconforming iron mesh had SHA256
`0d4449962b60757cba78e771a136b2805a4323e26a1ee6d7ae482751bb7e5a5f`.
Commit `caa4294ab1fc47df82dc4899c2fdf8fd71ce9f09` fixed the unmerged Cubit
partitions and replaced it with the 2408-HEX mesh `fdd13872...` recorded in
`case6_mesh_contract.json`. The selected input matches that post-fix artifact,
not the old 4768-HEX mesh with duplicated faces and hanging contacts.

The FEM input is a separate conforming TET generation route, SHA256
`dfc12b84...`, with 134686 elements. Its September 7 generation record references
the corrected iron asset manifest `b2bbe5db...`; the mesh bytes also equal the
September 4 local-gap artifact. The generation record, byte identity, and
present topology checks are separate evidence. They do not prove when those
bytes were originally generated or that their resolution is sufficient.

`mesh_revision_audit.json` records a fresh read-only check on Hibino:

- Both selected meshes: zero hanging facets and duplicated boundary-face pairs.
- FEM: 1393 Kelvin node pairs; maximum translation error `9.60e-15 m`.
- Periodic H1, orders 1 and 2: 1393 and 5566 constrained free DOFs.
- The **constant** trace squared-integral ratios differ from one by at most
  `3.6e-15`. This sanity check is not a proof for all nonconstant trace modes.

`fem_check.json` and `iron_check.json` are strict structural/label checks.
Topology, volume and area agreement do not prove geometric fidelity or field
convergence. The C-type N=6/12/24 gap family is a different geometry and must
not be substituted for this ESRF6 quadrupole.

## Input selection gate

Run the standard-library-only checker before launching this fixed-input lane:

```text
python validation_test/esrf_three_engine/check_case6_inputs.py --iron <selected-iron.vol> --fem <selected-fem.vol> --output <new-run>/input_identity.json
```

It rejects unregistered bytes regardless of filename. There is no force option.
A new mesh family needs a reviewed contract and new numerical evidence.
`input_identity.json` here was also measured during the running candidate audit;
it is not falsely labelled a pre-launch check. The two existing runs had already
checked their input hashes through their runner/preflight provenance.

Old scratch files are not deleted blindly because other tasks may reference
them. Their mere presence is not evidence that the current runner opened them.
