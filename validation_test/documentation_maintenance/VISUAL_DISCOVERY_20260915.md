# Visual-first discovery: ESIM, complex coil and interface coupling

Owner: documentation maintenance. User request: demonstrate capability through
figures/WebGUI, with MCP remaining the operating manual.

## Delivered scope

- ESIM: three fresh 100 A weak-coupled PEEC/BEM cases on mdx2; 10, 50 and
  100 kHz converge in 8, 7 and 7 iterations. Solve wall times are 239.49,
  169.04 and 156.98 seconds (not a performance comparison).
- Ten ESIM scenes: coil/workpiece assembly CAD, workpiece surface mesh,
  three accepted tangential fields, three local-loss diagnostics, surface
  resistance and the separately labelled heat-transfer artifact.
- Complex coil: four scenes (actual CAD, sampling mesh, field magnitude,
  field vectors), using the exact finite-segment CF source. Nine interior
  probes agree with the direct finite-segment routine to relative 5.77e-16;
  this checks implementation consistency, not the centreline approximation
  against a thick conductor. Historical notebook code/outputs are preserved.
- Mesh fusion: three scenes (two-material mesh, coupled solution, absolute
  error). Reuses the original Nitsche solver definitions at h=0.1, p=2.
  The geometric interface is conforming: no nonmatching-mesh claim.

The ESIM notebook puts the visual sequence before its derivation. The other
two notebooks insert scenes near the front without rerunning their historical
campaigns. The ESIM main benchmark links to the new spatial notebook.

## Important scope boundaries and diagnostics

- ESIM `qsurf_sol` uses a power-normalized incident Biot–Savart distribution.
  It is not the same local field as 0.5 Re(Zs)|Ht|². At 50 kHz their integrals
  are respectively 18.6768 W and 15.9352 W. The transfer integral matches BIE
  power to about 5.3e-15 relative by construction, not independent validation.
- Whole P1 mesh CAD check at 1% failed because unused air/coil regions have
  5.42%/8.44% volume discrepancy. For this boundary-only PEEC/BEM showcase,
  the declared whole-file exploratory tolerance is 10%, with independent
  strict workpiece volume/sibc area checks below 1% (0.188%/0.144%). This is
  not permission to use that mesh as a volume-FEM acceptance fixture.
- No temperature solve, resolved-volume eddy-current comparison, strong
  coupling certification or mesh-convergence claim is added.
- The first mdx attempt lacked OCP. A disposable system-site-packages venv
  supplied CAD/notebook dependencies; no global install was changed.
- The initial multipole coil-field candidate failed the near-field probe
  comparison (relative gap about 2.70e10). It was rejected. The delivered
  notebook explicitly uses the existing exact `h_segments_cf` route instead;
  no solver kernel was changed or tolerance relaxed for this check.
- Notebook execution needed `clipping=None`; this Netgen version rejects
  `clipping={'enable': False}`. All final cells executed without saved errors.
- Windows Jupyter emits selector-thread / unencrypted-loopback warnings.
  They are not solver failures. No services were exposed deliberately.

## Verification and recovery

Focused suite: 28 tests passed, including the final canonical citation assertion.
Saved widget states/17 referenced scenes and input/result hashes are checked.
References are generated from the canonical parent bibliography; the nonlinear
method uses `hollaus2026nonlinear`, not the obsolete `Hollaus2025` spelling.

Rendered visual QA remains pending: the available browser rejected local
artifact access under its URL policy. No alternative browser/localhost
workaround was used. Saved state audit is not visual acceptance.

Recovered roots on LAB:

- C:/temp/radia_esim_recovered_20260915
- C:/temp/radia_visual_recovered_20260915

Retain these and their LAB preparation directories under the documentation
maintenance owner until rendered visual QA, then remove those exact paths.
Remote C:/temp/radia_esim_webgui_20260915 and
C:/temp/radia_visual_discovery_20260915 are pending cleanup after evidence
commit/hash verification. The former owns the shared disposable venv.
