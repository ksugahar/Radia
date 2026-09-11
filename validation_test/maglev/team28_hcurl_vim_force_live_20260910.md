# TEAM28 Live Regression, 2026-09-10

The adjacent JSON is a new mdx2 solve, not a replay of the historical
acceptance record. It exercises the current Python changes against the HDiv
integration candidate's native kernel, not an installed release-wheel gate.

- Source baseline: `084ad5f00`, plus the runtime-safety changes committed with
  this record. SHA-256 identifies the exact executed source files below.
- Native source candidate: `966adc2e5` (native build predates the documentation
  merge; HDiv integration kernel unchanged).
- `_radia_pybind.pyd` SHA-256:
  `060e75808a1b4863aaa80976e7ca043031471f35672f834535ef9dc1d481e1f3`.
- Executed `_eddy_hybrid.py` SHA-256:
  `5b67bf0b1be833767427fbd6ac9d8c2682963f442ec48aacf9c1a9c8cce67267`.
- Executed `team28_hcurl_vim_force.py` SHA-256:
  `1728025a60674bdd63b8c0ddff40bb56349e94f612d7b61bebe33d30bbb80f80`.
- Python 3.12.10, NGSolve 6.2.2606, six NGSolve/OMP/MKL threads.
- Entry: `test_team28_hcurl_vim_force_live(Path('test-output'))`, with
  `PYTHONPATH` selecting the isolated `C:/temp/eddy-live-20260910` source tree.
  Installed editable packages were not changed.

Force magnitude error was 0.00512924, below the existing 0.01 limit. The
matrix-free reduced interaction had missing eigenvalue diagnostics; the runner
now probes its small live operator instead of relying on historical values.
The force and passivity checks all passed. This single mesh (`maxh=0.025`,
`outer_quad=4`) does not certify mesh/quadrature convergence or all geometries.

The Python distribution version reported in the JSON is source metadata, not
proof that this binary came from the corresponding published wheel.

## Remaining Review Scope

The epsilon-free claim still needs separation between analytic diagonal blocks
and sampled cross-interactions. Per-family projection tolerance, especially
pyramids versus polynomial affine families, also needs a separate numerical
change and validation. Neither issue is certified by this regression.
