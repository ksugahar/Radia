# TEAM28 Live Regression, 2026-10-03

The adjacent JSON is a fresh worker-b solve. It replaces the 2026-09-10 record,
whose stored keys had been renamed after the run instead of being re-run.

- Source: public main `d99efbc49` checkout; Python sources of `src/radia` are
  identical to the `v5.2.1` tag (`143a15e5a`).
- Runtime: the exact radia 5.2.1 candidate wheel (SHA-256
  `c35ae5d2d16da7fb4350a99002a52443c156323f8d5f689468ea99fa1287e8a8`,
  CI run 36997139820) installed in an isolated job venv, NGSolve 6.2.2607,
  Python 3.12.10. No shared installation was changed.
- Entry: `validation_test.maglev.team28_hcurl_vim_force.run([0.025], outer_quad=4)`,
  the same call as `test_team28_hcurl_vim_force_live`.

Force magnitude error 0.0051292, below the existing 0.01 limit; all checks,
including passivity of the HCurl-to-Foster handoff, passed. The 2026-09-10
run gave 0.0051292 and the same force to 3e-11 N. This single mesh
(`maxh=0.025`, `outer_quad=4`) does not certify mesh/quadrature convergence or
all geometries.
