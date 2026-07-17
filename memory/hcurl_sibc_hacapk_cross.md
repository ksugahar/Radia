# HCurl-SIBC HACApK cross composition

Decision recorded 2026-07-18.

## Rejected construction

The first cross-only prototype put all volume, bridge, and SIBC quadrature
points in one scalar HACApK matrix and returned zero whenever target and source
belonged to the same basis partition.  The resulting checkerboard zero pattern
is not geometric.  With small leaves, ACA selected unreliable pivots: on a
780-point, 26-mode synthetic probe, tightening `aca_eps` from `1e-7` to `1e-11`
still left about two percent reduced cross error.

Do not revive a partition-zero scalar H-matrix as the production cross path.

## Production construction

Build the ordinary stable full sampled Laplace Gram and project it once with
the three current-component CSR maps.  For `cross_only=True` composition,
subtract each small sampled reduced diagonal block and add the selected exact
or BEM diagonal operator.  Cross blocks stay in one HACApK action, while the
reduced correction is cheap and ACA sees the physical kernel without an
artificial sparsity pattern.  The same probe reached relative cross errors of
`3.96e-8` at `aca_eps=1e-9` and `3.78e-10` at `aca_eps=1e-11`.

## Function-space boundary

This construction combines only vector-current bases inside HCurl-VIM:
volume EVRS, conductor-cycle bridge, and surface-Omega/SIBC.  HDiv-MMM uses a
different BDM magnetic-charge Gram.  Couple HDiv and HCurl with a separate
rectangular field operator; never claim they are blocks of one isomorphic
H-matrix.
