# NGBEM Diagnostics Prune Record

Date: 2026-06-27

The old `examples/ngbem_diagnostics/*.py` files were diagnostic iterations for
the low-frequency eddy-current BEM/FEM-BEM line. They were pruned from examples
instead of being kept as user-facing docs; recover the old source from git
history if a forensic comparison is ever needed.

Keep these lessons, not the scratch scripts:

- Vector FEM/BEM loss checks must separate basis orientation, loop projection,
  vector-potential normalization, and SIBC power normalization. Mixing those
  checks in one script made the old diagnostics hard to trust.
- Thin-skin conductor loss belongs on the SIBC path. A volume FEM/BEM solve is
  useful when the skin depth is mesh-resolvable or magnetic material behavior
  matters, but it should not be the default reference for thin-skin sweeps.
- Shielded PEEC coupling should be validated through a canonical solver/test
  surface, not by keeping many one-off compare and diagnose scripts under
  examples.
- Future runnable checks should land under `validation_test/` with explicit
  assertions and compact fixtures. Human-facing derivation belongs in polished
  solver docs, not source archives.
