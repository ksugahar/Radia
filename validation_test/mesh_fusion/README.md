# Mesh-fusion validation

The five `results_phase*.json` files record the numerical evidence presented by
`docs/mesh_fusion/mesh_fusion.ipynb`:

1. heterogeneous-mesh Poisson convergence,
2. Nitsche-mortar Poisson convergence and interface continuity,
3. two-dimensional electromagnetic Nitsche coupling,
4. harmonic-mortar Fourier orthogonality, and
5. accelerator shim/yoke field convergence.

The notebook recomputes and embeds its public results without writing docs-local
JSON. This directory owns the checked numerical records.
