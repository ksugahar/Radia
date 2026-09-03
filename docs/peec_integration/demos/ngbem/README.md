# ngbem — NGSolve BEM-based PEEC Demos

Galerkin BEM-based PEEC solver demonstrations using ngsbem (NGSolve BEM).

## Files

| File | Description |
|------|-------------|
| `demo_ngbem_peec.py` | PEEC Loop-Star solver: matrix assembly, MQS/Full impedance, high-order convergence |
| `demo_ngbem_eddy.py` | FEM-BEM eddy current solver: Dirichlet BC / Calderon projector coupling |
| `demo_ngbem_coupled.py` | ngsbem coupled analysis demo |

## See also

- Public result notebook: [`../../public_demo.ipynb`](../../public_demo.ipynb)
- Public numerical showcase: [`../../peec_showcase.ipynb`](../../peec_showcase.ipynb)
- Validation modules: `validation_test/peec_integration/ngsbem_peec_demo/`
