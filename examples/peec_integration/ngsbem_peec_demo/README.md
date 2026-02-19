# ngsbem_peec_demo — PEEC Impedance Extraction Notebook

Jupyter notebook demonstrating PEEC impedance extraction using ngsbem (NGSolve BEM).

## Notebook

**[peec_impedance.ipynb](peec_impedance.ipynb)** — Full PEEC Loop-Star workflow:

1. Surface mesh generation (Netgen OCC)
2. BEM matrix assembly (inductance L, potential coefficients P)
3. MQS impedance sweep (R + jωL)
4. Full Loop-Star impedance with self-resonance (Schur complement)
5. High-order convergence study (p = 0, 1, 2)

## Prerequisites

- NGSolve with ngsbem
- Radia PEEC module (`src/radia/ngbem_peec.py`)

## Running

```bash
jupyter notebook peec_impedance.ipynb
```

## See also

- Source module: `src/radia/ngbem_peec.py`
- Demo scripts: [`../ngbem/`](../ngbem/)
