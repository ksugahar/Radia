---
name: cohomology-cuts
description: Compute magnetostatic cohomology cuts on checked Cubit/Netgen meshes with the native Radia de Rham engine.
---

# Native cohomology cuts for magnetostatics

The source mesh is a checked Cubit or Netgen `.vol` with a versioned label
contract. Validate that contract before initializing the solver. Gmsh is used
only for derived MSH 4.2 post-processing; it does not create or transfer solver
meshes, cuts or solver state.

A scalar potential alone cannot carry a nonzero circulation around a coil
hole. The formulation is

    H = -grad(phi) + sum_k NI_k h_k
    curl(h_k) = 0
    integral(loop_j, h_k . dl) = delta_jk

Here `NI_k` is ampere-turns. The number of independent generators is the first
Betti number of the field domain. A torus surface has two independent loops;
the field region must exclude current-carrying holes when computing cuts.

## Existing production route

`radia.cohomology.cohomology_basis` uses the combinatorial Hodge Laplacian
`L1 = G G^T + Curl^T Curl` on the NGSolve/Netgen mesh and returns native
lowest-order HCurl generators. `radia.cohomology_cut.CohomologyCutSolver`
consumes them through `setup_from_mesh(mesh, loop_circles=...)`.

After checking the `.vol` label contract:

```python
from ngsolve import Mesh, TaskManager
from radia.cohomology_cut import CohomologyCutSolver

mesh = Mesh("field.vol")
solver = CohomologyCutSolver()
with TaskManager():
    n_coils = solver.setup_from_mesh(mesh)
    # Check n_coils against the declared physical topology before excitation.
    solver.solve([100.0], {"iron": 1000.0})
    B = solver.get_B()
```

For multiple coils, declare `loop_circles` to associate unit circulations with
physical coil order. Never infer coil ordering from display entity tags.
Verify the discrete curl, unit circulations, topology and solver residual
using the owning validation cases. Do not change their acceptance criteria.
For a parameterized winding surface, analytic toroidal/poloidal generators
are the independent route; see the stream-function knowledge.

## References

- P. R. Kotiuga: cuts for magnetic scalar potentials.
- Z. Ren (2002): T-Omega dual/thick-cut formulation.
- M. Pellikka, S. Suuriniemi, L. Kettunen, C. Geuzaine (2013): homology and
  cohomology computation in finite element modeling, SIAM J. Sci. Comput.

Bibliographic provenance is distinct from an executable solver dependency.
