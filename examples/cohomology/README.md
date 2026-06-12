# Cohomology cut examples

Gmsh-free cohomology / homology **cut** demos for multiply-connected
magnetostatics, built on the pure-Python `radia.cohomology` engine (Kotiuga 1987 /
Bossavit / Pellikka 2013) — no Gmsh `computeHomology`, no `.msh → .vol` transfer.

A cohomology cut supplies the curl-free basis functions `h_k` (unit circulation
around the k-th current loop) that make the magnetic **scalar** potential
single-valued in a current-linking region, via the total-scalar T-Ω formulation

```
    H = -grad(phi) + sum_k NI_k h_k .
```

## Examples

| Script | What it shows |
|---|---|
| `tomega_wire.py` | A straight wire (current `I` along `z`) through an annular air region (`b1 = 1`). The cohomology cut carries the `NI` ampere-turns; the solved T-Ω field reproduces Ampère's wire field `H_φ(r) = I/(2πr)` and `∮ H·dl = I` (exact, by the unit-circulation cut). |

## Run

```bash
python tomega_wire.py
```

Requires `ngsolve` / `netgen` and the in-repo `radia.cohomology` /
`radia.cohomology_cut` engine. See the `cohomology-cuts` skill and
`src/radia/cohomology_cut.py` for the underlying solver.
