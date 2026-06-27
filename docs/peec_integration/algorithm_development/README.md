# algorithm_development — Dowell Continued-Fraction Derivation

Scripts for deriving the continued-fraction expansion of the Dowell skin-effect formula.

## Files

| File | Description |
|------|-------------|
| `derive_dowell_cf_algorithm.py` | Main derivation: z·coth(z) continued-fraction coefficients |
| `derive_dowell_cf.py` | Continued-fraction expansion exploration |
| `derive_dowell_pade.py` | Padé approximant derivation |
| `derive_dowell_prima.py` | PRIMA connection to Dowell formula |

## Key Result

```
z·coth(z) = 1 + w/(3 + w/(5 + w/(7 + w/(9 + ...))))
w = z², τ = d²·μ·σ/2
```
