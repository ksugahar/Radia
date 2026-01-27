# Nonlinear B-H Curve Simulation (20000 AT)

Electromagnet simulation with nonlinear iron (B-H curve).

## Reference

B-H curve data from: `S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/BHカーブ/iron.bh`

## Parameters

- **Material**: Nonlinear B-H curve (21 data points)
- **Coil**: 20000 AT (Ampere-turns)
- **Mesh**: 13 hexahedral elements (1/4 model)

## B-H Curve

| H (A/m) | B (T) |
|---------|-------|
| 0 | 0 |
| 82 | 1.14 |
| 898 | 1.59 |
| 4582 | 1.81 |
| 17736 | 2.01 |
| 68322 | 2.20 |
| 318000 | 2.56 |

## Results

| Point (m) | |B| (mT) |
|-----------|---------|
| (0, 0, 0) | 68.0 |
| (0, 0, 0.01) | 79.7 |
| (0, 0, 0.02) | 129.7 |

## Comparison with mu=1000

| Point | Nonlinear (20000 AT) | Linear mu=1000 (1000 A) | Ratio |
|-------|---------------------|-------------------------|-------|
| (0,0,0) | 68.0 mT | 3.4 mT | 20.0x |
| (0,0,0.01) | 79.7 mT | 4.0 mT | 19.9x |
| (0,0,0.02) | 129.7 mT | 6.5 mT | 20.0x |

Current ratio: 20000/1000 = 20x, field ratio ~20x (linear region)

## Run

```bash
python run_simulation.py
```
