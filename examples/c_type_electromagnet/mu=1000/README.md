# mu_r = 1000 Case

Electromagnet simulation with linear iron (mu_r = 1000).

## Parameters

- **Material**: Linear iron, mu_r = 1000
- **Coil**: 1000 A, center at (0, 0, 0.05) m
- **Mesh**: 13 hexahedral elements (1/4 model)

## Results

| Point (m) | |B| (mT) |
|-----------|---------|
| (0, 0, 0) | 3.38 |
| (0, 0, 0.01) | 3.98 |
| (0, 0, 0.02) | 6.48 |

## Run

```bash
python run_simulation.py
```
