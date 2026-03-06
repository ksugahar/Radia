# Electromagnet Example (mu=1000)

C-Type electromagnet model for IMA (Image Method of Analysis) verification against ELF.

## Model Description

- **Material**: Linear soft iron, mu_r = 1000
- **Mesh**: 52 hexahedral elements (full model)
- **Symmetry variants**: full, x-mirror, z-mirror, quarter
- **Current**: 2000 AT racetrack coil

## Verification Results

All models verified against ELF reference data (2026-02-01).

### Matrix Comparison

| Model | Elements | DOF | Max Diff | Rel Error | Status |
|-------|----------|-----|----------|-----------|--------|
| full | 52 | 312 | 7.99e-07 | 7.77e-06% | PASS |
| x-mirror | 26 | 156 | 7.97e-07 | 7.77e-06% | PASS |
| z-mirror | 26 | 156 | 8.01e-07 | 7.89e-06% | PASS |
| quarter | 13 | 78 | 8.04e-07 | 8.50e-06% | PASS |

### Field Comparison at (0,0,0)

| Model | ELF Bz (mT) | Radia Bz (mT) | Diff (%) | Status |
|-------|-------------|---------------|----------|--------|
| full | -228.12 | -227.85 | 0.12% | PASS |
| x-mirror | -228.12 | -227.85 | 0.12% | PASS |
| z-mirror | -228.12 | -227.85 | 0.12% | PASS |
| quarter | -228.12 | -227.85 | 0.12% | PASS |

## File Structure

```
mu=1000/
  README.md              # This file
  coil_model.py          # Racetrack coil model (common)
  full/
    verify_elf_radia.py  # Full model verification
  x-mirror/
    verify_elf_radia.py  # X-mirror (MIMA X) verification
  z-mirror/
    verify_elf_radia.py  # Z-mirror (MIMA -Z) verification
  quarter/
    verify_elf_radia.py  # Quarter model (MIMA X + MIMA -Z) verification
```

## Running Verification

```bash
# Full model
cd full && python verify_elf_radia.py

# X-mirror (symmetric)
cd x-mirror && python verify_elf_radia.py

# Z-mirror (antisymmetric)
cd z-mirror && python verify_elf_radia.py

# Quarter model
cd quarter && python verify_elf_radia.py
```

## IMA Symmetry Settings

| Model | Radia image= | ELF MIMA | Description |
|-------|--------------|----------|-------------|
| full | None | (none) | No symmetry |
| x-mirror | '+x' | MIMA X | Symmetric on x=0 |
| z-mirror | '-z' | MIMA -Z | Antisymmetric on z=0 |
| quarter | '+x-z' | MIMA X, MIMA -Z | Combined |

## ELF Reference Data

Location: `S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\`

| Subdirectory | Description |
|--------------|-------------|
| full/ | Full model (52 elements) |
| x-mirror/ | X-mirror model (26 elements) |
| z-mirror/ | Z-mirror model (26 elements) |
| quater/ | Quarter model (13 elements) |

Files in each subdirectory:
- `ELF_magic.meg` - Mesh geometry (nodes and elements)
- `ELF_magic.mat` - Interaction matrix (Fortran binary)
- `ELF_magic.mag` - Solution (field at observation points)
