# GPU-Accelerated Multi-Particle Trajectory Plan

**Status**: Future (PEEC priority first)
**Created**: 2026-02-28
**Hardware**: Quadro RTX 5000 (16GB VRAM, 3072 CUDA cores, Turing), CUDA driver 13.0

## Background

The goal is multi-particle trajectory calculation (single-pass, not circular accelerator)
with >10,000 dipole elements.

## Current State

| Component | Status | Notes |
|-----------|--------|-------|
| Particle trajectory | Single-particle | RK4/5 adaptive, calls `ComputeFieldForTrajectory()` per step |
| Fld (batch) | Working | TaskManager parallelized |
| nvcc | Not installed | CUDA Toolkit not present on machine |
| CuPy/Numba | Not installed | PyTorch CPU-only installed |

## Stage 1: Multi-Particle Batch Trajectory (CPU)

**Goal**: Batch all particle positions per RK step for efficient field evaluation.

**Why this comes first**: The current single-particle code calls `Fld()` one point at a time.
Fixing the algorithm is more impactful than changing hardware.

### Design

```
Pre-solved Radia model (fixed dipoles after Solve())
  -> For each RK4 step:
       Collect all N_particle positions
       Fld(batch)(positions)  -> TaskManager parallelized
       Update all particles with B_all
```

### Key Properties
- Particles are independent (no inter-particle forces)
- Source (dipoles) is static -> computed once
- Python-level implementation sufficient (field eval dominates, not loop overhead)

### Implementation Sketch

```python
def multi_particle_trajectory(radia_obj, particles, ds, n_steps):
    """
    particles: list of (energy_GeV, x, dx_ds, z, dz_ds)
    ds: arc length step
    """
    for step in range(n_steps):
        # Collect all particle positions
        points = np.array([[p.x, 0.0, p.z] for p in particles])

        # Batch evaluation - TaskManager parallelized
        B_all = np.asarray(rad.Fld(radia_obj, 'b', points))

        # RK4 update each particle
        for i, p in enumerate(particles):
            Bx, By, Bz = B_all[i*3:(i+1)*3]
            p.rk4_step(ds, Bx, By, Bz)
```

### Files to Modify
- New: `src/radia/particle_trajectory_batch.py` (Python multi-particle driver)
- Existing: `rad_field_unified.cpp` (batch evaluation)

## Stage 2: GPU Acceleration of P2P Kernel

**Goal**: GPU-accelerate the near-field direct computation (P2P) for large element counts.

**Prerequisite**: Stage 1 must be complete (batch evaluation provides GPU-worthy workload).

### Recommended: CuPy RawKernel

- No nvcc installation needed
- `pip install cupy-cuda12x` only
- Write dipole kernel in CUDA syntax as Python string

```python
import cupy as cp

p2p_kernel = cp.RawKernel(r'''
extern "C" __global__
void dipole_p2p(const double* src, const double* mom,
                const double* tgt, double* field,
                int nsrc, int ntgt) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= ntgt) return;
    double fx=0, fy=0, fz=0;
    double tx=tgt[3*i], ty=tgt[3*i+1], tz=tgt[3*i+2];
    for (int j=0; j<nsrc; j++) {
        double rx=tx-src[3*j], ry=ty-src[3*j+1], rz=tz-src[3*j+2];
        double r2 = rx*rx + ry*ry + rz*rz;
        double r_inv = rsqrt(r2);
        double r3_inv = r_inv * r_inv * r_inv;
        double mx=mom[3*j], my=mom[3*j+1], mz=mom[3*j+2];
        double m_dot_r = mx*rx + my*ry + mz*rz;
        double c = 3.0 * m_dot_r * r3_inv * r_inv * r_inv;
        fx += c*rx - mx*r3_inv;
        fy += c*ry - my*r3_inv;
        fz += c*rz - mz*r3_inv;
    }
    field[3*i]   = fx * 1e-7;  // mu_0/(4*pi)
    field[3*i+1] = fy * 1e-7;
    field[3*i+2] = fz * 1e-7;
}
''', 'dipole_p2p')
```

## Performance Estimates

| Scenario | Current (1-particle seq.) | Stage 1 (batch CPU) | Stage 2 (+ GPU P2P) |
|----------|:---:|:---:|:---:|
| 1000 particles, 10K dipoles | ~1000x Fld calls/step | 1x batch call/step | 1x GPU kernel |
| Speedup vs current | 1x | ~100-1000x | ~500-5000x |

## References

- CuPy RawKernel: https://docs.cupy.dev/en/stable/user_guide/kernel.html
