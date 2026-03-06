# GPU-Accelerated Multi-Particle Trajectory Plan

**Status**: Future (PEEC priority first)
**Created**: 2026-02-28
**Hardware**: Quadro RTX 5000 (16GB VRAM, 3072 CUDA cores, Turing), CUDA driver 13.0

## Background

Radia integrates ExaFMM-t for O(N log N) field evaluation via `FldBatch(method=1)`.
The goal is multi-particle trajectory calculation (single-pass, not circular accelerator)
with >10,000 dipole elements.

## Current State

| Component | Status | Notes |
|-----------|--------|-------|
| ExaFMM-t | CPU-only | OpenMP + AVX SIMD, no GPU in any ExaFMM variant |
| Particle trajectory | Single-particle | RK4/5 adaptive, calls `ComputeFieldForTrajectory()` per step |
| FldBatch | Working | OpenMP + FMM auto-switch (N>=1000) |
| nvcc | Not installed | CUDA Toolkit not present on machine |
| CuPy/Numba | Not installed | PyTorch CPU-only installed |

## Stage 1: Multi-Particle Batch Trajectory (CPU, ExaFMM)

**Goal**: Use ExaFMM's FMM tree properly by batching all particle positions per RK step.

**Why this comes first**: The current single-particle code calls `Fld()` one point at a time,
wasting FMM's O(N+M) advantage. Fixing the algorithm is more impactful than changing hardware.

### Design

```
Pre-solved Radia model (fixed dipoles after Solve())
  -> Build FMM source tree ONCE
  -> For each RK4 step:
       Collect all N_particle positions
       FldBatch(positions, method=FMM)  -> O(N_dipole + N_particle)
       Update all particles with B_all
```

### Key Properties
- Particles are independent (no inter-particle forces)
- Source tree (dipoles) is static -> build once, reuse every step
- Target positions change each step -> target tree rebuilt (cheap for clustered beam)
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

        # Batch FMM evaluation - O(N_dipole + N_particle)
        B_all = rad.FldBatch(radia_obj, 'b', points.flatten(), method=1)

        # RK4 update each particle
        for i, p in enumerate(particles):
            Bx, By, Bz = B_all[i*3:(i+1)*3]
            p.rk4_step(ds, Bx, By, Bz)
```

### Files to Modify
- New: `src/radia/particle_trajectory_batch.py` (Python multi-particle driver)
- Existing: `rad_field_unified.cpp` (optional: cache FMM source tree across FldBatch calls)

## Stage 2: GPU Acceleration of ExaFMM P2P Kernel

**Goal**: GPU-accelerate the near-field direct computation (P2P), which is 60-80% of FMM runtime.

**Prerequisite**: Stage 1 must be complete (batch evaluation provides GPU-worthy workload).

### Options

| Approach | nvcc Required | Effort | Performance |
|----------|:------------:|--------|-------------|
| **A. CUDA P2P kernel** | Yes | Medium | Best |
| **B. CuPy RawKernel** | No | Low | Good |
| **C. jaxFMM replacement** | No | High | Good (Laplace only, no dipole) |

### Recommended: Option B (CuPy RawKernel)

- No nvcc installation needed
- `pip install cupy-cuda12x` only
- Write P2P dipole kernel in CUDA syntax as Python string
- ExaFMM-t tree traversal stays on CPU
- P2P interaction list handed to GPU

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

### Integration Point

In `rad_exafmm.cpp`, the P2P direct computation loop (lines 155-200) would be
replaced by a GPU call when CuPy is available. The FMM tree (M2L, M2M, L2L)
stays on CPU since it's already fast.

## Performance Estimates

| Scenario | Current (1-particle seq.) | Stage 1 (batch FMM) | Stage 2 (+ GPU P2P) |
|----------|:---:|:---:|:---:|
| 1000 particles, 10K dipoles | ~1000x FMM calls/step | 1x FMM call/step | 1x FMM + GPU P2P |
| Speedup vs current | 1x | ~100-1000x | ~500-5000x |

## References

- ExaFMM-t: https://github.com/exafmm/exafmm-t
- CuPy RawKernel: https://docs.cupy.dev/en/stable/user_guide/kernel.html
- jaxFMM: https://pypi.org/project/jaxFMM/
- Yokota & Barba, "Treecode and FMM for N-Body Simulation with CUDA", GPU Computing Gems
