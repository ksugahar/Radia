# Beam Tracking Design Document

## Overview

This document describes the design for beam tracking in Radia, including:
1. Current B-field based Lorentz force integration
2. Future A-field based Hamiltonian formulation
3. GPU acceleration plans (2026 roadmap)

**Scope**: Magnetic field only (no electric field). Energy is conserved during tracking.

**Target Applications**:
- Undulators and wigglers (synchrotron light sources)
- Permanent magnet systems (insertion devices)
- Magnetic chicanes and dipole magnets

**Out of Scope**:
- RF cavity acceleration (requires electric field)
- Space charge effects
- Synchrotron radiation energy loss

For combined E+B tracking, use dedicated accelerator codes (elegant, OPAL, MAD-X, etc.).

## 1. Current Implementation (B-field based)

### 1.1 Physics

Relativistic Lorentz force equation:
```
dp/dt = q(v × B)
```

For a particle traveling primarily in the y-direction (longitudinal):
```
d²x/ds² = -(q/p) * [z'*By - (1+x'²)*Bz + x'*z'*Bx]
d²z/ds² =  (q/p) * [x'*By - (1+z'²)*Bx + x'*z'*Bz]
```

Where:
- `s` = longitudinal arc length (≈ y for small angles)
- `x' = dx/ds`, `z' = dz/ds` = transverse slopes
- `p` = particle momentum
- `q` = particle charge (electron: -e)

### 1.2 Implementation

File: `src/core/rad_particle_trajectory.cpp`

```cpp
class radTPrtclTrj {
    void TrjEqs(double s, double* F, double* dFds) {
        Field.P = TVector3d(F[0], s, F[2]);
        FldSrcPtr->B_genComp(&Field);  // Direct B-field calculation

        // Lorentz force equations
        dFds[0] = xd;
        dFds[1] = -Buf*(zd*rB.y - (1.+xdxd)*rB.z + xdzd*rB.x);
        dFds[2] = zd;
        dFds[3] = Buf*(xd*rB.y - (1.+zdzd)*rB.x + xdzd*rB.z);
    }
};
```

### 1.3 Integration Methods

| Method | Order | Adaptive | Use Case |
|--------|-------|----------|----------|
| `RungeKutta4` | 4 | No | Fast, fixed step |
| `RungeKuttaStep5` | 5 | Yes | High accuracy |
| `AutoPropagate` | 5 | Yes | Automatic step control |

### 1.4 Python API

```python
# Current API
rad.FldPtcTrj(obj, E_GeV, [x0, x'0, z0, z'0], [y0, y1], np)
# Returns: [[y, x, x', z, z'], ...] at np points
```

---

## 2. Future: A-field Based Formulation (Hamiltonian)

### 2.1 Motivation

For symplectic integration and GPU parallelization, the **Hamiltonian formulation** using vector potential A is preferred:

**Advantages**:
1. **Symplectic integrators**: Preserve phase space volume (important for long-term tracking)
2. **A is smoother than B**: Better convergence for numerical evaluation
3. **Natural for periodic systems**: Hamiltonian averaged over period gives tune/beta functions
4. **GPU parallelization**: A-field evaluation is embarrassingly parallel

### 2.2 Scope: Magnetic Field Only (No Electric Field)

**Design Decision**: The A-field formulation focuses on **pure magnetic field** environments (φ = 0).

**Rationale**:
- Target applications: undulators, wigglers, permanent magnet systems
- Electric fields (RF cavities, space charge) are handled by dedicated accelerator codes
- Simplifies implementation while covering primary Radia use cases

### 2.3 Hamiltonian Formulation (Magnetic Only)

Relativistic Hamiltonian with vector potential (no scalar potential):
```
H = c*sqrt(m²c² + (p - qA)²)

For ultra-relativistic particles (E >> mc²):
H ≈ |p - qA| ≈ p - qA_s + (px - qAx)²/(2p) + (pz - qAz)²/(2p)
```

Canonical equations of motion:
```
dx/ds = ∂H/∂px = (px - qAx)/p
dpx/ds = -∂H/∂x = q*∂As/∂x - q*(px-qAx)/p * ∂Ax/∂x - q*(pz-qAz)/p * ∂Az/∂x

dz/ds = ∂H/∂pz = (pz - qAz)/p
dpz/ds = -∂H/∂z = q*∂As/∂z - ...
```

**Note**: Without electric field, energy is conserved: dE/ds = 0.

### 2.4 Simplified Form (Paraxial Approximation)

For small transverse momenta (accelerator physics):
```
x' ≈ px/p - qAx/p
z' ≈ pz/p - qAz/p

dx'/ds ≈ (q/p) * [∂Ax/∂s - ∂As/∂x] = (q/p) * Bz
dz'/ds ≈ (q/p) * [∂Az/∂s - ∂As/∂z] = -(q/p) * Bx
```

This recovers the Lorentz force equation, but the A-based form enables symplectic integration.

### 2.5 Proposed API

```python
# New API (future)
tracker = rad.BeamTracker(obj, field_method='direct')

# B-field based (current, default)
traj = tracker.track_lorentz(E_GeV, [x0, x'0, z0, z'0], [s0, s1],
                              method='rk45', steps=1000)

# A-field based (Hamiltonian, magnetic only - no electric field)
traj = tracker.track_hamiltonian(E_GeV, [x0, px0, z0, pz0], [s0, s1],
                                  method='symplectic4', steps=1000)

# Batch tracking (GPU accelerated)
trajs = tracker.track_batch(E_GeV, initial_conditions, [s0, s1],
                            method='symplectic4', device='cuda')

# Note: Electric fields (RF cavities, etc.) are NOT supported.
# For combined E+B tracking, use dedicated accelerator codes (elegant, OPAL, etc.)
```

---

## 3. Symplectic Integrators

### 3.1 Why Symplectic?

Standard RK4/RK5 integrators do NOT preserve phase space volume. For long-term tracking (1000+ turns), this causes:
- Artificial damping/growth
- Incorrect tune calculation
- Wrong dynamic aperture

### 3.2 Proposed Integrators

| Method | Order | Stages | Accuracy | Speed |
|--------|-------|--------|----------|-------|
| Leapfrog (Verlet) | 2 | 2 | Low | Fast |
| Ruth 3rd order | 3 | 3 | Medium | Medium |
| Forest-Ruth 4th | 4 | 4 | High | Medium |
| Yoshida 6th | 6 | 8 | Very High | Slow |

### 3.3 Implementation Sketch

```cpp
// Symplectic 4th order (Forest-Ruth)
class SymplecticIntegrator4 {
    static constexpr double c1 = 0.6756035959798289;
    static constexpr double c2 = -0.1756035959798289;
    static constexpr double d1 = 1.3512071919596578;
    static constexpr double d2 = -1.7024143839193153;

    void step(double* q, double* p, double ds) {
        // q = (x, z), p = (px, pz)
        drift(q, p, c1 * ds);   // x += c1*ds * px/p
        kick(q, p, d1 * ds);    // px += d1*ds * (q*dAs/dx)
        drift(q, p, c2 * ds);
        kick(q, p, d2 * ds);
        drift(q, p, c2 * ds);
        kick(q, p, d1 * ds);
        drift(q, p, c1 * ds);
    }
};
```

---

## 4. Data Structures

### 4.1 Current

```cpp
class radTPrtclTrj {
    radTg3d* FldSrcPtr;  // Pointer to field source
    radTField Field;      // Field evaluation buffer
    // RK4/5 temporary arrays
};
```

### 4.2 Proposed

```cpp
namespace radia::beam {

// Field evaluator interface (magnetic only, no electric field)
class FieldEvaluator {
public:
    virtual void B(const Vec3d& r, Vec3d& B) = 0;
    virtual void A(const Vec3d& r, Vec3d& A) = 0;
    virtual void gradA(const Vec3d& r, Mat3x3& dA) = 0;
    // Note: No scalar potential phi - magnetic fields only
};

// Direct evaluator (current)
class DirectFieldEvaluator : public FieldEvaluator {
    radTg3d* source_;
};

// Particle state (magnetic field only - energy conserved)
struct ParticleState {
    double s;           // Longitudinal position
    double x, z;        // Transverse positions
    double px, pz;      // Canonical momenta (or x', z' for Lorentz)
    double E;           // Energy [GeV] - constant (no electric field)
};

// Tracker
class BeamTracker {
    std::unique_ptr<FieldEvaluator> field_;

public:
    // Lorentz force (B-based)
    std::vector<ParticleState> track_lorentz(
        const ParticleState& initial,
        double s_end, int steps,
        IntegratorType type = IntegratorType::RK45);

    // Hamiltonian (A-based, symplectic)
    std::vector<ParticleState> track_hamiltonian(
        const ParticleState& initial,
        double s_end, int steps,
        SymplecticOrder order = SymplecticOrder::Fourth);

    // Batch tracking (GPU)
    std::vector<std::vector<ParticleState>> track_batch(
        const std::vector<ParticleState>& particles,
        double s_end, int steps,
        Device device = Device::CPU);
};

} // namespace radia::beam
```

---

## 5. Migration Path

### 5.1 Phase 1: Refactor Current Code (2026 Q1)

1. Extract `FieldEvaluator` interface
2. Create `DirectFieldEvaluator` from existing code
3. Add A-field evaluation to `radTg3d`
4. Unit tests for B = curl(A)

### 5.2 Phase 2: Add Symplectic Integrators (2026 Q2)

1. Implement `SymplecticIntegrator` base class
2. Add Forest-Ruth 4th order
3. Add Yoshida 6th order for validation
4. Benchmark against RK45

### 5.3 Phase 3: GPU Acceleration (2026 Q3-Q4)

1. Implement GPU dipole kernel (CuPy RawKernel)
2. Implement GPU symplectic kernel
3. Multi-particle batch tracking
4. Performance optimization

---

## 6. Compatibility

### 6.1 Backward Compatibility

The existing API will be preserved:
```python
rad.FldPtcTrj(obj, E, initial, range, np)  # Still works
```

### 6.2 New API (Additive)

```python
# New module
from radia.beam import BeamTracker

tracker = BeamTracker(obj)
traj = tracker.track(...)
```

---

## 7. References

1. Forest & Ruth, "Fourth-order symplectic integration", Physica D 43 (1990)
2. Yoshida, "Construction of higher order symplectic integrators", Phys. Lett. A 150 (1990)
3. Wolski, "Beam Dynamics in High Energy Particle Accelerators"

---

**Last Updated**: 2026-01-14
**Author**: Claude Code (Radia Development)
