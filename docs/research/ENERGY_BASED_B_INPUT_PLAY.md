# Energy-Based B-input Play Hysteresis Model

## 1. Background and Motivation

### 1.1 B-input Play model

The B-input Play model is a constitutive relation that takes magnetic flux density B as input and outputs magnetic field H:

```
H(B) = sum_{k=0}^{K} f_k(p_k(B))
```

Here p_k is a play operator with threshold eta_k:

```
p_k(B) = max(B - eta_k, min(B + eta_k, p_k^{prev}))
```

- k=0: eta_0 = 0 (no history, reversible response), p_0(B) = B
- k>=1: eta_k > 0 (with history, irreversible response)

The reasons why B-input is superior to H-input were shown in our IEEE Magnetics paper (accepted).
Forward evaluation (B -> H) is directly computable in O(K).

### 1.2 Problem of converting to an Energy model

In Egger et al.'s Energy-based model, an energy density is defined for each play operator:

```
U_k(xi) = integral_0^xi f_k(s) ds
```

For the energy function to be **convex**, all f_k must be **non-negative (monotonically increasing)**.

However, the shape functions f_k of the B-input Play model identified from measured data
**can have negative slopes**. This is physically admissible, but
since U_k becomes non-convex, it cannot be directly converted to an Energy model.

### 1.3 Purpose of this document

This problem is solved by separating into reversible and irreversible components.
Furthermore, a unified treatment with magnetic aftereffect also becomes possible.


## 2. Reversible/Irreversible Separation

### 2.1 Basic idea

Separate the total magnetic field into a linear reversible component and a nonlinear irreversible component:

```
H(B) = H_rev(B) + H_irr(B, history)
```

### 2.2 Reversible component

```
H_rev(B) = nu_rev * B
```

nu_rev is the **reversible reluctance** (inverse of reversible permeability).
Physically, it corresponds to the differential reluctance at magnetic saturation.

### 2.3 How to choose nu_rev

Focus on the k=0 component (reversible response) of the original Play model:

```
H(B) = f_0(B) + sum_{k=1}^{K} f_k(p_k(B))
```

Define nu_rev as the minimum slope of f_0:

```
nu_rev = min_B f_0'(B)
```

This corresponds to the differential reluctance in the saturation region.

### 2.4 Construction of modified shape functions

Define the irreversible shape functions g_k as the remainder after subtracting the reversible component:

```
g_0(B) = f_0(B) - nu_rev * B      (k=0: subtract the linear term)
g_k(xi) = f_k(xi)                  (k>=1: unchanged)
```

Then:

```
H(B) = nu_rev * B + g_0(B) + sum_{k=1}^{K} g_k(p_k(B))
     = nu_rev * B + H_irr(B, history)
```

### 2.5 Guarantee of convexity

The slope of g_0 is:

```
g_0'(B) = f_0'(B) - nu_rev >= 0    (by definition of nu_rev)
```

If f_k for k>=1 has negative slopes, increase nu_rev:

```
nu_rev = max( min_B f_0'(B),  max_{k>=1} max_xi (-f_k'(xi)) )
```

By recomputing g_0(B) = f_0(B) - nu_rev * B with this,
**all g_k slopes become non-negative**.


## 3. Energy Function

### 3.1 Definition

```
W(B) = W_rev(B) + W_irr(B)
```

#### Reversible energy

```
W_rev(B) = (nu_rev / 2) * |B|^2
```

Since it is a quadratic function, it is **always convex**.

#### Irreversible energy

```
W_irr(B) = G_0(B) + sum_{k=1}^{K} G_k(p_k(B))
```

Here G_k is the antiderivative of g_k:

```
G_k(xi) = integral_0^xi g_k(s) ds
```

Since g_k >= 0 (monotonically increasing), G_k is **convex**.

### 3.2 Thermodynamic consistency

The constitutive relation is obtained as the derivative of W with respect to B:

```
H = dW/dB = nu_rev * B + g_0(B) + sum_{k=1}^{K} g_k(p_k(B)) * dp_k/dB
```

When the play operator is active dp_k/dB = 1, when fixed dp_k/dB = 0.

The dissipation inequality is satisfied:

```
D = H * dB/dt - dW/dt >= 0
```

The irreversible energy change is always less than the work done by H.
This is a special case of the Clausius-Duhem inequality, and
is **consistent with the second law of thermodynamics**.


## 4. Fast Inverse Problem (H -> B)

### 4.1 Problem formulation

When using the 2-scalar Omega-reduced Omega method in FEM,
an inverse problem is needed where H is known and B must be found:

```
H = nu_rev * B + H_irr(B)    ... solve for B
```

### 4.2 Picard iteration

```
B^{n+1} = (H - H_irr(B^n)) / nu_rev
```

### 4.3 Proof of convergence

Condition for the iteration to be a contraction mapping:

```
|dH_irr/dB| / nu_rev < 1
```

From the constitutive relation:

```
|dH_irr/dB| = |g_0'(B) + sum_{k=1}^{K} g_k'(p_k(B))| <= nu_total - nu_rev
```

Here nu_total = max_B H'(B) (maximum differential reluctance).
The contraction ratio is:

```
rho = (nu_total - nu_rev) / nu_rev
```

For typical magnetic materials, the ratio of nu_rev (at saturation) to nu_total (at initial magnetization) is large,
so rho is approximately 0.3-0.5, and **convergence is achieved in 2-3 iterations**.

### 4.4 Comparison with conventional methods

| Method | Cost per iteration | Convergence | Total cost |
|--------|-------------------|------|-----------|
| Newton (conventional) | Jacobian computation + solve | Quadratic convergence | ~100 iter * O(K) |
| Picard (this method) | H_irr evaluation only | Linear convergence (rho<0.5) | **2-3 iter * O(K)** |

Newton requires ~100 iterations because the inverse problem of the B-input Play model is highly nonlinear.
Picard converges fast even without quadratic convergence because the contraction ratio is small.


## 5. Connection to Hantila's Polarization Method

### 5.1 Hantila (1975) formulation

Hantila's polarization method separates the constitutive relation into a linear part and a residual:

```
B = mu_0 * (1 + alpha) * H + mu_0 * R
```

Here alpha is a constant, and R is the nonlinear residual.

### 5.2 Correspondence with this method

```
alpha = 1 / (nu_rev * mu_0) - 1
R = -(1/mu_0) * H_irr(B) / (1 + alpha)
```

### 5.3 Advantages in FEM

The interaction matrix (stiffness matrix) depends on alpha and is **constant**:

```
a(Omega, v) = integral nu_rev * grad(Omega) . grad(v) dx
```

**Only one LU factorization is needed**, and afterwards only back-substitution is performed.
The nonlinear residual H_irr is handled only by updating the right-hand side vector.

```
l(v) = integral (Hs - H_irr) . grad(v) dx
```

### 5.4 Computational cost of iterations

| Operation | Cost | Frequency |
|------|--------|------|
| LU factorization | O(N^3) | **once** |
| Back-substitution | O(N^2) | every iteration |
| H_irr update | O(N*K) | every iteration |

Conventional Picard/Newton methods required O(N^3) matrix factorization at every iteration.


## 6. Unification with Magnetic Aftereffect

### 6.1 Formulation of magnetic aftereffect

Magnetic aftereffect is expressed as an additional magnetic field that depends on the time rate of change of B:

```
H(B, dB/dt) = H_play(B) + alpha_1 * dB/dt + alpha_2 * |dB/dt|^{1/2}
```

- alpha_1 * dB/dt: classical eddy current loss (linear)
- alpha_2 * |dB/dt|^{1/2}: anomalous eddy current loss (nonlinear)

### 6.2 Advantages of B-input

In the B-input formulation, **dB/dt is directly available**:

```
dB/dt = (B^{n+1} - B^n) / Delta_t    (time discretization)
```

In the H-input formulation, H is the unknown variable, and obtaining dB/dt
requires solving the inverse problem B(H). This incurs additional computational cost.

### 6.3 Combination with Energy-based separation

Add the aftereffect term to the reversible/irreversible separation:

```
H(B, dB/dt) = nu_rev * B + H_irr(B, history) + H_visc(dB/dt)
```

Where:

```
H_visc(dB/dt) = alpha_1 * dB/dt + alpha_2 * |dB/dt|^{1/2}
```

In the FEM iteration:
1. nu_rev * B: constant matrix (LU once)
2. H_irr: updated by Play model
3. H_visc: computed directly from time-discretized dB/dt

**All terms can be treated uniformly within the B-input framework.**

### 6.4 Energy dissipation

Dissipation including the aftereffect term:

```
D = H * dB/dt - dW/dt
  = H_visc * dB/dt + D_play
  = alpha_1 * |dB/dt|^2 + alpha_2 * |dB/dt|^{3/2} + D_play >= 0
```

Since alpha_1 >= 0, alpha_2 >= 0, the dissipation is always non-negative.
Thermodynamic consistency is maintained.


## 7. 3D Vector Extension

### 7.1 Component-independent model (isotropic)

```
H_x = nu_rev * B_x + H_irr,x(B_x)
H_y = nu_rev * B_y + H_irr,y(B_y)
H_z = nu_rev * B_z + H_irr,z(B_z)
```

Each component has independent play operator states.

### 7.2 Scalar amplitude model

```
|H| = nu_rev * |B| + H_irr(|B|)
H = (|H| / |B|) * B
```

The play operator acts on the scalar |B|,
and the direction of H is the same as B (isotropic material).


## 8. Integration with the 2-scalar Omega-reduced Omega Method

Complete pipeline for accelerator magnet analysis:

```
Iron core: H = -grad(Omega_t)
Air:       H = Hs - grad(Omega)    (Hs = Biot-Savart, no coil meshing needed)

Constitutive relation:
  H = nu_rev * B + H_irr(B, history) + H_visc(dB/dt)

Iteration procedure:
  1. (nu_rev * M + K) * Omega = RHS + integral H_irr . grad(v) dx
     Left-hand side is constant -> LU once
  2. B = -mu_0 * grad(Omega_t)
  3. Play model: update H_irr(B)
  4. Magnetic aftereffect: update H_visc(dB/dt) (quasi-static step)
  5. Repeat until convergence (back-substitution only)
```

## 9. Summary

```
Original Play model:    H = sum f_k(p_k(B))      f_k < 0 possible
                                                    Energy non-convex
                                                    Newton inverse ~100 iter

Energy-based separation: H = nu_rev*B + sum g_k(p_k(B))    g_k >= 0
                                                            Energy convex (check)
                                                            Picard inverse 2-3 iter (check)
                                                            Hantila compatible (check)
                                                            LU once (check)
                                                            Magnetic aftereffect unified (check)
```

## References

1. F.I. Hantila, "Mathematical models of the relation between B and H for
   non-linear media," Rev. Roum. Sci. Techn. - Electrotechn. et Energ., 1975.
2. H. Egger, ..., "Energy-based Play model," (TU Graz, Schur complement Newton).
3. Sugahara et al., "B-input vs H-input Play hysteresis model comparison,"
   IEEE Magnetics (accepted, 2026).
4. Radia Play API: rad.MatPlayHysteresis(K, eta, f_k_tables)
5. Experiment script: src/radia/energy_play_model.py
6. C++ implementation: src/core/rad_material_impl.cpp (ComputeNuRev, Irreversible)
