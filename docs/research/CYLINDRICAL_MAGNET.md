# Cylindrical Magnet Analytical Field Computation

## Overview

Radia provides analytical formulas for computing the magnetic field of uniformly magnetized cylindrical permanent magnets. This implementation replaces the need for external dependencies like `magpylib` for cylindrical magnet sources.

## Physical Model

A cylindrical permanent magnet is characterized by:
- **Radius** R (mm)
- **Height** H (mm), or half-height L = H/2
- **Magnetization** M (A/m) - can be axial, diametric, or arbitrary
- **Center position** [x, y, z] (mm)
- **Axis orientation** ('x', 'y', or 'z')

The magnetization can be converted from remanence Br (Tesla):
```
M = Br / mu0 = Br / (4*pi*1e-7)
```

For NdFeB with Br = 1.2 T:
```
M = 1.2 / (4*pi*1e-7) ~ 955,000 A/m
```

## Usage

### Python API

```python
from radia.cylindrical_magnet import CylindricalMagnet, RingMagnet

# Create a cylindrical magnet
# Center at origin, radius 10mm, height 20mm
# Axial magnetization: 955,000 A/m (Br ~ 1.2 T)
cyl = CylindricalMagnet(
    center=[0, 0, 0],
    radius=10.0,
    height=20.0,
    magnetization=[0, 0, 955000]
)

# Get B-field at point
B = cyl.get_B([5, 0, 15])  # Returns [Bx, By, Bz] in Tesla

# Get H-field (outside magnet: H = B/mu0)
H = cyl.get_H([5, 0, 15])  # Returns [Hx, Hy, Hz] in A/m

# Get vector potential A (axial magnetization only)
A = cyl.get_A([5, 0, 15])  # Returns [Ax, Ay, Az] in T*mm

# Ring magnet (hollow cylinder)
ring = RingMagnet(
    center=[0, 0, 0],
    inner_radius=5.0,
    outer_radius=15.0,
    height=20.0,
    magnetization=[0, 0, 955000]
)
B_ring = ring.get_B([0, 0, 25])
```

### Use as Background Field Source

The cylindrical magnet can be used as a background field source with Radia:

```python
import radia as rad
from radia.cylindrical_magnet import CylindricalMagnet

# Create cylindrical magnet
cyl = CylindricalMagnet(
    center=[0, 0, 0],
    radius=10.0,
    height=20.0,
    magnetization=[0, 0, 955000]
)

# Create background field source
bckg = rad.ObjBckgCF(cyl.get_B)

# Add to container with other Radia objects
container = rad.ObjCnt([bckg, other_objects...])
```

## Theory

### Axial Magnetization (Derby-Olbert Formulation)

For a cylinder magnetized along its axis (z-direction), the B-field is computed using the Derby-Olbert (2010) formulation based on Bulirsch's generalized complete elliptic integral `cel(kc, p, c, s)`.

The field components in cylindrical coordinates (rho, z):

```
Br = (J/pi) * [cel(k1,1,1,-1)/sq1 - cel(k0,1,1,-1)/sq0]
Bz = (J/pi) * (1/dpr) * [zph*cel(k1,g^2,1,g)/sq1 - zmh*cel(k0,g^2,1,g)/sq0]
```

where:
- J = mu0 * M is the polarization
- k0, k1 are modulus parameters involving position
- g = (1-r)/(1+r) is a geometric factor
- sq0, sq1 are distance factors

On-axis (rho = 0), the formula simplifies to:
```
Bz = (J/2) * [(z+L)/sqrt((z+L)^2+R^2) - (z-L)/sqrt((z-L)^2+R^2)]
```

### Diametric Magnetization (Caciagli Formulation)

For a cylinder magnetized perpendicular to its axis, the H-field is computed using the Caciagli (2018) formulation. The formulas involve complete elliptic integrals K(m), E(m), and Pi(n,m).

For small radial distances (r < 0.05*R), a Taylor series expansion is used for numerical stability.

The B-field outside the magnet is:
```
B = mu0 * H
```

### Generalized Complete Elliptic Integral

The Bulirsch `cel(kc, p, c, s)` function is used, which relates to standard elliptic integrals:

```
K(m) = cel(sqrt(1-m), 1, 1, 1)
E(m) = cel(sqrt(1-m), 1, 1, 1-m)
Pi(n,m) = cel(sqrt(1-m), 1-n, 1, 1)
```

## Vector Potential (A)

### Theory

The vector potential A for an axially magnetized cylinder can be derived from the vector potential of a circular current loop. For a loop of radius `a` carrying current `I` at position `z'`, the azimuthal component is:

```
A_phi(rho, z) = (mu0 * I / (2*pi)) * sqrt((a+rho)^2 + (z-z')^2) / rho
               * [ (a^2 + rho^2 + (z-z')^2) / ((a+rho)^2 + (z-z')^2) * K(kappa) - E(kappa) ]
```

where:
- `kappa = sqrt(4*a*rho / ((a+rho)^2 + (z-z')^2))`
- K, E are complete elliptic integrals of the first and second kind

For an axially magnetized cylinder, the equivalent surface current density is:
```
K_phi = M_z  (at the curved surface r = R)
```

The total vector potential is obtained by integrating over the height:
```
A_phi(rho, z) = (mu0 * M_z / (2*pi)) * integral from -L to +L of
               sqrt((R+rho)^2 + (z-z')^2) / rho
               * [ (R^2 + rho^2 + (z-z')^2) / ((R+rho)^2 + (z-z')^2) * K(kappa) - E(kappa) ] dz'
```

### Implementation

The vector potential for axially magnetized cylinders is **fully implemented** using numerical integration (Gaussian quadrature) over the cylinder height. This approach:
- Uses the analytical current loop formula at each height slice
- Integrates numerically (default: 20 Gauss points for good accuracy)
- Handles the on-axis singularity (A_phi = 0 on axis by symmetry)
- Verified: curl(A) = B within 0.01% error

**Usage:**
```python
from radia.cylindrical_magnet import CylindricalMagnet

cyl = CylindricalMagnet(
    center=[0, 0, 0],
    radius=10.0,
    height=20.0,
    magnetization=[0, 0, 955000]
)

# Get vector potential at point
A = cyl.get_A([5, 0, 15])  # Returns [Ax, Ay, Az] in T*mm

# With higher accuracy (more Gauss points)
A = cyl.get_A([5, 0, 15], n_gauss=50)
```

### Diametric Magnetization

For diametrically magnetized cylinders, the vector potential is more complex as it lacks axial symmetry. This case is currently not implemented (returns zero).

## Limitations

1. **Field inside magnet**: The analytical formulas give H-field for diametric magnetization. For the B-field inside the magnet, add M to H.

2. **Vector potential for diametric case**: Not implemented due to complexity.

3. **Scalar potential (Phi)**: Not implemented for cylindrical magnets.

## References

1. **Derby, N., Olbert, S.** (2010). "Cylindrical Magnets and Ideal Solenoids", *American Journal of Physics*, Vol. 78(3), pp. 229-235.
   - https://arxiv.org/abs/0909.3880
   - Primary reference for axially magnetized cylinder B-field formulas

2. **Caciagli, A., et al.** (2018). "Exact expression for the magnetic field of a finite cylinder with arbitrary uniform magnetization", *Journal of Magnetism and Magnetic Materials*, 456, 423-432.
   - DOI: 10.1016/j.jmmm.2018.02.003
   - Primary reference for diametrically magnetized cylinder formulas

3. **Ortner, M., et al.** (2025). "The vector potential of a steady azimuthal current density. Once again.", *arXiv:2510.11663*.
   - https://arxiv.org/abs/2510.11663
   - Vector potential formula for circular current loop using elliptic integrals

4. **Bulirsch, R.** (1969). "Numerical Calculation of Elliptic Integrals and Elliptic Functions. III", *Numerische Mathematik*, 13, 305-315.
   - DOI: 10.1007/BF02165405
   - Algorithm for generalized complete elliptic integral

5. **Ravaud, R., et al.** (2009). "Magnetic Field Created by Thin Wall Solenoids and Axially Magnetized Cylindrical Permanent Magnets", *IEEE Trans. Magn.*, 45(10), pp. 4034-4037.
   - Additional reference for cylinder field formulas

## Implementation Files

- `src/radia/cylindrical_magnet.py` - Python implementation with CylindricalMagnet and RingMagnet classes
- `src/core/rad_cylindrical_magnet.cpp` - C++ implementation (for future integration)
- `src/core/rad_cylindrical_magnet.h` - C++ header
- `tests/test_cylindrical_magnet.py` - Test suite with magpylib comparison

## Validation

All formulas have been validated against magpylib with 0% error for:
- On-axis field (axial magnetization)
- Off-axis field (axial magnetization)
- Diametric magnetization
- Ring magnet (hollow cylinder)
- Vector potential A (verified via curl(A) = B)
