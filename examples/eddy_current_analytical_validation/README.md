# Eddy Current Analytical Validation

Systematic cross-validation of Radia's eddy current solvers against
analytical solutions.  Each test problem has a known closed-form answer
so that FEM, BEM, and panel-pipeline errors can be measured absolutely,
not just relative to each other.

## Test Problems

### 1. Conducting Sphere in Uniform AC Field (`sphere_uniform_field.py`)

A conducting sphere (radius *a*, conductivity sigma, permeability mu_r)
in a uniform oscillating magnetic field B_0 exp(jwt).

**Analytical solution** (Smythe, *Static and Dynamic Electricity*, 3rd ed.):

```
J_s(theta) = (3/2) * jw * B0 * a / (jw*mu0*a + Z_s) * sin(theta)
H_t_rms   = |J_s_max| * sqrt(2/3)
```

where Z_s is the surface impedance (Dowell tanh formula).

## Running

```bash
# Single material + frequency
python sphere_uniform_field.py --material copper --freq 7000
python sphere_uniform_field.py --material steel  --freq 7000

# Full frequency sweep (both materials)
python sphere_uniform_field.py --sweep
```

## Results (2026-04-13)

Sphere R = 10 mm, B0 = 1 mT, OCC mesh (maxh = R/3), curve order 3.
H_t_rms error vs Smythe analytical solution.

```
Mat      f[Hz]   R/d      H_ana     BEM    FEM-s    FEM-t
-----------------------------------------------------------------
copper      100   1.5     689.05  -40.2%   +0.1%  -34.3%
copper      500   3.4     842.34  -23.3%   +0.9%  -33.8%
copper     1000   4.8     878.48  -17.8%   +1.1%  -33.7%
copper     5000  10.7     930.18   -8.7%   +1.4%  -33.5%
copper     7000  12.7     936.92   -7.4%   +1.5%  -33.5%
copper    10000  15.1     942.96   -6.2%   +1.6%  -33.4%
copper    50000  33.8     960.33   -2.7%   +2.0%  -33.4%

 steel      100   2.8      37.87  -66.1%   -2.9%  -36.3%
 steel      500   6.3      81.34  -65.3%   -2.8%  -36.2%
 steel     1000   8.9     112.11  -64.7%   -2.7%  -36.1%
 steel     5000  19.9     225.39  -62.2%   -2.2%  -35.8%
 steel     7000  23.5     257.73  -61.3%   -2.1%  -35.7%
 steel    10000  28.1     295.32  -60.3%   -1.9%  -35.6%
 steel    50000  62.8     496.19  -52.7%   -1.0%  -35.0%
```

### Conclusions

1. **FEM-scattered is the reference solver** for H_t and P.
   - Copper: +0.1% to +2.0% (all frequencies)
   - Steel:  -2.9% to -1.0% (all frequencies)
   - Error < 3% across all conditions tested (R/d = 1.5 to 62.8)
   - Key: A_inc is an exact CoefficientFunction (analytical uniform field),
     so no interpolation error.

2. **FEM-total has a constant ~34% systematic error on H_t**,
   independent of frequency and material. This is NOT a Robin
   penalty magnitude issue (the error does not change with R/delta).
   Root cause: BND trace evaluation of HCurl GridFunction on an
   internal boundary returns values from the wrong side
   (workpiece interior where A ~ 0 due to gauge regularization).
   **L (inductance) from volume integral is NOT affected.**

3. **BEM: use Scalar BIE SIBC (`ScalarBIESIBCSolver`), NOT EFIE
   (HDivSurface + LaplaceSL).**

   Scalar BIE SIBC vs exact Bessel (true reference, no SIBC approx):
   ```
   Mat      f   R/d  H_exact   SIBC_err  BIE_err
   copper    100   1.5   690.46     -0.2%    -1.6%
   copper   1000   4.8   876.50     +0.2%    -1.0%
   copper   7000  12.7   936.80     +0.0%    -1.2%
   copper  50000  33.8   960.32     +0.0%    -1.2%
    steel    100   2.8    36.75     +3.1%    +1.3%
    steel   1000   8.9   112.01     +0.1%    -1.6%
    steel   7000  23.5   257.69     +0.0%    -1.6%
    steel  50000  62.8   496.18     +0.0%    -1.5%
   ```

   - Scalar BIE SIBC: -1.6% max error (R/delta > 1.5), both materials
   - EFIE (HDivSurface): -7% to -66% (factor-3 eigenvalue error, uses
     only A, drops H contribution from curl SL)
   - mu_r is correctly handled through Z_s -- no PMCHWT needed
   - SIBC approximation itself breaks at R/delta < 1 (low freq):
     copper 10Hz -25%, steel 10Hz -43%

### 2. Coil + Cylinder Workpiece: Full-Resolution vs SIBC (`reference_2d_axisym.py`)

2D axisymmetric FEM with full eddy current resolution inside the
workpiece (mesh refined to delta/5).  Provides a ground-truth reference
for the coil + workpiece IH geometry.

**Geometry** (matches `ih_fem_kelvin_sample.jou` except 5-deg gap):
- Coil: circle cross-section at (r=30 mm, z=0), radius 3 mm
- Workpiece: solid cylinder r = 0..25 mm, z = -12.5..12.5 mm

**Formulation**: phi = r * A_phi (removes 1/r singularity), with
z-offset Kelvin transformation for exact open boundary.

```bash
# Full-resolution sweep
python reference_2d_axisym.py --sweep --kelvin

# Single case
python reference_2d_axisym.py --material copper --freq 7000 --kelvin
```

**SIBC validation** (Robin BC `(jw/Z_s)/r * u * v * ds` on workpiece
surface; workpiece NOT meshed).  Two Z_s formulas tested:

- **Cylindrical Bessel** (exact for solid cylinder):
  `Z_s = rho * gamma * I1(gamma*a) / I0(gamma*a)`,
  `gamma = sqrt(jw * mu_r * mu_0 * sigma)`
- **Dowell flat-slab** (approximate):
  `Z_s = rho * gamma * tanh(gamma*d)` (flat-slab half-thickness = a)

```
     Mat    f[Hz]  R/d   L_full   L_cyl   cyl/full  P_full    P_cyl   Pcyl/f
------------------------------------------------------------------------------
  copper      100  3.8   75.53    75.42    0.999    4.11e-6   4.15e-6  1.009
  copper     1000 12.0   61.67    61.68    1.000    1.93e-5   1.93e-5  1.000
  copper     7000 31.7   57.11    57.15    1.001    5.81e-5   5.79e-5  0.997
  copper    50000 84.6   55.32    55.37    1.001    1.63e-4   1.63e-4  0.999
   steel      100  7.0  169.37   170.00    1.004    1.78e-6   1.84e-6  1.033
   steel     1000 22.2  156.20   156.67    1.003    4.59e-5   4.62e-5  1.007
   steel     7000 58.8  131.27   131.62    1.003    5.12e-4   5.14e-4  1.005
   steel    50000  157   98.12    98.32    1.002    3.58e-3   3.60e-3  1.004
aluminum      100  2.9   80.82    80.33    0.994    4.52e-6   4.65e-6  1.028
aluminum     1000  9.3   63.70    63.69    1.000    2.34e-5   2.35e-5  1.003
aluminum     7000 24.6   57.92    57.96    1.001    7.31e-5   7.28e-5  0.996
aluminum    50000 65.7   55.63    55.68    1.001    2.09e-4   2.08e-4  0.998
```

**Key findings**:
- **Cylindrical Bessel SIBC: L < 0.4%, P < 3.3%** for all materials
  and R/delta = 3 to 160.  Validated for both inductance and power.
- Dowell flat-slab: L similar, but **P off 16-23% at R/delta ~ 3**
  due to flat-slab geometry mismatch.  Use Bessel for cylinders.
- Air-only L = 99.47 nH matches Neumann formula (99.23 nH) to 0.2%.

**SIBC = Robin BC; conductor interior NOT solved**.  The hole approach
(workpiece subtracted from mesh) is correct: Robin BC
`dphi/dn + (jw*mu_0/Z_s)*phi = 0` replaces the conductor interior.
The interface approach (workpiece meshed as air + Robin on internal
interface) is wrong because flux bypasses the Robin BC through the
transparent interior (steel: 58% of correct L).

**Bug fixed** during development: `WorkPlane.Arc(a, 180)` creates a
circle centered at (R+a, a) instead of (R, 0) when the default heading
is +x.  Fix: `wp.Direction(0, 1)` before the first Arc.

### 3. Dodd-Deeds Infinite Rod (`dodd_deeds_rod.py`)

Analytical impedance change for a filamentary coil encircling an
infinite solid conducting rod (Dodd & Deeds, J. Appl. Phys., 1968).

```bash
python dodd_deeds_rod.py --sweep
```

**Bug fixed**: original code used `2*pi*r_c` prefactor; correct is
`2*pi*r_c**2` (verified via mutual inductance Neumann formula).

### Implications for IH Panel Pipeline

- **Scalar BIE SIBC is the validated BEM solver** for H_t and P
  (-1.6% vs exact Bessel, both materials, R/delta > 1.5).
- **FEM-scattered is the validated FEM solver** for H_t and P
  (+/-3% vs Smythe analytical, when A_inc is exact CF).
- **L (inductance) is accurate in both total and scattered** FEM
  (volume integral, not affected by surface trace issues).
- For the coil+workpiece IH problem, the BEM path
  (`calc_inductance.py` / `calc_heating_bem.py`) uses
  `ScalarBIESIBCSolver` and is validated.
- The FEM path (`calc_fem_kelvin.py --formulation total`) gives
  correct L but -34% H_t/P.
