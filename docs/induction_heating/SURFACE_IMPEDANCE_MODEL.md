# User-designed surface impedance Zs(H, T)

`radia.panels.em_table.SurfaceImpedanceModel` accepts a user-designed material
law. It does not prescribe a steel curve or fit experimental measurements.

## Contract

- H: non-negative tangential magnetic-field **peak** amplitude, A/m.
- T: temperature in degrees Celsius.
- Zs: complex ohms, with finite components and Re(Zs) >= 0.
- Frequency: one explicitly declared positive frequency in Hz per model.
- Inputs broadcast as NumPy arrays. Return the broadcast shape, or a scalar.
- Outside the declared domain, evaluation raises an error. No silent clamping.
- Mean heat flux is `q = 0.5 * Re(Zs) * H**2` in W/m2.

```python
import numpy as np
from radia.panels.em_table import SurfaceImpedanceModel

# Synthetic demonstration only: replace this expression with your model.
def impedance(H, T):
    return 0.01 * (1 + T / 1000) / (1 + H / 10000) * (1 + 1j)

model = SurfaceImpedanceModel(
    evaluate=impedance,
    frequency_hz=150000,
    h_bounds=(0, 100000),
    temperature_bounds_c=(20, 800),
    name="synthetic-demonstration-v1",
)
H = np.array([0, 1000, 10000])
Zs = model.impedance(H, 20, frequency_hz=150000)
q = model.heat_flux(H, 20, frequency_hz=150000)
```

Existing `EMTable` data can be wrapped with
`SurfaceImpedanceModel.from_table(table, name="material-revision")`.
The table uses bilinear interpolation in log(H) and T, with at least two
strictly increasing nodes per axis. Unlike the older raw interpolation API,
this wrapper rejects out-of-range queries. Positive H nodes are required;
there is no implicit extension to H=0. Provide an explicit analytic model
when a zero-field limit is needed. Table heat is calculated from interpolated
Zs and the queried H, not independently interpolated from stored q_surf.

## Spatial field and conditional eddy-current updates

The constitutive law is evaluated locally on the workpiece surface:

`Zs(x,y,z) = model(|Ht(x,y,z)|, T(x,y,z))`.

Store its complex values at the surface integration points with their mesh,
region and quadrature mapping. Never compare only the mean impedance: two
fields with the same mean can yield different current and heating patterns.

The coupled solver must reuse its last **converged** eddy-current solution
when both the sampled Zs field and all other electromagnetic inputs are
unchanged. Do not rerun the eddy-current solve merely because a thermal time
step elapsed. In particular, unchanged temperature and excitation after
nonlinear convergence allow the electromagnetic state to be reused.
Temperature changes that leave Zs and other EM coefficients unchanged also
permit reuse. The thermal problem continues advancing with the saved heat
load; its own coefficients and boundary conditions are updated independently.

Reuse requires the same mesh/geometry, surface mapping, frequency, source
phasors, electromagnetic boundary conditions, volume material coefficients,
material-model revision and solver/accuracy settings. A changed current or
frequency invalidates solution reuse even if temperature is unchanged.
If only the right-hand side changes, reuse of a valid factorization or
preconditioner is a separate optimization; it does not make the old field a
valid solution. No automatic I-squared scaling for nonlinear materials.

During an ESIM iteration, changing Ht can change Zs even at constant T.
Reevaluate the constitutive law using the current field and continue the
nonlinear convergence check. Never label an unconverged iterate reusable.
After an actual solve, accept the field and its heat only after the shared
true-residual and nonlinear gates pass. Save the accepted state and record
whether the next EM evaluation reused it or solved anew, with the reason.

Initially use exact equality of the finite sampled complex Zs arrays and
identical input identities. Approximate reuse based on tolerances requires
separate field/power error validation; do not silently introduce a temperature
threshold. Array ordering changes and remeshing invalidate the cache.

This is the integration contract for the coupled solver. Material evaluation
is implemented; the solver-level reuse cache is not wired into the BEM or
Simulink drivers yet.

## Integration status

This interface validates material models independently of any solver. It does
**not yet enable** local ESIM in the weak PEEC+BEM solver or nonlinear behavior
in the native Simulink IH runtime. Those routes retain their current guards
until weighted surface assembly, reciprocal reaction, heat transfer and
coupled updates are verified together. Do not remove guards to run this model.
A thermal material update must refresh the electromagnetic solution when the
field changes; scaling a fixed heat pattern by I squared is not that update.

Experiments can supply a replacement expression or table later. Numerical
checks use synthetic models and independent reference problems; experimental
agreement is a separate validation claim. Save model revision, input table,
frequency, current convention and temperature units with each study.
