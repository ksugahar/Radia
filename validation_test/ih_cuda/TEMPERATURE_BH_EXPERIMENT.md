# Temperature-dependent B-H experimental workflow

This research workflow separates measured material properties, electromagnetic
loss prediction and thermal calibration. The original public tetrahedron uses
synthetic material values; it is a software demonstration, not experimental
agreement or a validated material model.

## Material input and assumptions

Use SI units: temperature in K, H in A/m, B in T, resistivity in ohm m,
frequency in Hz and surface impedance in ohms. Record specimen preparation,
orientation, dimensions and uncertainty privately with the experimental data.
Use a closed magnetic-path specimen where practical, or explicitly account for
demagnetizing field when interpreting the applied field as internal H.

The current amplitude B(H,T) model is single valued. It does not represent
remanence, coercivity, minor-loop history or hysteresis-loop loss. Do not insert
an entire measured hysteresis loop into its monotone amplitude table. Retain
measured loops and their cycle loss separately; disagreement in loss may expose
this modeling limitation even if amplitude permeability agrees.

Temperature interpolation uses the same increasing H knots at every measured
temperature and element. Tables start at H=B=0, have differential permeability
at least mu0, and connect continuously to B=mu0*(H+Ms) at the last knot.
The first segment determines initial permeability. No temperature extrapolation
is allowed. Evaluate temperatures through an explicit element map, not by
assuming FE coefficients are nodal temperatures. Different elements may select
different measured material tables.

## Measurements

1. Measure B-H loops over the intended temperature, frequency and excitation
   amplitude ranges. Record synchronized drive current and sense voltage,
   winding turns, effective path length and cross-sectional area, acquisition
   rate, offsets, phase calibration, waveform and temperature stability.
2. Derive declared amplitude B-H data and independently retain loop area/cycle
   loss. Repeat increasing and decreasing temperature sequences to identify
   irreversible changes and repeatability. Include held-out conditions.
3. Measure resistivity versus temperature independently. Record thermal
   conductivity, volumetric heat capacity and their provenance; identify
   convection/radiation/contact losses using separate cooling experiments.
4. On the simple coupled specimen, record current, voltage, frequency, phase,
   active/reactive input power, coil temperature and specimen temperatures at
   documented locations. Calibrate temperature sensors/emissivity and separate
   coil losses from specimen heating. Preserve the raw synchronized traces.

## Validation sequence

First verify the material table and ESIM cell against an independent CPU route.
Then hold temperature fixed and compare impedance and specimen loss across
frequency and amplitude. Next compare thermal evolution using independently
measured heating or prescribed heat, before enabling electromagnetic feedback.
Finally run closed coupling: accepted T -> element material -> Zs -> complete
reaction-field solve -> coherent heat -> backward-Euler T update.

Choose experimental agreement limits from uncertainty and the intended use
before fitting. Fit using one subset and test held-out temperatures/amplitudes.
Report amplitude/phase, power and temperature errors separately. Numerical
convergence (linear true residual <=1e-6 and the tighter configured constitutive
criterion) does not establish experimental accuracy.

Check timestep and surface/cell mesh refinement separately. Save initial state,
drive, boundary conditions, material-table hash, source/binary hashes, residuals,
energy-balance errors and accepted state trajectory in JSON. Report failed
candidates as failures; retain the previous accepted temperature, heat and
electromagnetic state. Reset and repeated Simulink runs must be reproducible.

The public flat-P1 genus-zero surface and supplied coil impedance are bounded
research assumptions. A thin surface thermal layer is not a resolved bulk
temperature field. Validate any real specimen's geometry, skin-depth regime,
bulk heat transport and boundary conditions separately. Experimental data and
real-device models remain in authorized private storage unless explicitly
cleared for publication.
