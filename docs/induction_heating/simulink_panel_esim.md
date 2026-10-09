# Frozen panel ESIM at the Simulink geometry boundary

Geometry Update has three surface-impedance modes: `uniform`, `element-file`,
and `per-panel-esim`. Panel ESIM requires a B-H table (H in A/m, B in T), half
thickness in metres, peak reference current in amperes, outer relative tolerance
and iteration budget. The adaptive `table` evaluator is the default; `direct`
remains selectable. Both require final direct certification on every panel.

This is a **frozen constitutive response at the assembly reference current**.
The outer nonlinear solve runs at explicit geometry update or initialization;
it never runs at each simulation step. The native input `current_A` is the
peak phasor amplitude/envelope, not instantaneous samples of a carrier sine.
The reference phase is zero radians; frequency is the Geometry Update frequency.
The converged panel impedance stays fixed while native fields scale with the
current amplitude and heat scales with its square. This does not reproduce
nonlinear ESIM at arbitrary time-varying drive amplitudes.

The allowed drive band defaults to plus/minus 10% and can be tightened down to
zero; it cannot be widened beyond 10%. Outside that band, including zero drive, the Eddy block raises an error.
To use another operating amplitude, change the reference current and explicitly
rebuild/reinitialize. A drive band is a use restriction, not an accuracy guarantee.

The immutable B-H snapshot, ordered surface identity, frequency, iterations,
residual and direct certificate are bound to `converged-panel-zs.json` using
`radia.panel_surface_impedance.v1`. Fresh and reused updates verify its hash and
certificate; material, reference current or frequency changes invalidate reuse.
The named MATLAB entry point is `radia.simulink.assembleIHPanelESIMFromGeometry`.
Genus-1/strong coupling uses validated Python routes only; unsupported backend
combinations raise rather than falling back to uniform impedance.

For example, assemble at `esim_reference_current_A=2` and connect a Constant
with value `2` to Eddy input 1, labelled **peak current amplitude (envelope)**.
A slowly varying envelope between 1.8 and 2.2 A remains inside the default band.
Do not connect `2*sin(2*pi*f*t)` or a soft-start ramp from zero: zero/low samples
raise `IHESIMDriveBand` with an explicit envelope/zero-start error. There is no
startup grace interval. For another operating envelope, explicitly rebuild at
its reference amplitude before using it.
