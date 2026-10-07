Per-panel ESIM acceleration is opt-in for the weak P1 genus-0 workpiece path:

```
--impedance-model esim --esim-per-panel --esim-panel-evaluator table
```

The default `direct` mode solves one finite cell per panel per outer iteration.
`table` constructs an
adaptive linear table in log field amplitude using the same finite-cell solver.
It checks real and imaginary midpoint errors, each normalized by impedance
magnitude, against `esim_tol/10`. These probes are empirical error controls;
they do not establish a global interpolation error bound.

Fields outside the current table trigger direct range extension and refinement,
with no extrapolation. The existing 0.001 A/m field floor is retained. A changed
cell configuration invalidates the evaluator; tables are never reused across
runs. Nonconverged/nonfinite/nonpassive cells and an exhausted construction
budget (`--esim-table-max-cells`, default 4096) fail loudly.

Every mode re-solves the electromagnetic field at the accepted impedance and
then **directly evaluates every panel** to certify the final constitutive
mismatch against `esim_tol`. Cache/table results never count as certification.
The final O(N) cell cost remains. `esim_panel_evaluation` separates outer,
table-construction, seed and certification calls/times, cache/interpolation
hits and the material/configuration hash. Cell execution is serial.

`benchmark_panel_esim.py` compares both modes on one self-authored cylinder
and analytic saturating material law. It checks nonuniform impedance, final
fields/power, constitutive certification, outer iterations and Anderson history.
Run repeated cold comparisons on an admitted idle compute host with identical
thread budgets; wall-time improvement is a measurement, not a CI pass criterion.
This stage does not extend genus-1 or strong-coupling support.

The [recorded comparison](results/panel_esim_acceleration_lab.json) used 240
panels and two repeats on LAB with one thread. Direct evaluation took
64.78/64.69 s and 3121 direct cell solves per run; table evaluation took
24.15/24.38 s and 1185 direct solves. Each mode retained 240 direct final
certification solves and 12 outer iterations. The largest table-versus-direct
relative Zs error was 1.03e-4 and relative power error 2.35e-5.

These are cold per-run tables/caches with assembly and field export included;
native imports and the shared mesh were prepared before the timed runs. The
prescribed filament source bypasses the coil solver. These timings apply to
this example and fixed runtime, with no process-startup or general speed claim.

Table evaluation remains opt-in: a finite construction budget can fail even
when direct evaluation is available. A small final constitutive residual is
not a general bound on nonlinear solution error or a proof of identical branch
selection. The recorded agreement verifies this example, with direct evaluation
retained as the reference for other material laws and solver controls.
