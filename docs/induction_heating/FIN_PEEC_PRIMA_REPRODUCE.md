# Reproduce fin-volume PEEC with PRIMA

The [fin beak STEP](../../tests/coil_from_cad/fixtures/fin_beak.step)
is a finned conductor validation geometry. Use a Radia CPU installation with the
native PEEC extension and STEP reader, and run the examples from the repository
root. The example excitations are illustrative. Branch currents, computed
heating and solver logs are generated locally.

## Compare unreduced and reduced coil solves

```sh
python -m radia.panels.calc_inductance --coil-solver fin-volume --coil-step tests/coil_from_cad/fixtures/fin_beak.step --coil-only --frequency 222000 --current 10 --fin-n-stations 25 --fin-n-lanes 24 --fin-n-layers 6
python -m radia.panels.calc_inductance --coil-solver fin-volume --coil-step tests/coil_from_cad/fixtures/fin_beak.step --coil-only --frequency 222000 --current 10 --fin-n-stations 25 --fin-n-lanes 24 --fin-n-layers 6 --fin-prima-order 67
```

The first command solves the full PEEC system. The second explicitly selects
PRIMA on the same lattice. Check the reported branch count, energy check and
full-order KCL/KVL residuals, then compare coil resistance and inductance.
PRIMA restores the full branch currents and raises if the existing residual
gate fails. The fin sweep is limited to 12,000 branches.

## Reuse a basis for four inputs

Use one terminal-current input and three Cartesian uniform-field inputs. The
example drives the coil with 10 A peak and a uniform field of 25 A/m peak per
axis, in phase. Use a nonzero terminal current. Set the field vector to zero for a coil-only drive on the four-input basis.

```python
from ngsolve import TaskManager
from radia.fin_volume import solve_fin_volume_prima_sweep

with TaskManager():
    sweep = solve_fin_volume_prima_sweep(
        "tests/coil_from_cad/fixtures/fin_beak.step",
        frequencies=(111000, 166500, 222000, 333000, 444000),
        expansion_frequency=222000,
        prima_order=67,
        sigma=5.8e7,
        current_A=10.0,
        uniform_field_A_per_m=(25.0, 25.0, 25.0),
        n_stations=25,
        n_lanes=24,
        n_layers=6,
        max_branches=12000,
    )
```

The result includes assembly and basis times and one solution per frequency.
The chosen PRIMA order is specific to this example. To include workpiece
heating, supply a checked genus-zero workpiece VOL with its boundary label and
material parameters to the weak scalar P1 SIBC route. The Python entry
`radia.panels.calc_inductance.run_fin_volume_prima_heating_sweep` takes the
usual induction arguments, the frequency sequence, the expansion frequency,
the real uniform field vector and an output directory. That entry includes the
uniform field in the direct workpiece excitation as well as the coil solve.

This geometry and workflow support reproducibility and comparison of the
reduced and unreduced calculations. **Absolute workpiece heating accuracy is
unqualified.**
