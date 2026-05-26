# Engineering applications

| Script | What it does |
|---|---|
| `pwm_transient_field.py` | PWM (30 kHz square, 5 periods) transient on cylinder ── port current vs analytic step-response superposition (rel-err $7\times 10^{-8}$); simultaneously reconstructs internal $H_z(r,t)$, $J_\varphi(r,t)$ vs independent FDM PDE (rel-err $7\times 10^{-4}$) |
| `pde_pwm_reference.py` | Independent reference: 1D radial finite-difference PDE solver (LSODA, 99 DOF) for PWM transient — ground truth for `pwm_transient_field.py` |
| `realistic_inputs.py` | Sine sweep + dead-time PWM driving inputs |
| `ih_workpiece_extract.py` | Induction-heating workpiece: NGSolve FEM frequency sweep → hierarchical Cauer extraction |
| `production_ngsolve_cube_kelvin.py` | Production NGSolve FEM cube + Kelvin transformation (full geometry, BDDC preconditioned) |
