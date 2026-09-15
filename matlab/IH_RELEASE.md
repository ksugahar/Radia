# Radia IH Simulink Preview

## Validated scope and review limitations

The geometry assembler generates a frozen unit-current EM loss distribution;
runtime heat scales with current squared. It currently emits one Eddy unknown
and disables rotation. It does **not** generate temperature-dependent BH
operators or a nonlinear BH iteration. The generic small dense runtime and
hand-authored slope configurations are not evidence for those production features.
Slope assembly costs O(n_temperature * n_unknown^2), followed by an O(n_unknown^3)
dense solve when temperature changes; no large-scale performance claim is made.

`heat_projection` is a diagonal quadratic loss representation, not a general
FEM/BEM loss tensor. Overlapping basis coefficients require cross terms and must
not be supplied as if their squared magnitudes alone represented local loss.
The single frozen-mode assembler does not omit such cross terms: the full EM
solution is evaluated before its loss field is stored.

For periodic sample configurations, temperature state and incoming Thermal heat
are in material coordinates. Only Eddy maps to/from the stationary source frame;
Thermal never rotates its stored state. Fractional-angle linear interpolation
preserves a weighted integral after correction, not hotspot shape or peak.
Repeated interpolation is diffusive. Nearly cancelling integrals that cannot be
corrected safely raise an error instead of amplifying the field.

Dense flattened inputs are explicitly row-major ABI data. MATLAB callers should
pass a two-dimensional matrix to the validator or flatten with
`reshape(A.',1,[])`, never `A(:)`. A vector alone contains no recoverable layout
information; the current ABI cannot detect a caller's incorrect flattening.

`radia_ih.slx` is the first native induction-heating runtime preview. It contains
separate readable Level-2 MATLAB Eddy and Thermal S-Functions backed by checked
`radia_mex` object handles, explicit current/angle/ambient source blocks, a
visible temperature-feedback loop, and an `IH Parameters` block. Initialization
uses the configured Python installation only to locate NGSolve's shared runtime
DLLs. Eddy and Thermal steps execute in MEX and never call Python per step.
This package targets 64-bit Windows and is generated and gated with MATLAB and
Simulink R2026a.

After extracting the release, start MATLAB in the extracted directory and run:

```matlab
model = install_radia_ih();
```

The saved model contains a one-DOF diagnostic configuration so installation
can be checked immediately. It is not a physical design result. For an actual
case, set `IH Parameters > IH configuration MAT/JSON` to a file containing
`config` or `radia_ih_config`. Use `radia.simulink.makeIHNativeConfig` when
constructing the checked native configuration. That function requires
preassembled Eddy/Thermal operators and successful strict-label `.vol` reports.

The `Geometry Update` block provides `Browse...` controls for the workpiece
Netgen `.vol`/`.vol.gz` mesh and the coil `.step`/`.stp` (PEEC) or
`.vol`/`.vol.gz` (BEM-A) input. After
the assemble command and its output MAT/JSON configuration are set once,
selecting a replacement file is enough: the next diagram update or simulation
start detects path and content changes, rebuilds the operators, and reloads the
configuration. `Rebuild now` forces the same operation explicitly.

The repository now provides an explicit initialization-time operator assembly
route; its live scope and input contract are owned by the Radia IH MCP manual.
This archive alone does not certify an arbitrary CAD/VOL-to-heat workflow.
An unmocked production geometry acceptance remains separate from the verified
preassembled native runtime. LUT and lumped-state-space IH helpers are
not included in this package. Preparing a physical configuration currently also
requires the separately installed `check-vol` command for strict mesh-label
reports; the diagnostic model does not require it.

The thermal assembly supports order 1 and 2, including a separate 2D
axisymmetric thermal mesh with the physical `2*pi*r` weight. Order 2 uses
NGSolve H1 coefficients (P2 on triangles, Q2 on quadrilaterals), while the
surface heat handoff remains P1. Temperature ports in this mode carry FE
coefficients, not nodal Kelvin values. The monitor evaluates physical
temperatures at mapped quadrature points; its extrema are sampled, not a
certified global bound. The mean uses the FE integral. Periodic rotation of
this coefficient representation is unsupported and fails explicitly.

The discrete update order is fixed:

```text
Eddy at T(t)
  -> Eddy heat mapped into workpiece material coordinates
  -> Thermal update to T(t + dt), without rotating stored temperature
```

For a linear magnetic law, a current-only change rescales heat density without
refactorizing the Eddy operator. A temperature-dependent operator causes the
Eddy solve to update. Nonlinear BH iteration is not present in this preview and
`bh_mode="nonlinear"` fails during configuration instead of silently using a
linear solve.
