# Element-impedance independent FEM comparison plan

Status: proposed before numerical execution. This case must not be described as
passing until both routes, refinements and every required gate have run. Review
this plan before admitting a heavy job or changing its thresholds.

## Scope and ownership

The temperature-BH GPU candidate is `266580d23fb7b799c943020a285efdc8d8cc97eb`.
It is an implementation reference, not the final integrated source or accepted
wheel. Record the final source, wheel and actual loaded MEX hashes before runs.
Native implementation, adapters and packaging retain their existing owner.
This directory owns the independent comparison driver and public case only.

There are three separate verdicts:

1. Independent Python ESIM/direct-equation parity for B(H,T), rho(T), element
   Zs, electromagnetic state and thermal increments. The existing CPU oracle
   shares BEM geometry assembly; it does not establish independent FEM accuracy.
2. This prescribed, heterogeneous, **same-boundary-law** BEM versus FEM test.
   Both routes solve anew; recalculating loss from a frozen field is excluded.
   The supplied complex Zs has no declared bulk-material equivalent. Passing
   this test does not validate ESIM, skin-resolved volume physics or B(H,T).
3. Exact installed GPU artifact and official Simulink/Toolkit acceptance,
   including model checks, save/close/reopen, repeat, errors and context cleanup.

Full nonlinear volumetric FEM and experimental accuracy remain unverified until
separately performed. No previous one-hole or uniform-Zs result substitutes for
this case. The GPU candidate supports closed flat-P1 genus-zero surfaces and
bounded complete systems (the selected block preconditioner allows at most 512
unknowns). The public temperature-BH example uses a surface thermal layer,
supplied coil impedance and a single-valued amplitude law, without hysteresis
history or minor-loop loss. Neither sample certifies private engineering cases.

## Geometry and prescribed data

`case.json` is authoritative. Its tetrahedron is an original elementary public
case, with four explicitly ordered outward triangles and an external closed
square filament. Indices in that file are one-based. Patch e follows triangle e;
every refined child inherits its parent's Zs. Never permute material rows or
replace the four complex impedances with an average. Save hashes of the case,
refined vertices/triangles, parent-patch map and all actual sample points.

Two uniform BEM refinements retain the exact planar tetrahedron: 16 and 64
triangles, respectively. The geometry is not rounded or changed between routes.
The coil is prescribed at 1 A peak and 1 kHz, exp(+i omega t). Do not infer an
experimental operating instruction or a bulk conductivity from these numbers.

FEM builds its air domain independently from the four original planar faces,
subtracting the tetrahedron from a sphere. Netgen/OCC is an explicit independent
mesh-generation alternative for this elementary comparison; no BEM assembly
or BEM fields are imported into FEM. Keep four distinct boundary labels and
`air`/`outer` labels; save and check each `.vol` before solving. Curve only the
outer sphere; the tetrahedron remains exactly planar. Use caller-owned NGSolve
TaskManager and sparse Cholesky, without a PARDISO fallback.

## Independent FEM equations and conventions

Use a scattered magnetic vector potential A_s in air, HCurl order 2, and
nu0=1/mu0. The incident A_i and H_i come from independent Gauss integration of
the declared filament segments. Set tangential A_s=0 on the outer sphere.
On the body, n_air points into the conductor, n_body=-n_air. Assemble

    integral_air nu0 curl(A_s).curl(v) + nu0*gauge A_s.v
      + integral_body (i*omega/Zs) A_s,t.v_t
    = integral_body [n_air cross H_i - (i*omega/Zs) A_i,t].v_t .

The boundary coefficient is piecewise complex and uses the actual patch Zs.
The total tangential field is H_t=(i*omega/Zs) n_body cross (A_i+A_s)_t.
The independent reaction impedance uses composite filament subsegments of
maximum length 2 mm and the declared Gauss order on each subsegment:
Z_reaction=(i*omega/I) integral_coil A_s.dl. Compare this reaction contribution,
not a large supplied coil resistance that can mask a body-coupling error.

For BEM use its declared physical map H_t=field_gradient*phi; verify the source
orientation and phasor convention before comparison. No fitted phase, sign,
amplitude or spatial registration is permitted. Surface power density is
q=0.5*Re(Zs)*sum(abs(H_t)^2), face power is q*area, and total power is its sum.
GPU heat output is face-integrated W, not W/m2; retain both quantities.

Compare complex vector fields at the same recorded physical quadrature points,
with positive area weights. FEM probes approach from the air side; halve their
offset to bound evaluation error. Compare vector differences before taking a
norm, preserving phase. No raw-DOF comparison is permitted.

## Predeclared acceptance gates

For scalar relative errors use abs(a-b)/max(abs(reference),absolute_floor).
For spatial errors use the area-weighted L2 norm of the difference divided by
max(the reference L2 norm,absolute_floor). H norms have units A and q norms
have units W/m after area weighting. The comparison-plan JSON written before
execution must include the quantities, norms, regions, normalizations, floors,
tolerances below, exact commands, versions and source/input identities.

| Quantity | Norm / region | Absolute floor | Limit |
| --- | --- | --- | --- |
| Total body power | scalar relative, all four patches | 1e-18 W | 2% |
| Complex body reaction impedance | complex magnitude relative, prescribed filament port | 1e-12 ohm | 3% |
| Surface q | area-weighted L2, all four patches | 1e-18 W/m | 10% |
| Complex tangential H | area-weighted vector L2, all four patches | 1e-12 A | 5% |
| Linear true residual | norm(b-Ax)/max(norm(b),1e-30), each route | 1e-30 in RHS units | 1e-6 |
| FEM mesh effect | relative power and reaction impedance, same radius | scalar floors above | 0.5% each |
| FEM exterior-domain effect | relative power and reaction impedance, same fine surface size | scalar floors above | 0.5% each |
| BEM refinement effect | relative power and reaction impedance, declared two levels | scalar floors above | 1% each |
| Filament quadrature effect | relative power and reaction impedance, orders 16/24 | scalar floors above | 0.1% each |
| Air-side probe offset effect | area-weighted H and q L2, offsets 1e-9/5e-10 m | spatial floors above | 0.1% each |
| EM power identity | abs(P-0.5*Re(Z_reaction)*I^2)/max(P,1e-18), each route | 1e-18 W | 1% |

The power budget assigns 0.5% to FEM mesh, 0.5% to exterior truncation, and
1% to BEM discretization. Quadrature and trace-offset checks must be below
0.1% so they cannot dominate these components. Spatial limits are separately
declared because differentiation and local loss amplify discretization error;
they cannot be replaced by total-power agreement. If these controls fail,
refine within the supported scope or report failure; never increase limits
after seeing results. The nonlinear coupled sample additionally retains its
predeclared CPU parity, constitutive, thermal-conservation and lifecycle gates.

## Student handoff and release checklist

- Identify the exact accepted version/artifact and actual MATLAB `which` paths.
  Older wheels without the IH payload are not a temperature-BH delivery.
- Include a runnable generic model and explicit row-to-triangle/region mapping,
  face areas and SI units. Explain nodal surface temperature, integrated heat,
  complex Zs observation, and the one-update time alignment of accepted Zs.
- Material edits must recompute initial permeability and the continuous Ms
  tail from B(H,T), use a covering Kelvin axis, and preserve geometry identity.
  Stop/restart and reassemble after geometry/operator changes.
- Supply repeatable baseline output, a bounded material-change exercise,
  out-of-range error/recovery steps and saved numerical/lifecycle evidence.
- Distinguish prescribed element-Zs input from Zs derived from nonlinear
  B(H,T). State which mode each example executes and whether a FEM comparison
  is same-boundary-law or volumetric; never average the material to match FEM.
- Keep private measurements, real-device geometry and personal handoffs outside
  the public repository. Public docs, commits and artifacts contain no names.
- Codex implements and validates. Reciprocal Claude review, integration, CI
  and publication are coordinated only through repository-management leads.
