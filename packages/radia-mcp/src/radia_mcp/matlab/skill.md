# Radia MATLAB
MathWorks' official MATLAB MCP Server is the execution substrate. The MATLAB
Agentic Toolkit owns generic MATLAB workflows and the Simulink Agentic Toolkit
owns generic model inspection, editing, checking, and testing. This
package adds the Radia/NGSolve MEX capability contract, a table-backed
Optuna-like `radia.optuna` layer with `SimulinkRunner`, and 43 generic ML/RL
gates backed by 86 self-contained MATLAB functions. The Simulink surface also
contains the validated 50 Hz TEAM 28 six-stage CLN force LUT and a numeric
HCurl Eddy Bubble/CLN reduced-state bridge. The bridge consumes trusted R/L/P
matrices, uses the existing `hybrid_vim.solve` MEX kernel for harmonic solves,
and builds a passive discrete state-space model; NGSolve remains the owner of
mesh assembly, high-order geometry, and DoF orientation. Acoustic FEM/BEM
remains separate.

For a reusable file exchange, call `radia.vim.ExportHCurlEddyCLNJSON` from
the NGSolve/Python side and load it with
`radia.simulink.loadHCurlEddyCLNModel`. An optional reduced force operator is
evaluated by `radia.simulink.evaluateHCurlEddyCLNForce`. The Team28 export
driver is `validation_test/maglev/team28_hcurl_vim_force.py --export-model`;
the first p=6 case becomes a MATLAB-readable exchange file. This is a
numeric reduced-model bridge, not a MATLAB reimplementation of NGSolve mesh
assembly.

The first high-level Python-free assembly path is
`radia.ngsolve.hcurl_eddy_cln_model`. It calls the MEX-native HCurl response
reduction on a `.vol` mesh, returns `M_r = V' M V`, `K_r = V' K V`, and
`P_r = V' ports`, and forms a local diffusion CLN model with user-supplied
conductivity and reluctivity. This is the MATLAB/Simulink route for the
high-order local FE projection; it is intentionally distinct from the full
HCurl-VIM Laplace/BEM inductance and rationalized SIBC path.

For moving coupling, use `radia.vim.ExportHCurlEddyCLNFamilyJSON` with a
strictly ordered height coordinate and a common reduced state basis. MATLAB
loads it with `radia.simulink.loadHCurlEddyCLNFamily`, evaluates the selected
snapshot through `radia.simulink.interpolateHCurlEddyCLNFamily`, and exposes
the result through `buildHCurlEddyCLNFamilyModel`. The default policy is
linear interpolation with an error outside the sampled height range; clamping
and PCHIP are explicit opt-ins. Python/NGSolve is therefore a preprocessing
dependency for assembling the family, not a runtime dependency for MATLAB or
Simulink.

Use `matlab_radia_mex_contract` to inspect the current C++ `radia_mex` command
list and the explicit non-parity boundaries. Use
`matlab_optuna_simulink_contract` for the MATLAB optimization and Simulink
workflow. The MEX boundary shares numerical contracts, not Python object
identity: NGSolve continues to own meshes, spaces, Piola maps, and FE
orientation.

For analytic electromagnetic sensitivities, use
`radia.topopt.optimizeAdjoint` with `Solver="mma"` or `Solver="sqp"`.
Objective and constraint gradients are mandatory and use a
design-by-constraint Jacobian. `radia.topopt.checkAdjointGradient` is an
explicit QA diagnostic, not a production finite-difference fallback. The
native HCurl material-topology path is
`radia.topopt.optimizeHCurlActivationAdjoint`; it consumes the existing MEX
complex-adjoint result and enforces the material-volume inequality in MATLAB.
