function tests = test_python_fallback
% Verify the explicit in-process Python fallback contract and conversion.
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(testDir));
addpath(fullfile(repoRoot, "matlab"));
end

function testAcousticFallbackReportsBackendAndDll(testCase)
zeta = [0.1+0.2i, -0.4+0.3i];
result = radia.python.acoustic("cq", "bdf_delta", {zeta}, ...
    Keywords=struct("method", "BDF2"));
expected = 1.5 - 2*zeta + 0.5*zeta.^2;
verifyEqual(testCase, result.value, expected, RelTol=0, AbsTol=2e-15);
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.python.execution_mode, "InProcess");
verifyTrue(testCase, endsWith(lower(result.python.library), "python312.dll"));
end

function testFamilyWrapperRejectsWrongModule(testCase)
call = @() radia.python.acoustic("not_acoustic", "run", {});
verifyError(testCase, call, "MATLAB:validators:mustBeMember");
end

function testNumpyScalarConvertsToMatlabNumeric(testCase)
result = radia.internal.callPython("numpy", "float64", {1.25});
verifyEqual(testCase, result.value, 1.25, "AbsTol", 0);
verifyClass(testCase, result.value, "double");
end

function testIsochronousTopoptHasNamedFallback(testCase)
result = radia.python.isochronousTopopt( ...
    "density_to_s", {[0, 1], 100});
verifyEqual(testCase, result.value, [1e6, 0.01], "RelTol", 1e-14);
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.python.execution_mode, "InProcess");
end

function testBeamCurvilinearHasNamedFallback(testCase)
result = radia.python.beamCurvilinear( ...
    "bishop_rmf_frame", ...
    {[0, 0, 0; 0, 0, 1], [0, 0, 1; 0, 0, 1]});
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.beam_curvilinear");
verifyEqual(testCase, double(result.value.arc_length_m), [0, 1], ...
    "AbsTol", 0);
verifyEqual(testCase, double(result.value.horizontal), ...
    [1, 0, 0; 1, 0, 0], "AbsTol", 2e-15);
end

function testAcceleratorTaylorTopoptHasCubicNamedFallback(testCase)
result = radia.python.acceleratorTaylorTopopt( ...
    "third_order_taylor_map_from_multipoles", ...
    {[0, 0, 0, 0, 0, 0, 3], 0.08, 3}, ...
    Keywords=struct("maximum_step_m", 0.01));
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.accelerator_taylor_topopt");
transfer = result.value;
verifyEqual(testCase, size(double(transfer.R)), [6, 6]);
U = double(transfer.U);
verifyEqual(testCase, size(U), [6, 6, 6, 6]);
verifyGreaterThan(testCase, abs(U(2, 1, 1, 3)), 0.4);
end

function testAcceleratorTaylorTopoptHasSymplecticKANFallback(testCase)
result = radia.python.acceleratorTaylorTopopt( ...
    "Symplectic2x2KAN", {0.4, -0.2, 0.7});
verifyEqual(testCase, result.backend, "python-fallback");
block = double(result.value.matrix);
verifyEqual(testCase, size(block), [2, 2]);
verifyEqual(testCase, det(block), 1, "AbsTol", 2e-14);
end

function testAcceleratorLieTopoptHasCanonicalNamedFallback(testCase)
result = radia.python.acceleratorLieTopopt( ...
    "fourth_order_lie_map_from_multipoles", ...
    {[0, 0, 0, 0, 0, 0, 3], 0.08, 3}, ...
    Keywords=struct("maximum_step_m", 0.01));
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.accelerator_lie_topopt");
transfer = result.value;
verifyEqual(testCase, size(double(transfer.R)), [6, 6]);
verifyEqual(testCase, size(double(transfer.f3)), [6, 6, 6]);
verifyEqual(testCase, size(double(transfer.f4)), [6, 6, 6, 6]);
verifyEqual(testCase, size(double(transfer.V)), [6, 6, 6, 6, 6]);
verifyEqual(testCase, size(double(transfer.f5)), [6, 6, 6, 6, 6]);
verifyLessThan(testCase, ...
    double(transfer.factorization.relative_reconstruction_error), 1e-10);
verifyLessThan(testCase, ...
    double(transfer.factorization.reconstructed_symplectic_residual.maximum), ...
    1e-10);
end

function testAcceleratorLieTopoptAcceptsIndependentReferenceCurvature(testCase)
result = radia.python.acceleratorLieTopopt( ...
    "canonical_body_hamiltonian_jet", ...
    {[0.4, 0, 0, 0, 0, 0, 0], 3}, ...
    Keywords=struct( ...
        "reference_curvature_per_m", 0.03));
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.accelerator_lie_topopt");
jet = result.value;
verifyEqual(testCase, size(double(jet.H5)), [6, 6, 6, 6, 6]);
H2 = double(jet.H2);
fieldCurvature = 0.4 / 3;
verifyEqual(testCase, H2(1, 1), 0.03 * fieldCurvature, ...
    "AbsTol", 2e-15);
verifyEqual(testCase, H2(1, 6), -0.03, "AbsTol", 2e-15);
verifyEqual(testCase, H2(6, 1), -0.03, "AbsTol", 2e-15);
end

function testStreamFunctionHasNamedFallback(testCase)
result = radia.python.streamFunction( ...
    "abe_reduce_node_potential_scales", ...
    {eye(2), [1; 3]});
verifyEqual(testCase, result.value(:), [1/3; 1], "AbsTol", 1e-14);
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.stream_function");
end

function testTopologyFamiliesHaveNamedFallbacks(testCase)
shape = radia.python.sheetMetalOptimization( ...
    "route_mesh_update", {[1, 1], [1, 1], [0.01, 0.02]});
verifyEqual(testCase, shape.backend, "python-fallback");
verifyEqual(testCase, shape.module, "radia.sheet_metal_optimization");
topology = radia.python.topologyOptimization( ...
    "solve_lp_update", {[0.5], [-1], [1], 1}, ...
    Keywords=struct("move_limit", 0.1));
verifyEqual(testCase, topology.backend, "python-fallback");
verifyEqual(testCase, topology.module, "radia.topology_optimization");
end

function testHarmonicBalanceHasNamedFallback(testCase)
result = radia.python.harmonicBalance("periodic_phase", {8});
verifyEqual(testCase, result.value, (0:7) * pi / 4, "AbsTol", 2e-15);
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.harmonic_balance");
end

function testXsuiteBridgeHasNamedFallback(testCase)
boundary = radia.python.xsuiteBridge( ...
    "AxisAlignedBox", {[-1, -1, -1], [1, 1, 1]});
x = [0; 0; 2];
y = zeros(3, 1);
z = zeros(3, 1);
result = radia.python.xsuiteBridge( ...
    "first_box_exit_events", {x, y, z, boundary.value});
verifyEqual(testCase, result.backend, "python-fallback");
verifyEqual(testCase, result.module, "radia.xsuite_bridge");
verifyEqual(testCase, result.value.particle_index, 0);
verifyEqual(testCase, result.value.step_index, 2);
verifyEqual(testCase, result.value.position_m, {2, 0, 0});
end

function testExternalPythonUsesChildSafeDllPathAndRestoresMatlabPath(testCase)
if ~ispc
    return
end
radia.setup(Force=true);
activePath = string(getenv("PATH"));
[status, childPath] = radia.internal.runPythonProcess("echo %PATH%");
verifyEqual(testCase, status, 0);
normalized = lower(replace(string(strtrim(childPath)), "/", "\"));
verifyFalse(testCase, contains(normalized, ...
    "\lib\site-packages\ngsolve_openblas"));
verifyFalse(testCase, contains(normalized, ...
    "\lib\site-packages\netgen"));
verifyEqual(testCase, string(getenv("PATH")), activePath);
end
