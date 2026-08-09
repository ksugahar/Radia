function tests = test_team28_coilbuilder_dynamic
%TEST_TEAM28COILBUILDERDYNAMIC Unit and artifact gates for the moving model.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = string(fileparts(mfilename("fullpath")));
repoRoot = string(fileparts(fileparts(here)));
addpath(fullfile(repoRoot, "matlab"));
testCase.TestData.RepoRoot = repoRoot;
end

function testValidatedCoilBuilderLUT(testCase)
lut = radia.simulink.makeTeam28CoilBuilderLUT();
verifyEqual(testCase, lut.schema, "radia.team28.cln_lut.v1");
verifyTrue(testCase, lut.source_validation_passed);
verifyEqual(testCase, lut.family_state_order, 3);
verifyEqual(testCase, lut.family_snapshot_count, 25);
verifyEqual(testCase, lut.frequency_Hz, 50.0, "AbsTol", 0);
verifyEqual(testCase, lut.physical_force_factor, 1.0, "AbsTol", 0);
verifyLessThan(testCase, ...
    lut.coilbuilder_vector_potential_max_relative_l2, 1.0e-3);
verifyLessThan(testCase, ...
    lut.coilbuilder_flux_density_max_relative_l2, 5.0e-3);
verifyLessThan(testCase, lut.lift_curve_global_relative_error, 0.02);
verifyEmpty(testCase, regexp(lut.source_file, ...
    '^[A-Za-z]:[\\/]|^[/\\]{2}', "once"));
verifyTrue(testCase, startsWith(replace(lut.source_file, "\", "/"), ...
    "validation_test/maglev/"));
verifyTrue(testCase, all(diff(lut.height_offset_m) > 0));
[force_N, lift_N] = radia.simulink.evaluateTeam28CLNForce(lut, 0.0, 20.0);
verifyEqual(testCase, force_N, -lift_N, "AbsTol", 1e-12);
verifyEqual(testCase, lift_N, 1.1019289804974595, "AbsTol", 1e-12);
end

function testMechanicalContractHasStableEquilibrium(testCase)
lut = radia.simulink.makeTeam28CoilBuilderLUT();
plant = radia.simulink.makeTeam28CoilBuilderPlant(lut);
verifyEqual(testCase, plant.schema, ...
    "radia.team28.coilbuilder_mechanical.v1");
verifyGreaterThan(testCase, plant.disk_mass_kg, 0);
verifyGreaterThan(testCase, ...
    plant.linearized_restoring_stiffness_N_per_m, 0);
verifyGreaterThan(testCase, plant.viscous_damping_Ns_per_m, 0);
verifyEqual(testCase, plant.damping_ratio, 0.35, "AbsTol", 0);
verifyEqual(testCase, plant.excitation_frequency_hz, 50.0, "AbsTol", 0);
verifyEqual(testCase, plant.electromagnetic_model_class, ...
    "fixed_frequency_cycle_averaged_force_height_lut");
verifyEqual(testCase, plant.height_coupling, "quasi_steady_interpolation");
verifyFalse(testCase, plant.electromagnetic_state_transient_included);
verifyFalse(testCase, plant.motional_emf_included);
verifyFalse(testCase, plant.damping_identified_from_measurement);
verifyEqual(testCase, plant.force_family_snapshot_count, 25);
verifyEqual(testCase, plant.eddy_state_order, 3);
verifyLessThan(testCase, abs(plant.equilibrium_absolute_bottom_m - ...
    0.0110555366063325), 1e-12);
end

function testSavedModelUsesLevel2BlockAndUpdates(testCase)
hasSimulink = exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5;
if ~hasSimulink
    testCase.assumeFail("Simulink is unavailable on this MATLAB runtime.");
end
modelFile = fullfile(testCase.TestData.RepoRoot, "matlab", ...
    "radia_team28_coilbuilder_dynamic.slx");
verifyTrue(testCase, isfile(modelFile));
load_system(modelFile);
cleanup = onCleanup(@() closeIfLoaded("radia_team28_coilbuilder_dynamic"));
block = "radia_team28_coilbuilder_dynamic/" + ...
    "CoilBuilder HCurl Eddy-Bubble Lift";
verifyEqual(testCase, string(get_param(block, "BlockType")), ...
    "M-S-Function");
verifyEqual(testCase, string(get_param(block, "FunctionName")), ...
    "radia_team28_cln_lut_sfunction");
verifyEqual(testCase, string(get_param(block, "Parameters")), ...
    "radia_team28_coilbuilder_lut");
set_param("radia_team28_coilbuilder_dynamic", ...
    "SimulationCommand", "update");
clear cleanup
closeIfLoaded("radia_team28_coilbuilder_dynamic");
end

function testDynamicResultArtifact(testCase)
resultFile = fullfile(testCase.TestData.RepoRoot, "validation_test", ...
    "maglev", "team28_coilbuilder_dynamic_simulink_results.json");
verifyTrue(testCase, isfile(resultFile));
result = jsondecode(fileread(resultFile));
names = string(fieldnames(result));
verifyEqual(testCase, names(1), "radia_version");
verifyTrue(testCase, result.pass);
verifyTrue(testCase, result.checks.validation_passed);
verifyLessThan(testCase, result.errors.terminal_height_abs_m, 0.1e-3);
verifyLessThan(testCase, ...
    result.errors.terminal_force_balance_abs_N, 0.02);
end

function closeIfLoaded(modelName)
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end
