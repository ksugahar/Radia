function tests = test_circuit_field_simulink_mex
tests = functiontests(localfunctions);
end

function testSingleBranchExactZOH(testCase)
model = radia.simulink.makeCircuitFieldStateSpace(4, 2, 2, ...
    SampleTime_s=0.1);

verifyEqual(testCase, model.schema, "radia.circuit-field.state-space.v1");
verifyEqual(testCase, model.backend, "native-mex-sfunction");
verifyEqual(testCase, model.inductance_H, 1, "AbsTol", 1e-12);
verifyEqual(testCase, model.A, -2, "AbsTol", 1e-12);
verifyEqual(testCase, model.B, 1, "AbsTol", 1e-12);
verifyEqual(testCase, model.Ad, exp(-0.2), "AbsTol", 1e-12);
verifyEqual(testCase, model.Bd, (1 - exp(-0.2)) / 2, "AbsTol", 1e-12);
verifyFalse(testCase, model.python_per_step);
verifyFalse(testCase, model.field_factorization_per_step);
end

function testTwoBranchCommonVoltageDimensions(testCase)
K = [4, 1; 1, 3];
source = eye(2);
model = radia.simulink.makeCircuitFieldStateSpace(K, source, [1; 2], ...
    SampleTime_s=1e-4, VoltageInputMode="common");

verifySize(testCase, model.Ad, [2, 2]);
verifySize(testCase, model.Bd, [2, 1]);
verifySize(testCase, model.Cd, [4, 2]);
verifySize(testCase, model.Dd, [4, 1]);
verifyLessThan(testCase, max(real(eig(model.A))), 0);
verifyEqual(testCase, model.mex_s_function, "radia_state_space_mex_sfunction");

perBranch = radia.simulink.makeCircuitFieldStateSpace(K, source, [1; 2], ...
    SampleTime_s=1e-4, VoltageInputMode="per_branch");
verifySize(testCase, perBranch.Bd, [2, 2]);
end

function testBuilderUsesExistingNativeSFunction(testCase)
assumeTrue(testCase, exist("new_system", "file") == 2 || ...
    exist("new_system", "builtin") == 5);
name = "radia_circuit_field_mex_test";
cleanup = onCleanup(@() closeIfLoaded(name));
model = radia.simulink.makeCircuitFieldStateSpace(4, 2, 2);
radia.simulink.buildCircuitFieldStateSpaceModel(name, model, ...
    Save=false, Open=false);

verifyEqual(testCase, string(get_param(name + "/Circuit Field ROM", ...
    "FunctionName")), "radia_state_space_mex_sfunction");
contract = getVariable(get_param(name, "ModelWorkspace"), ...
    "radia_circuit_field_contract");
verifyEqual(testCase, string(contract.backend), "native-mex-sfunction");
verifyFalse(testCase, contract.python_per_step);
clear cleanup
end

function closeIfLoaded(name)
if bdIsLoaded(name)
    close_system(name, 0);
end
end
