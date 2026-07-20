function tests = test_rl_workflow
tests = functiontests(localfunctions);
end

function testIHEnvironmentResetStepAndTermination(testCase)
plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=2, SampleTime_s=0.1);
environment = radia.rl.makeIHEnvironment(plant, ...
    TargetTemperature_K=303.15, MaxSteps=3);

[observation, resetInfo] = environment.reset();
verifySize(testCase, observation, [4, 1]);
verifyEqual(testCase, resetInfo.step, 0);
verifyEqual(testCase, environment.StepCount, 0);

[nextObservation, reward, isDone, stepInfo] = environment.step(100);
verifySize(testCase, nextObservation, [4, 1]);
verifyTrue(testCase, isfinite(reward));
verifyFalse(testCase, isDone);
verifyEqual(testCase, stepInfo.power_W, 100);
verifyGreaterThan(testCase, nextObservation(1), observation(1));

environment.step(100);
[~, ~, isDone] = environment.step(100);
verifyTrue(testCase, isDone);
verifyEqual(testCase, environment.StepCount, 3);
verifyError(testCase, @() environment.step(100), "radia:rl:EpisodeDone");
end

function testGenericEnvironmentAndSeed(testCase)
state = 0;
environment = radia.rl.Environment( ...
    ResetFcn=@resetFcn, StepFcn=@stepFcn, MaxSteps=2);
[observation, info] = environment.reset();
verifyEqual(testCase, observation, 0);
verifyEqual(testCase, info.label, "reset");
[observation, reward, isDone, info] = environment.step(2);
verifyEqual(testCase, observation, 2);
verifyEqual(testCase, reward, -2);
verifyFalse(testCase, isDone);
verifyEqual(testCase, info.label, "step");
environment.seed(11);

    function [observation, info] = resetFcn()
        state = 0;
        observation = state;
        info = struct("label", "reset");
    end

    function [observation, reward, isDone, info] = stepFcn(action)
        state = state + action;
        observation = state;
        reward = -abs(action);
        isDone = false;
        info = struct("label", "step");
    end
end

function testFunctionEnvAdapterContract(testCase)
hasRL = exist("rlFunctionEnv", "file") == 2 || ...
    exist("rlFunctionEnv", "builtin") == 5;
if ~hasRL
    testCase.assumeFail("Reinforcement Learning Toolbox is not installed.");
    return
end

plant = radia.simulink.makeIHPlant( ...
    HeatCapacity_J_per_K=10, ThermalConductance_W_per_K=2, SampleTime_s=0.1);
environment = radia.rl.makeIHEnvironment(plant, MaxSteps=2);
observationInfo = rlNumericSpec([4, 1]);
actionInfo = rlNumericSpec([1, 1], LowerLimit=0, UpperLimit=1e6);
rlEnvironment = environment.toFunctionEnv(observationInfo, actionInfo);
verifyClass(testCase, rlEnvironment, "rl.env.rlFunctionEnv");
end
