function tests = test_rl_workflow
tests = functiontests(localfunctions);
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

environment = makeCounterEnvironment();
observationInfo = rlNumericSpec([1, 1]);
actionInfo = rlNumericSpec([1, 1], LowerLimit=0, UpperLimit=1e6);
rlEnvironment = environment.toFunctionEnv(observationInfo, actionInfo);
verifyClass(testCase, rlEnvironment, "rl.env.rlFunctionEnv");
end

function environment = makeCounterEnvironment()
state = 0;
environment = radia.rl.Environment( ...
    ResetFcn=@resetCounter,StepFcn=@stepCounter,MaxSteps=2);

    function [observation,info] = resetCounter()
        state = 0;
        observation = state;
        info = struct("step",0);
    end

    function [observation,reward,isDone,info] = stepCounter(action)
        state = state + action;
        observation = state;
        reward = -abs(action);
        isDone = false;
        info = struct("step",state);
    end
end
