function rlEnvironment = toFunctionEnv(environment, observationInfo, actionInfo)
%TOFUNCTIONENV Convert a Radia environment to MATLAB rlFunctionEnv.

arguments
    environment (1,1) radia.rl.Environment
    observationInfo
    actionInfo
end

rlEnvironment = environment.toFunctionEnv(observationInfo, actionInfo);
end
