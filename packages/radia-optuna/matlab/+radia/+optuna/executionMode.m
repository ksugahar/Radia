function mode = executionMode(requested)
%EXECUTIONMODE Get or select TPE CPU/GPU execution in this MATLAB process.
%   executionMode("cpu") uses OpenMP host kernels in the shared-kernel MEX binary.
%   executionMode("gpu") requires CUDA and never falls back to CPU.
%   The initial mode is "gpu". Clear/reload of optuna_mex resets it.
%   This process-wide setting is not Study state; set it in each worker.
arguments
    requested (1,1) string {mustBeMember(requested,["cpu","gpu"])} = "gpu"
end
if nargin == 0
    mode = string(radia.optuna.internal.NativeKernels.call("api.execution"));
else
    mode = string(radia.optuna.internal.NativeKernels.call("api.execution",char(requested)));
end
end
