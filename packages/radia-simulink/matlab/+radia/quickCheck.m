function report = quickCheck()
%QUICKCHECK Verify the CUDA MEX runtime without CPU NGSolve execution.
radia.gpu.requireRuntime();
device = gpuDevice();
report = struct("ready", true, "backend", "mexcuda", ...
    "device", string(device.Name));
end
