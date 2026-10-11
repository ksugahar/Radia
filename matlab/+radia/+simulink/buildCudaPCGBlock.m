function blockPath = buildCudaPCGBlock(modelName, options)
%BUILDCUDAPCGBLOCK Add an experimental CUDA-only typed CSR PCG block.
arguments
    modelName (1,1) string
    options.SystemExpression (1,1) string = "cudaSystem"
    options.Tolerance (1,1) double {mustBePositive} = 1e-10
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 10000
    options.SampleTime_s (1,1) double {mustBePositive} = 1e-3
    options.OutputDirectory (1,1) string = "C:\temp\radia_gpu_pcg"
end
if ~bdIsLoaded(modelName), new_system(modelName); end
blockPath = modelName+"/CUDA PCG";
parameters = sprintf('%s,%.17g,%.17g,%.17g,''%s''', ...
    options.SystemExpression,options.Tolerance,options.MaxIterations, ...
    options.SampleTime_s,replace(options.OutputDirectory,"'","''"));
add_block('simulink/User-Defined Functions/Level-2 MATLAB S-Function',char(blockPath), ...
    'FunctionName','radia_gpu_pcg_sfun','Parameters',parameters);
mask = Simulink.Mask.create(blockPath);
mask.Description = 'FP64 CUDA CSR Jacobi-PCG. Matrix upload at Start; RHS upload and solution download each sample. No CPU fallback.';
mask.Display = 'disp(''CUDA PCG'')';
end
