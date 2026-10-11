function binary = buildOptunaMex(outputDirectory)
%BUILDOPTUNAMEX Compile CUDA numerical kernels and CPU host control together.
% No NGSolve, TaskManager, Python, Radia solver, or external CUDA override.
arguments
    outputDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')), 'matlab')
end
package = fileparts(mfilename('fullpath'));
configuration = fullfile(matlabroot, 'toolbox', 'parallel', 'gpu', 'extern', ...
    'src', 'mex', 'win64', 'nvcc_msvcpp2022.xml');
assert(isfile(configuration), 'radia:optuna:Compiler', ...
    'This target requires the MATLAB-supported MSVC 2022 mexcuda configuration.');
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
mexcuda('-v', '-R2018a', '-f', configuration, ...
    'NVCCFLAGS=-std=c++20 --fmad=false -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80', ...
    '-DNOMINMAX', '-output', 'optuna_mex', '-outdir', outputDirectory, ...
    fullfile(package, 'native', 'optuna_mex.cpp'), ...
    fullfile(package, 'native', 'optuna_cuda.cu'));
binary = fullfile(outputDirectory, "optuna_mex." + mexext);
assert(isfile(binary), 'radia:optuna:Build', 'mexcuda did not produce optuna_mex.');
end
