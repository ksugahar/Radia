function build_mexcuda(repoRoot)
%BUILD_OPTUNA_MEXCUDA Independent Optuna gateway using MATLAB's CUDA toolchain.
% MATLAB selects its supported host compiler and bundled CUDA toolkit.
arguments
    repoRoot (1,1) string
end
if ~ispc || computer('arch') ~= "win64"
    error("radia:optuna:Platform", "The distributed gateway requires Windows x64.");
end
if ~license('test','Distrib_Computing_Toolbox')
    error("radia:optuna:CUDARequired", "Parallel Computing Toolbox is required.");
end
sourceRoot = fullfile(repoRoot,"packages","radia-optuna","native");
outputRoot = fullfile(repoRoot,"packages","radia-optuna","matlab");
buildRoot = fullfile(repoRoot,"build-optuna-cuda");
if ~isfolder(buildRoot), mkdir(buildRoot); end
previous = pwd;
cleanup = onCleanup(@()cd(previous));
cd(buildRoot);
fprintf("Optuna mexcuda build: MATLAB %s (%s), %s\n",version,version('-release'),matlabroot);
disp(mex.getCompilerConfigurations('C++','Selected'));
% Strict double precision and no FMA contraction preserve oracle tolerances.
mexcuda('-v','-R2018a','-output','optuna_mex','-outdir',outputRoot, ...
    '-DRADIA_OPTUNA_MEX_ONLY','-DRADIA_OPTUNA_CUDA','-DNOMINMAX', ...
    'NVCC_FLAGS=--std=c++17 --fmad=false --compiler-options=/openmp', ...
    fullfile(sourceRoot,'optuna_mex.cpp'),fullfile(sourceRoot,'optuna_cuda.cu'));
end
