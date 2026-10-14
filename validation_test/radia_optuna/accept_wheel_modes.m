function receipt=accept_wheel_modes(installedMatlabPath, testsDirectory, outputPath)
%ACCEPT_WHEEL_MODES Verify the installed distribution on both explicit routes.
arguments
    installedMatlabPath (1,1) string
    testsDirectory (1,1) string
    outputPath (1,1) string
end
addpath(installedMatlabPath,'-begin');
addpath(testsDirectory,'-begin');
assert(startsWith(string(which('radia.optuna.Study')),installedMatlabPath+filesep));
assert(startsWith(string(which('optuna_mex')),installedMatlabPath+filesep));
suite=testsuite(fullfile(testsDirectory,'test_optuna_upstream_oracle.m'));
selected=suite(~cellfun(@isempty,regexp({suite.Name},'TPE|Random|Sobol|Mixed|MultiObjective|SingleDistribution')));
history=testsuite(fullfile(testsDirectory,'test_optuna_table.m'));
history=history(contains({history.Name},'NamedHistory'));
extra=testsuite(fullfile(testsDirectory,'test_optuna_cuda.m'));
radia.optuna.executionMode('cpu');
cpuBefore=optuna_mex('api.info');
cpuTests=run([selected(:);history(:)]);
assert(all([cpuTests.Passed]) && numel(cpuTests)==16);
assert(startsWith(string(which('radia.optuna.Study')),installedMatlabPath+filesep));
cpuStandalone=test_standalone_simulink(installedMatlabPath);
cpuAfter=optuna_mex('api.info');
assert(cpuAfter.cuda_device==-1 && cpuAfter.host_to_device_bytes==0 && cpuAfter.device_to_host_bytes==0);
assert(cpuAfter.density_launches==0 && cpuAfter.cpu_density_calls>cpuBefore.cpu_density_calls);
radia.optuna.executionMode('gpu');
gpuBefore=optuna_mex('api.info');
gpuTests=run([selected(:);extra(:)]);
assert(all([gpuTests.Passed]) && numel(gpuTests)==21);
assert(startsWith(string(which('radia.optuna.Study')),installedMatlabPath+filesep));
gpuOnlyBefore=optuna_mex('api.info');
gpuStandalone=test_standalone_simulink(installedMatlabPath);
gpuAfter=optuna_mex('api.info');
assert(gpuAfter.model_build_launches>gpuBefore.model_build_launches && gpuAfter.model_radix_sorts>gpuBefore.model_radix_sorts);
assert(gpuAfter.cpu_preparation_calls==gpuOnlyBefore.cpu_preparation_calls);
manifest=suite(strcmp({suite.Name},'test_optuna_upstream_oracle/testOraclePolicyManifestIsComplete'));
assert(numel(manifest)==1); manifestResult=run(manifest); assert(all([manifestResult.Passed]));
d=gpuDevice;
receipt=struct('schema','radia.validation.optuna-installed-cpu-gpu.v1','ok',true, ...
    'matlab',version,'gpu',d.Name,'compute_capability',d.ComputeCapability, ...
    'toolkit',d.ToolkitVersion,'driver',d.GraphicsDriverVersion, ...
    'oracle','Optuna 5.0.0','cpu_passed',numel(cpuTests),'gpu_passed',numel(gpuTests), ...
    'manifest_passed',1,'cpu_before',cpuBefore,'cpu_after',cpuAfter, ...
    'gpu_before',gpuBefore,'gpu_after',gpuAfter, ...
    'cpu_simulink',cpuStandalone,'gpu_simulink',gpuStandalone);
fid=fopen(outputPath,'w');assert(fid~=-1);closeFile=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));
radia.optuna.executionMode('cpu');clear optuna_mex
end
