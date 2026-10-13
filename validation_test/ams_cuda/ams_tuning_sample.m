function ams_tuning_sample(job,binaryPath,output,matlabPath)
% Numerical timing excludes fixture loading and output serialization.
if nargin<4;matlabPath=fullfile(job,'matlab');end
addpath(matlabPath);addpath(binaryPath,'-begin');g=gpuDevice;wait(g);
initial=radia_cuda_ams('resident_status');labels=["back","front","left","right","top","bottom"];
before=g.AvailableMemory;timer=tic;
[nonlinear,nonlinearInfo]=radia.gpu.solveHCurl(fullfile(job,'model.vol'),labels, ...
    Gamma=4,CurrentDensity=[20 0 0],SetupThreads=1,LinearTolerance=1e-12,ReturnOnGPU=false);
wait(g);nonlinearSeconds=toc(timer);timer=tic;
[linear,linearInfo]=radia.gpu.solveHCurl(fullfile(job,'model.vol'),labels, ...
    Gamma=0,CurrentDensity=[1 0 0],SetupThreads=1,LinearTolerance=1e-12,ReturnOnGPU=false);
wait(g);linearSeconds=toc(timer);timer=tic;
system=radia.gpu.HCurlSystem(fullfile(job,'model.vol'),labels,SetupThreads=1);
clean=onCleanup(@()delete(system));wait(g);preparedSeconds=toc(timer);
retainedBefore=g.AvailableMemory;scales=[1 2 .5 1 0 1];reuse=zeros(numel(linear),numel(scales));
reuseSeconds=zeros(numel(scales),1);reuseInfo=cell(numel(scales),1);
for k=1:numel(scales)
    wait(g);timer=tic;
    [reuse(:,k),reuseInfo{k}]=system.solve(CurrentScale=scales(k),Tolerance=1e-12,ReturnOnGPU=false);
    wait(g);reuseSeconds(k)=toc(timer);
    assert(norm(reuse(:,k)-scales(k)*linear)<=1e-9*max(1,norm(scales(k)*linear)));
end
retainedAfter=g.AvailableMemory;setupInfo=system.SetupInfo;clear clean;
try;system.solve();error('radia:test:missingError','Released system accepted solve');
catch ex;assert(strcmp(ex.identifier,'radia:cuda:released')||contains(ex.identifier,'InvalidHandle'),ex.identifier);end
assert(isequal(radia_cuda_ams('resident_status'),initial),'radia:test:leak','AMS state leaked');
% Read the loaded Windows module, not merely MATLAB's search-path resolution.
code=sprintf(['$ProgressPreference=''SilentlyContinue'';$ErrorActionPreference=''Stop'';' ...
    '@(Get-Process -Id %d -Module | Where-Object {$_.ModuleName -eq ''radia_cuda_ams.mexw64''} | ' ...
    'ForEach-Object {[PSCustomObject]@{path=$_.FileName;sha256=(Get-FileHash -LiteralPath $_.FileName -Algorithm SHA256).Hash.ToLower()}}) | ConvertTo-Json -Compress'],feature('getpid'));
encoded=matlab.net.base64encode(unicode2native(code,'UTF-16LE'));
[status,loaded]=system(['powershell -NoProfile -EncodedCommand ' encoded]);
assert(status==0,'radia:test:module','Could not inspect loaded native module');
nativeModule=jsondecode(loaded);assert(isscalar(nativeModule),'radia:test:module','Mixed AMS native modules loaded');
report=struct('nonlinear_seconds',nonlinearSeconds,'linear_fresh_seconds',linearSeconds, ...
    'prepared_setup_seconds',preparedSeconds,'reused_seconds',reuseSeconds,'nonlinear',nonlinearInfo, ...
    'linear',linearInfo,'reuse', {reuseInfo},'setup',setupInfo,'scales',scales, ...
    'device_free_before',before,'device_free_retained_before',retainedBefore, ...
    'device_free_retained_after',retainedAfter,'device_free_after',g.AvailableMemory, ...
    'loaded_binary',which('radia_cuda_ams'),'native_module',nativeModule,'matlab',version,'device',g.Name,'registry_restored',true);
save(output,'nonlinear','linear','reuse','report');
clear system;clear radia_cuda_ams;rmpath(binaryPath);
end
