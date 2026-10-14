function accept_payload(job)
% Fresh MATLAB process; no Python adapter or Engine execution in this gate.
addpath(job,fullfile(job,'matlab'));d=gpuDevice();
assert(d.AvailableMemory>max(2*2^30,.2*d.TotalMemory));
native_hcurl(fullfile(job,'model.vol'),job,fullfile(job,'matlab'),fullfile(job,'hcurl.mat'),4);
fmm_physics(fullfile(job,'physics-input.mat'),job,fullfile(job,'matlab'),fullfile(job,'physics.mat'));
fmm_numeric(fullfile(job,'numeric-input.mat'),job,fullfile(job,'matlab'),fullfile(job,'numeric.mat'));
fmm_boundary(fullfile(job,'boundary.vol'),job,fullfile(job,'matlab'),fullfile(job,'boundary.mat'));
assert(radia_cuda_fmm('live_count')==0);
assert(all(structfun(@(x)x==0,radia_cuda_ams('resident_status'))));
pid=feature('getpid');
shell=fullfile(getenv('SystemRoot'),'System32','WindowsPowerShell','v1.0','powershell.exe');
command=sprintf('"%s" -NoProfile -ExecutionPolicy Bypass -File "%s" -ProcessId %d -OutputPath "%s"', ...
    shell,fullfile(job,'modules.ps1'),pid,fullfile(job,'modules.json'));
[status,message]=system(command);assert(status==0,message);
report=struct('accepted',true,'matlab',version,'gpu',d.Name, ...
    'compute_capability',d.ComputeCapability,'driver',d.DriverVersion, ...
    'toolkit',d.ToolkitVersion,'pid',pid,'python_adapter_used',false, ...
    'timing_classification','Same-artifact numerical and DLL acceptance; diagnostic timings only');
file=fopen(fullfile(job,'matlab-result.json'),'w');clean=onCleanup(@()fclose(file));
fprintf(file,'%s',jsonencode(report,PrettyPrint=true));
end
