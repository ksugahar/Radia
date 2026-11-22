function native_hcurl(volPath,binaryPath,matlabPath,outputPath,gamma)
% End-to-end .vol -> native FE -> CUDA assembly -> native AMS -> GPU solve.
if nargin<5;gamma=0;end
current=1;if gamma>0;current=20;end
addpath(binaryPath);addpath(matlabPath);device=gpuDevice();wait(device);
before=double(device.AvailableMemory);initial=radia_cuda_ams('resident_status');total=tic;
labels=["back","front","left","right","top","bottom"];
[fem,assembly,model]=radia_cuda_ams_validation('assemble_hcurl',char(volPath), ...
    cellstr(labels),'.*',1,1,[current 0 0],gamma,true);
[solution,solverReport]=radia.gpu.solveHCurl(volPath, ...
    ["back","front","left","right","top","bottom"],Gamma=gamma,CurrentDensity=[current 0 0],ReturnOnGPU=false);
rhs=gather(fem.rhs);wait(device);elapsed=toc(total);
mass_values=gather(fem.mass.values);
maps=struct('dofs',gather(fem.dofs),'curls',gather(fem.curls),'volume',gather(fem.volume));
assert(isequal(radia_cuda_ams('resident_status'),initial),'radia:cuda:leak','AMS contexts leaked');
after=double(device.AvailableMemory);
warmed=zeros(3,1);
for k=1:3
    wait(device);timer=tic;
    [repeat,~]=radia.gpu.solveHCurl(volPath, ...
        ["back","front","left","right","top","bottom"],Gamma=gamma,CurrentDensity=[current 0 0],ReturnOnGPU=false);
    wait(device);warmed(k)=toc(timer);
    assert(norm(repeat-solution)<=1e-10*max(norm(solution),realmin),'radia:cuda:repeat','Repeated solve changed solution');
end
assert(isequal(radia_cuda_ams('resident_status'),initial),'radia:cuda:leak','Repeated AMS contexts leaked');
report=struct('total_seconds',elapsed,'assembly',assembly,'solver',solverReport, ...
    'gpu_name',device.Name,'available_before',before,'available_after',after, ...
    'retained_device_memory_sample_bytes',max(0,before-after),'matlab_version',version, ...
    'timing_classification','100 build-host diagnostic; not formal performance acceptance', ...
    'warmed_transfer_inclusive_seconds',warmed,'warmed_transfer_inclusive_median',median(warmed), ...
    'cold_production_seconds',solverReport.total_seconds, ...
    'total_seconds_definition','Validation includes a separate assembly export in addition to the production solve');
linearAMS=radia.gpu.HCurlSystem(volPath,labels,CurrentDensity=[current 0 0]);
linearIC=radia.gpu.HCurlICSystem(volPath,labels,CurrentDensity=[current 0 0]);
[ams1,amsReport1]=linearAMS.solve(ReturnOnGPU=false);
[ams2,amsReport2]=linearAMS.solve(CurrentScale=.5,ReturnOnGPU=false);
[ams0,amsReport0]=linearAMS.solve(CurrentScale=0,ReturnOnGPU=false);
[ic1,icReport1]=linearIC.solve(ReturnOnGPU=false);
[ic2,icReport2]=linearIC.solve(CurrentScale=.5,ReturnOnGPU=false);
assert(all(ams0==0)&&amsReport0.linear.true_relative_residual==0);
assert(linearAMS.SetupInfo.ams_setup.hierarchy_built==1&& ...
    amsReport1.linear.hierarchy_built==0&&amsReport2.linear.hierarchy_built==0&& ...
    amsReport1.linear.hierarchy_refreshed==0&&amsReport2.linear.hierarchy_refreshed==0);
assert(icReport1.factor_reused&&icReport2.factor_reused);
icSetup=linearIC.SetupInfo.ic_setup;
assert(icSetup.symbolic_download_bytes==0&&icSetup.symbolic_upload_bytes==0);
for setup={linearAMS.SetupInfo,linearIC.SetupInfo}
    assert(setup{1}.assembly.matrix_h2d_bytes==0&&setup{1}.assembly.matrix_d2h_bytes==0);
end
linearChecks=struct('ams_setup',linearAMS.SetupInfo,'ic_setup',linearIC.SetupInfo, ...
    'ams_first',amsReport1,'ams_reuse',amsReport2,'ic_first',icReport1,'ic_reuse',icReport2);
delete(linearAMS);delete(linearIC);
assert(isequal(radia_cuda_ams('resident_status'),initial));
for command=["setup","nonlinear_raw","create"]
    rejected=false;
    try;radia_cuda_ams(char(command));catch e;rejected=contains(e.message,'Unknown AMS command');end
    assert(rejected);
end
% Bounded cases only. A larger validator must select checked HDF5 beforehand.
save(outputPath,'solution','rhs','model','mass_values','maps','report','linearChecks','ams1','ams2','ic1','ic2');
end
