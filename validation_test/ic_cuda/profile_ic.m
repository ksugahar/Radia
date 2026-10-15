function profile_ic(volPath,payload,icPath,output)
% Explicit diagnostics, separate from end-to-end acceptance measurements.
addpath(payload);addpath(icPath,'-begin');g=gpuDevice;wait(g);
labels=["back","front","left","right","top","bottom"];
timer=tic;
[~,fem,assembly]=radia.gpu.assembleHCurl(volPath,labels,Gamma=0);
source=radia_cuda_ams('nonlinear_create',fem);
sourceCleanup=onCleanup(@()radia_cuda_ams('nonlinear_release',source));
[~,values,~]=radia_cuda_ams('nonlinear_evaluate',source,zeros(size(fem.rhs),'like',fem.rhs));
operator=fem.mass;operator.values=values;operator.free=fem.free;clear sourceCleanup;
wait(g);assemblySeconds=toc(timer);timer=tic;
[h,setup]=radia_cuda_ic('create',operator,1.05);
cleanup=onCleanup(@()radia_cuda_ic('release',h));wait(g);setupSeconds=toc(timer);
seconds=zeros(3,1);linear=cell(3,1);
for k=1:3
    wait(g);timer=tic;
    [x,linear{k}]=radia_cuda_ic('solve',h,fem.rhs,1e-12,20000);
    solution=gather(x);wait(g);seconds(k)=toc(timer);
end
profile=radia_cuda_ic('profile',h,fem.rhs,30);
report=struct('schema','radia.gpu.ic-phase-profile.v1','assembly_seconds',assemblySeconds, ...
    'factor_setup_seconds',setupSeconds,'solve_seconds',seconds,'linear',{linear}, ...
    'setup',setup,'phases',profile,'assembly',assembly,'device',g.Name, ...
    'matlab',version,'binary',which('radia_cuda_ic'), ...
    'timing_scope','Explicit microprofile; excluded from six-route acceptance timing');
save(output,'solution','report','-v7');
fid=fopen(output+".json",'w');assert(fid~=-1);fileCleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
end
