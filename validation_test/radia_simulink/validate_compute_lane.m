function validate_compute_lane(payload,optunaPayload,output)
% Compute-only use of exact installed wheels. No compiler or Python solver.
restoredefaultpath;
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'validation_test','radia_simulink'));
if ~isfolder(output), mkdir(output); end
check_installed(payload,fullfile(output,'installed'),false);
addpath(payload,'-begin');
g=gpuDevice; assert(contains(g.Name,'A100'),'A100 compute capability required');
assert(string(g.ComputeCapability)=="8.0",'Unexpected compute capability');
assert(g.TotalMemory>=39e9,'Expected A100 40 GB or greater');
records=cell(3,1);
for k=1:3
    n=[32768 131072 524288]; n=n(k);
    % Strictly SPD tridiagonal, analytical x=1. CSR uses zero-based int32.
    rowptr=int32([0;2;2+3*(1:n-2)';3*n-2]);
    columns=zeros(3*n-2,1,'int32'); values=zeros(3*n-2,1);
    cursor=1;
    for row=1:n
        cols=max(1,row-1):min(n,row+1); count=numel(cols);
        columns(cursor:cursor+count-1)=int32(cols-1);
        values(cursor:cursor+count-1)=-1; values(cursor+find(cols==row)-1)=4;
        cursor=cursor+count;
    end
    rhs=2*ones(n,1);rhs([1 end])=3;
    times=zeros(4,1); cold=0;
    for repeat=1:4
        wait(g); timer=tic;
        [x,info]=radia.gpu.csrPCG(gpuArray(rowptr),gpuArray(columns), ...
            gpuArray(values),gpuArray(rhs),gpuArray(.25*ones(n,1)),1e-10,1000);
        x=gather(x);wait(g);times(repeat)=toc(timer);
        if repeat==1,cold=times(repeat);end
        applied=4*x-[0;x(1:end-1)]-[x(2:end);0];
        residual=norm(applied-rhs)/norm(rhs);
        assert(residual<=1e-6 && info.true_relative_residual<=1e-6);
        assert(norm(x-1)/sqrt(n)<1e-8);
    end
    records{k}=struct('rows',n,'cold_seconds',cold, ...
        'warmed_seconds',times(2:end),'warmed_median_seconds',median(times(2:end)), ...
        'transfer_inclusive',true,'true_relative_residual',residual,'solver',info);
end
receipt=struct('schema','radia.gpu.compute-validation.v1','passed',true, ...
    'device',g.Name,'compute_capability',g.ComputeCapability,'memory_bytes',g.TotalMemory, ...
    'driver',g.GraphicsDriverVersion,'toolkit',g.ToolkitVersion,'matlab',version, ...
    'cases',{records},'source_commit',getenv('GITHUB_SHA'),'backend','gpu');
fid=fopen(fullfile(output,'csr.json'),'w');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));fclose(fid);
addpath(fullfile(root,'validation_test','radia_optuna'));
validate_cuda(root,optunaPayload,fullfile(output,'optuna.json'),131072);
end
