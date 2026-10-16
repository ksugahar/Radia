function result=checkTransactionSnapshot(testArtifact,productionArtifact,outdir,blockUnits)
% Limited RTX transaction/schema acceptance; not B1 or release qualification.
% Requires owned Engine/admission. Exact test+production receipts are mandatory.
% Criteria: all raw accepted bytes+step/epoch/valid equal after normal rejects;
% certified CSR reconstruction <=1e-12 layout / <=1e-8 residual agreement;
% live_count restored after 8 cycles, warm memory loss <=64 MiB (sample only).
arguments
 testArtifact (1,:) char
 productionArtifact (1,:) char
 outdir (1,:) char
 blockUnits cell
end
assert(evalin('base','exist(''radia_transaction_owned_session'',''var'')==1 && isequal(radia_transaction_owned_session,true)'),...
 'Fresh operation-owned launcher session required');
assert(~isfolder(outdir),'New evidence directory required');mkdir(outdir);
result=struct('accepted',false,'b1_closed',false,'release_accepted',false,...
 'timing_classification','transaction/schema diagnostics; not performance acceptance',...
 'fault_kind','host exceptions and simulated poison; no actual asynchronous CUDA fault',...
 'memory_limit_bytes',64*1024^2,'cycles',8);
h=[];artifact='';owner=containers.Map({'h'},{[]});cleanup=onCleanup(@()destroyOwned(owner));
try
 testReceipt=jsondecode(fileread(fullfile(testArtifact,'build.json')));
 productionReceipt=jsondecode(fileread(fullfile(productionArtifact,'build.json')));
 assert(testReceipt.accepted&&productionReceipt.accepted);
 assert(testReceipt.transaction_test_build&&~testReceipt.distribution_eligible);
 assert(~productionReceipt.transaction_test_build&&productionReceipt.distribution_eligible);
 assert(~testReceipt.source_dirty&&~productionReceipt.source_dirty);
 assert(strcmp(testReceipt.source_commit,productionReceipt.source_commit));
 result.source_commit=testReceipt.source_commit;
 result.test_mex_sha256=fileHash(fullfile(testArtifact,'radia_cuda_ih.mexw64'));
 result.production_mex_sha256=fileHash(fullfile(productionArtifact,'radia_cuda_ih.mexw64'));
 assert(strcmp(result.test_mex_sha256,testReceipt.binary_sha256));
 assert(strcmp(result.production_mex_sha256,productionReceipt.binary_sha256));
 switchArtifact(testArtifact);
 [cfg,assembly]=makeTemperatureBHConfig();
 fields={'coil';'port';'incident_potential';'incident_gauge';'body_potential';'body_gauge';'mass_flux';'reaction_potential';'reaction_gauge'};
 assert(isequal(fieldnames(assembly.ordered_unknown_indices),fields),'Reviewed ordered block names required');
 dims=zeros(1,numel(fields));ordered=[];
 for k=1:numel(fields)
  indices=assembly.ordered_unknown_indices.(fields{k});dims(k)=numel(indices);ordered=[ordered,indices]; %#ok<AGROW>
 end
 assert(isequal(ordered,1:size(cfg.eddy,1)),'Complete ordered unknown partition required');
 assert(numel(blockUnits)==numel(dims),'Reviewed physical block units required');
 assert(isequal(blockUnits(:),{'A';'V';'A';'A';'A';'A';'A/m';'V';'V'}),'Reviewed unknown units required');
 cfg.validation_work_snapshot=true;
 cfg.validation_work_coil_port_indices=[assembly.unknown_indices.coil,assembly.unknown_indices.port];
 cfg.validation_work_block_dims=dims;cfg.validation_work_block_units=blockUnits;
 result.block_names=fields;result.block_dims=dims;result.block_units=blockUnits;
 h=ownedCreate(cfg);radia_cuda_ih('update',h,1e4,300);
 before=radia_cuda_ih('transaction_test_capture',h);save(fullfile(outdir,'before.mat'),'before');
 keys={'em','temperature','heat','cells','loss_weights','matrix_values','snapshot_bank'};
 hashes=struct();for k=1:numel(keys),hashes.(keys{k})=byteHash(before.(keys{k}));end
 result.before_accepted_byte_sha256=hashes;
 for stage=1:7
  radia_cuda_ih('transaction_test_failure',h,double(stage));
  if stage==7,failure=expectFailure(@()radia_cuda_ih('reset',h),'Injected');
  else,failure=expectFailure(@()radia_cuda_ih('update',h,1e4,300),'Injected');end
  after=radia_cuda_ih('transaction_test_capture',h);same(before,after);save(fullfile(outdir,sprintf('stage-%d.mat',stage)),'before','after','stage','failure');
 end
 failure=expectFailure(@()radia_cuda_ih('update',h,1e4,1e12),'temperature');
 after=radia_cuda_ih('transaction_test_capture',h);same(before,after);
 save(fullfile(outdir,'thermal-rejection.mat'),'before','after','failure');
 radia_cuda_ih('transaction_test_constitutive_limit',h,2);
 failure=expectFailure(@()radia_cuda_ih('update',h,1e4,300),'constitutive');
 after=radia_cuda_ih('transaction_test_capture',h);same(before,after);
 save(fullfile(outdir,'constitutive-rejection.mat'),'before','after','failure');
 radia_cuda_ih('transaction_test_constitutive_limit',h,300);
 radia_cuda_ih('reset',h);after=radia_cuda_ih('transaction_test_capture',h);
 assert(after.steps==0&&~after.snapshot_valid&&after.snapshot_generation==before.snapshot_generation+1);
 expectSnapshotFailure(h);
 stale=h;radia_cuda_ih('destroy',h);h=[];owner('h')=[];expectSnapshotFailure(stale);
 % Explicit non-symmetric, non-Hermitian perturbation tests serialization only.
 cfg.eddy(2,3)=cfg.eddy(2,3)+(1e-6+2e-6i);
 h=ownedCreate(cfg);radia_cuda_ih('update',h,1e4,300);
 snapshot=radia_cuda_ih('accepted_work_snapshot',h);
 n=double(snapshot.n_rows);a=sparse(n,n);z=complex(snapshot.A_certified_re,snapshot.A_certified_im);
 for row=1:n
  range=double(snapshot.csr_rowptr(row))+1:double(snapshot.csr_rowptr(row+1));
  a(row,double(snapshot.csr_col(range))+1)=z(range);
 end
 zs=complex(snapshot.Zs_re,snapshot.Zs_im);
 expected=cfg.eddy+reshape(cfg.operator_impedance_map*zs,n,n).';
 assert(norm(a-a.','fro')>1e-12&&norm(a-a','fro')>1e-12);
 layoutError=norm(a-expected,'fro')/max(norm(expected,'fro'),realmin);
 assert(layoutError<=1e-12);
 x=complex(snapshot.x_re,snapshot.x_im);b=complex(snapshot.rhs_re,snapshot.rhs_im);
 residual=norm(a*x-b)/max(norm(b),1*(norm(b)==0));
 assert(residual<=1e-6&&abs(residual-snapshot.certified_true_residual_native)<=1e-8);
 result.layout_relative=layoutError;result.recomputed_certified_residual=residual;
 save(fullfile(outdir,'snapshot.mat'),'snapshot','cfg','assembly');
 [~,~,onInfo]=radia_cuda_ih('outputs',h);result.on_info=onInfo;
 radia_cuda_ih('transaction_test_poison',h);
 expectSnapshotFailure(h);expectFailure(@()radia_cuda_ih('update',h,1e4,300),'poisoned');
 for command={'outputs','reset','electromagnetic_state','element_impedance','diagnostic_operator','transaction_test_capture'}
  expectFailure(@()poisonCommand(command{1},h),'poisoned');
 end
 radia_cuda_ih('destroy',h);h=[];owner('h')=[];assert(radia_cuda_ih('live_count')==0);
 memoryBefore=gpuDevice().AvailableMemory;memorySamples=zeros(1,8);
 for k=1:8
  h=ownedCreate(cfg);radia_cuda_ih('update',h,1e4,300);
  radia_cuda_ih('destroy',h);h=[];owner('h')=[];assert(radia_cuda_ih('live_count')==0);
  assert(~mislocked('radia_cuda_ih'));memorySamples(k)=gpuDevice().AvailableMemory;
 end
 memoryAfter=gpuDevice().AvailableMemory;assert(memoryBefore-memoryAfter<=64*1024^2);
 result.warm_available_memory_before=memoryBefore;result.warm_available_memory_after=memoryAfter;
 result.cycle_available_memory=memorySamples;
 result.memory_strictly_decreasing=all(diff(memorySamples)<0);
 assert(~result.memory_strictly_decreasing,'Strict per-cycle memory loss requires investigation');
 switchArtifact(productionArtifact);cfg.validation_work_snapshot=false;
 h=ownedCreate(cfg);radia_cuda_ih('update',h,1e4,300);
 [~,~,offInfo]=radia_cuda_ih('outputs',h);result.off_info=offInfo;
 assert(offInfo.snapshot_bytes==0&&offInfo.atomicity_fix_bytes>0);
 assert(onInfo.snapshot_bytes>0&&onInfo.atomicity_fix_bytes==offInfo.atomicity_fix_bytes);
 assert(isfinite(offInfo.atomicity_stage_sync_s)&&offInfo.atomicity_stage_sync_s>=0);
 expectFailure(@()radia_cuda_ih('transaction_test_poison',h),'Unknown');
 expectSnapshotFailure(h);radia_cuda_ih('destroy',h);h=[];owner('h')=[];
 result.matlab=version;device=gpuDevice();result.gpu=device.Name;
 result.driver=device.GraphicsDriverVersion;result.cuda=device.ToolkitVersion;
 result.stage_rejections=7;result.thermal_rejected=true;result.esim_rejected=true;
 result.accepted=true;
catch err
 result.error=getReport(err,'extended','hyperlinks','off');
end
finish();
 function same(a,b)
  for q=1:numel(keys),assert(isequal(a.(keys{q}),b.(keys{q})),keys{q});end
  for key={'steps','snapshot_step','snapshot_valid','snapshot_generation','handle_generation','device_uuid'}
   assert(isequal(a.(key{1}),b.(key{1})),key{1});
  end
 end
 function switchArtifact(path)
  assert(isempty(h));clear radia_cuda_ih;
  if ~isempty(artifact),rmpath(artifact);end
  addpath(path);artifact=path;
  assert(strcmpi(char(java.io.File(which('radia_cuda_ih')).getCanonicalPath()),char(java.io.File(fullfile(path,'radia_cuda_ih.mexw64')).getCanonicalPath())));
 end
 function handle=ownedCreate(config)
  handle=radia_cuda_ih('create',config);owner('h')=handle;
 end
 function finish()
  destroyOwned(owner);
  f=fopen(fullfile(outdir,'result.json'),'w');assert(f>=0);writer=onCleanup(@()fclose(f));fprintf(f,'%s\n',jsonencode(result));
 end
end
function poisonCommand(command,h)
if strcmp(command,'outputs')
 [~,~,~]=radia_cuda_ih(command,h);
elseif strcmp(command,'reset')
 radia_cuda_ih(command,h);
else
 unused=radia_cuda_ih(command,h); %#ok<NASGU>
end
end
function message=expectFailure(action,word)
message='returned normally';failed=false;
try
 action();
catch err
 message=['exception identifier: ',err.identifier,'; message: ',err.message];
 failed=contains(lower(err.message),lower(word));
end
assert(failed,['Expected failure: ',word,'; actual message: ',message]);
end
function expectSnapshotFailure(h)
failed=false;try,s=radia_cuda_ih('accepted_work_snapshot',h);catch,failed=true;end %#ok<NASGU>
assert(failed);
end
function digest=fileHash(path)
f=fopen(path,'rb');assert(f>=0);closer=onCleanup(@()fclose(f));digest=byteHash(fread(f,Inf,'*uint8'));
end
function digest=byteHash(bytes)
md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(uint8(bytes(:)),'int8'));digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function destroyOwned(owner)
handle=owner('h');
owner('h')=[];
assert(isempty(owner('h')));
if ~isempty(handle)
 try
  radia_cuda_ih('destroy',handle);
 catch err
  if ~contains(err.message,'stale'),warning('radia:ih:TestCleanup','%s',err.message);end
 end
end
end
