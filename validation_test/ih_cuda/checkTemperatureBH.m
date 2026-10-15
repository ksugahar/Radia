function result=checkTemperatureBH(artifact,outdir,matlabRoot)
% Bounded complete reaction + surface thermal research acceptance.
if ~isfolder(outdir),mkdir(outdir);end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin<3,matlabRoot=fullfile(root,'packages','radia-gpu','matlab');end
addpath(matlabRoot,fileparts(mfilename('fullpath')));addpath(artifact);
result=struct('accepted',false,'research_only',true,'experimental_validation',false,'toolkit_accepted',false);
try
 [cfg,assembly]=makeTemperatureBHConfig();save(fullfile(outdir,'input.mat'),'cfg','assembly');
 assert(radia_cuda_ih('live_count')==0);
 for invalidCase=1:4
  bad=cfg;
  if invalidCase==1,bad.bh_values_T(1,1,2)=NaN;end
  if invalidCase==2,bad.bh_knots(2)=0;end
  if invalidCase==3,bad.bh_values_T(1,1,end)=0;end
  if invalidCase==4,bad=rmfield(bad,'bh_values_T');end
  failed=false;
  try
   badHandle=radia_cuda_ih('create',bad);radia_cuda_ih('destroy',badHandle);
  catch
   failed=true;
  end
  assert(failed&&radia_cuda_ih('live_count')==0);
 end
 bad=cfg;bad.initial_temperature_K=700*ones(4,1);
 badHandle=radia_cuda_ih('create',bad);badCleaner=onCleanup(@()radia_cuda_ih('destroy',badHandle));
 failed=false;
 try
  radia_cuda_ih('update',badHandle,1e4,300);
 catch
  failed=true;
 end
 [badT,badQ,badInfo]=radia_cuda_ih('outputs',badHandle);
 assert(failed&&all(badT==700)&&all(badQ==0)&&badInfo.steps==0);clear badCleaner;
 badHandle=radia_cuda_ih('create',cfg);badCleaner=onCleanup(@()radia_cuda_ih('destroy',badHandle));
 failed=false;
 try
  radia_cuda_ih('update',badHandle,1e4,1e8);
 catch
  failed=true;
 end
 [badT,badQ,badInfo]=radia_cuda_ih('outputs',badHandle);
 assert(failed&&isequal(badT,cfg.initial_temperature_K)&&all(badQ==0)&&badInfo.steps==0&&all(radia_cuda_ih('element_impedance',badHandle)==0));clear badCleaner;
 h=radia_cuda_ih('create',cfg);cleaner=onCleanup(@()radia_cuda_ih('destroy',h));
 t=zeros(6,4);q=t;z=complex(t);balance=zeros(5,1);
 for k=0:5
  [temp,heat,info]=radia_cuda_ih('outputs',h);t(k+1,:)=temp;q(k+1,:)=heat;z(k+1,:)=radia_cuda_ih('element_impedance',h);
  if k>0
   previous=t(k,:)';
   r=assembly.thermal_capacity*(temp-previous)/cfg.sample_time_s+cfg.stiffness*temp+cfg.convection*(temp-300)-cfg.heat_projection*heat;
   balance(k)=norm(r)/max(norm(cfg.heat_projection*heat)+norm(cfg.convection*(temp-300)),1e-20);
   assert(all(isfinite([info.eddy_true_residual,info.constitutive_true_residual,info.thermal_true_residual,balance(k)])));
   assert(info.eddy_true_residual<=1e-6 && info.constitutive_true_residual<=1e-8 && info.thermal_true_residual<=1e-6 && balance(k)<=1e-6);
  end
  if k<5,radia_cuda_ih('update',h,1e4,300);end
 end
 em=radia_cuda_ih('electromagnetic_state',h);save(fullfile(outdir,'accepted.mat'),'t','q','z','em');rejected=false;
 try
  radia_cuda_ih('update',h,realmax,300);
 catch
  rejected=true;
 end
 [afterT,afterQ,afterInfo]=radia_cuda_ih('outputs',h);
 assert(rejected&&isequal(afterT,t(end,:)')&&isequal(afterQ,q(end,:)')&&afterInfo.steps==5&&isequal(radia_cuda_ih('element_impedance',h),z(end,:).')&&isequal(em,radia_cuda_ih('electromagnetic_state',h)));
 radia_cuda_ih('reset',h);[resetT,resetQ,resetInfo]=radia_cuda_ih('outputs',h);
 assert(isequal(resetT,cfg.initial_temperature_K)&&all(resetQ==0)&&resetInfo.steps==0);clear cleaner;
 model='ih_temperature_bh_public';path=fullfile(outdir,[model '.slx']);
 path=radia.simulink.buildTemperatureBHResearchModel(cfg,path,CurrentA=1e4,AmbientK=300,Samples=5,Overwrite=true);
 load_system(path);closer=onCleanup(@()close_system(model,0));
 one=sim(model,'ReturnWorkspaceOutputs','on');two=sim(model,'ReturnWorkspaceOutputs','on');
 assert(isequal(one.temperature.Data,two.temperature.Data)&&isequal(one.heat.Data,two.heat.Data));
 assert(norm(one.temperature.Data-t,'fro')<1e-9&&norm(one.heat.Data-q,'fro')<1e-9);
 impedanceData=complex(one.zsReal.Data,one.zsImag.Data);save(fullfile(outdir,'simulink-diagnostic.mat'),'impedanceData','z');
 assert(isequal(one.zsReal.Data,two.zsReal.Data)&&isequal(one.zsImag.Data,two.zsImag.Data)&&norm(impedanceData-z,'fro')<1e-15);
 set_param([model '/Current'],'Value','realmax');failed=false;
 try
  sim(model,'ReturnWorkspaceOutputs','on');
 catch
  failed=true;
 end
 assert(failed&&radia_cuda_ih('live_count')==0);
 set_param([model '/Current'],'Value','1e4');save_system(model,path);
 assert(radia_cuda_ih('live_count')==0);
 assert(norm(z(end,:)-z(2,:))>0);
 result.impedance_real_ohm=real(z);result.impedance_imag_ohm=imag(z);
 save(fullfile(outdir,'accepted.mat'),'t','q','z','em');
 result.accepted=true;result.temperature_K=t;result.heat_W=q;result.info=info;result.thermal_balance_relative=balance;
 result.failed_candidate_preserved=true;result.reset=true;result.repeat_identical=true;result.save_close_reopen=true;
 result.simulink_failure_released=true;
 result.invalid_material_rejected=true;result.out_of_temperature_range_preserved=true;
 clear closer;
catch e
 result.error=getReport(e,'extended','hyperlinks','off');
end
result.matlab=version;f=fopen(fullfile(outdir,'result.json'),'w');assert(f>=0);writer=onCleanup(@()fclose(f));fprintf(f,'%s',jsonencode(result));
end
