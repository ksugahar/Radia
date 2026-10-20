function hdivGPUWorkflowSFunction(block)
% Shared GPU Newton owner; one explicit rising trigger schedules one solve.
% Outputs are held and become visible on the next discrete sample (Update).
% No coil voltage/inductance model, Python call, or CPU fallback is implemented.
block.NumDialogPrms=1;block.DialogPrmsTunable={'Nontunable'};
c=radia.simulink.validateHDivGPUWorkflowConfig(block.DialogPrm(1).Data);
block.NumInputPorts=2;block.NumOutputPorts=5;
for k=1:2
    block.InputPort(k).Dimensions=1;block.InputPort(k).Complexity='Real';
    block.InputPort(k).DirectFeedthrough=false;
end
block.InputPort(1).DatatypeID=8;block.InputPort(2).DatatypeID=0;
for k=1:5
    block.OutputPort(k).DatatypeID=0;block.OutputPort(k).Complexity='Real';
end
block.OutputPort(1).Dimensions=c.NumDoF;block.OutputPort(2).Dimensions=3;
for k=3:5,block.OutputPort(k).Dimensions=1;end
block.OutputPort(5).DatatypeID=6;
block.SampleTimes=[double(c.SampleTime) 0];
block.SimStateCompliance='DisallowSimState';
block.RegBlockMethod('PostPropagationSetup',@dworks);
block.RegBlockMethod('Start',@start);
block.RegBlockMethod('InitializeConditions',@initialize);
block.RegBlockMethod('Outputs',@outputs);
block.RegBlockMethod('Update',@update);
block.RegBlockMethod('Terminate',@terminate);
end
function dworks(block)
c=block.DialogPrm(1).Data;block.NumDworks=6;
names={'previous_trigger','coefficients','mean_magnetization','residual','iterations','status'};
dims=[1 c.NumDoF 3 1 1 1];types=[8 0 0 0 0 6];
for k=1:6
    block.Dwork(k).Name=names{k};block.Dwork(k).Dimensions=dims(k);
    block.Dwork(k).DatatypeID=types(k);block.Dwork(k).Complexity='Real';
    block.Dwork(k).UsedAsDiscState=true;
end
end
function initialize(block)
block.Dwork(1).Data=false;block.Dwork(2).Data=zeros(block.DialogPrm(1).Data.NumDoF,1);
block.Dwork(3).Data=zeros(3,1);block.Dwork(4).Data=NaN;
block.Dwork(5).Data=0;block.Dwork(6).Data=int32(0);
end
function start(block)
c=radia.simulink.validateHDivGPUWorkflowConfig(block.DialogPrm(1).Data);
assert(isempty(get_param(block.BlockHandle,'UserData')), ...
    'radia:simulink:HDivOwner','The block already has user data; refusing to overwrite');
if ~isfolder(c.RunDirectory),mkdir(c.RunDirectory);end
c.RunDirectory=tempname(c.RunDirectory);mkdir(c.RunDirectory);
try
    owner=radia.gpu.PreparedHDivNonlinear(c.MeshPath,c.Order,c.AppliedField,c.BHTable, ...
        Assembly=c.Assembly,HostSetup=c.HostSetup,Preconditioner=c.Preconditioner);
    try
        set_param(block.BlockHandle,'UserData', ...
            struct('schema','radia.hdiv-gpu-sfun.owner.v1','owner',owner,'runDirectory',c.RunDirectory), ...
            'UserDataPersistent','off');
    catch e
        delete(owner);rethrow(e);
    end
    journal(c,'prepared',struct('backend','gpu','preconditioner',c.Preconditioner));
catch e
    terminate(block);journal(c,'failed',struct('identifier',e.identifier,'message',e.message));rethrow(e);
end
end
function outputs(block)
for k=1:5,block.OutputPort(k).Data=block.Dwork(k+1).Data;end
end
function update(block)
trigger=logical(block.InputPort(1).Data);previous=logical(block.Dwork(1).Data);
block.Dwork(1).Data=trigger;
if ~trigger || previous,return;end
c=block.DialogPrm(1).Data;
try
    state=get_param(block.BlockHandle,'UserData');
    assert(isstruct(state) && strcmp(state.schema,'radia.hdiv-gpu-sfun.owner.v1') && ...
        isvalid(state.owner),'radia:simulink:HDivOwner','Missing owned GPU workflow');
    c.RunDirectory=state.runDirectory;
    answer=state.owner.solve(double(block.InputPort(2).Data),block.Dwork(2).Data);
    assert(numel(answer.coefficients)==c.NumDoF,'radia:simulink:HDivDimensions', ...
        'Checked FE coefficient count differs from the declared Simulink port');
    block.Dwork(2).Data=answer.coefficients;block.Dwork(3).Data=answer.mean_magnetization;
    block.Dwork(4).Data=answer.nonlinear_relative_residual;
    block.Dwork(5).Data=answer.nonlinear_iterations;block.Dwork(6).Data=int32(2);
    journal(c,'passed',struct('backend','gpu','field_scale',answer.field_scale, ...
        'coefficients',answer.coefficients,'mean_magnetization',answer.mean_magnetization, ...
        'true_relative_residual',answer.nonlinear_relative_residual, ...
        'newton_iterations',answer.nonlinear_iterations, ...
        'setup_upload_bytes',answer.material_stats.setup_upload_bytes));
catch e
    block.Dwork(6).Data=int32(-1);terminate(block);
    journal(c,'failed',struct('identifier',e.identifier,'message',e.message));rethrow(e);
end
end
function terminate(block)
state=get_param(block.BlockHandle,'UserData');
if isempty(state),return;end
assert(isstruct(state) && isfield(state,'schema') && ...
    strcmp(state.schema,'radia.hdiv-gpu-sfun.owner.v1'), ...
    'radia:simulink:HDivOwner','Refusing to clear unrelated block data');
set_param(block.BlockHandle,'UserData',[],'UserDataPersistent','off');
delete(state.owner);
end
function journal(c,status,value)
value.status=status;value.schema='radia.hdiv-gpu-sfun.result.v1';
tmp=[tempname(c.RunDirectory) '.json'];fid=fopen(tmp,'w');
assert(fid>=0,'radia:simulink:HDivArtifact','Cannot open owned result');
try
    fprintf(fid,'%s',jsonencode(value));assert(fclose(fid)==0);fid=-1;
    [ok,msg]=movefile(tmp,fullfile(c.RunDirectory,'result.json'),'f');
    assert(ok,'radia:simulink:HDivArtifact','%s',msg);
catch e
    if fid>=0,fclose(fid);end
    rethrow(e);
end
fid=fopen(fullfile(c.RunDirectory,'run.log'),'a');
assert(fid>=0,'radia:simulink:HDivArtifact','Cannot append owned log');
fprintf(fid,'%s %s\n',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),status);fclose(fid);
end
