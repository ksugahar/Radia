function ihResidentSFunction(block)
%IHRESIDENTSFUNCTION Retain a checked coupled IH CUDA context across steps.
% Fixed and explicit element-ESIM operator frontends share this state owner.
% Complete geometry/reaction assembly remains a separate research gate.
setup(block);
end
function setup(block)
block.NumDialogPrms=1;
block.DialogPrmsTunable={'Nontunable'};
config=block.DialogPrm(1).Data;
block.NumInputPorts=2;
block.NumOutputPorts=2;
for k=1:2
    block.InputPort(k).Dimensions=1;
    block.InputPort(k).DatatypeID=0;
    block.InputPort(k).Complexity='Real';
    block.InputPort(k).DirectFeedthrough=false;
    block.OutputPort(k).DatatypeID=0;
    block.OutputPort(k).Complexity='Real';
end
block.OutputPort(1).Dimensions=numel(config.initial_temperature_K);
block.OutputPort(2).Dimensions=size(config.loss_projection,1)/config.loss_components;
block.SampleTimes=[config.sample_time_s,0];
block.SimStateCompliance='DisallowSimState';
block.RegBlockMethod('PostPropagationSetup',@postPropagationSetup);
block.RegBlockMethod('Start',@start);
block.RegBlockMethod('InitializeConditions',@initializeConditions);
block.RegBlockMethod('Outputs',@outputs);
block.RegBlockMethod('Update',@update);
block.RegBlockMethod('Terminate',@terminate);
end
function postPropagationSetup(block)
block.NumDworks=2;
for k=1:2
    block.Dwork(k).Name=sprintf('ih_cuda_handle_%d',k);
    block.Dwork(k).Dimensions=1;
    block.Dwork(k).DatatypeID=7;
    block.Dwork(k).Complexity='Real';
    block.Dwork(k).UsedAsDiscState=false;
end
end
function start(block)
terminate(block);
h=radia_cuda_ih('create',block.DialogPrm(1).Data);
block.Dwork(1).Data=uint32(bitand(h,uint64(4294967295)));
block.Dwork(2).Data=uint32(bitshift(h,-32));
end
function h=handle(block)
h=bitor(uint64(block.Dwork(1).Data),bitshift(uint64(block.Dwork(2).Data),32));
end
function initializeConditions(block)
h=handle(block);
if h~=0, radia_cuda_ih('reset',h); end
end
function outputs(block)
[t,q,~]=radia_cuda_ih('outputs',handle(block));
block.OutputPort(1).Data=t;
block.OutputPort(2).Data=q;
end
function update(block)
try
    radia_cuda_ih('update',handle(block),double(block.InputPort(1).Data),double(block.InputPort(2).Data));
catch exception
    terminate(block);
    rethrow(exception);
end
end
function terminate(block)
h=handle(block);
if h==0, return; end
block.Dwork(1).Data=uint32(0);
block.Dwork(2).Data=uint32(0);
radia_cuda_ih('destroy',h);
end
