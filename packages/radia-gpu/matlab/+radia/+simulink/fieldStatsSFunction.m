function fieldStatsSFunction(block)
%FIELDSTATSSFUNCTION Reduce an N-wide field vector to min / mean / max.
%   Readable Level-2 MATLAB S-Function behind the Field Stats block
%   (top-level wrapper radia_field_stats_sfun).  Input: one dynamically
%   sized real vector (a distributed field such as temperature or heat
%   density).  Outputs: the three scalars min(u), mean(u), max(u).
%
%   The reduction is deliberately .m code instead of a primitive-block
%   web (MinMax / Sum / Width / Divide) so the computation is diffable,
%   testable, and editable without opening Simulink.  The mean is the
%   arithmetic mean over DOFs -- a display aid; volume-weighted means
%   stay owned by the result artifacts.
setup(block);
end

function setup(block)
block.NumInputPorts = 1;
block.NumOutputPorts = 3;
block.SetPreCompInpPortInfoToDynamic;
block.InputPort(1).DatatypeID = 0;
block.InputPort(1).Complexity = "Real";
block.InputPort(1).DirectFeedthrough = true;
for k = 1:3
    block.OutputPort(k).Dimensions = 1;
    block.OutputPort(k).DatatypeID = 0;
    block.OutputPort(k).Complexity = "Real";
end
block.SampleTimes = [-1 0];
block.SimStateCompliance = "DefaultSimState";
block.RegBlockMethod("SetInputPortDimensions", @setInputPortDimensions);
block.RegBlockMethod("SetInputPortSamplingMode", @setInputPortSamplingMode);
block.RegBlockMethod("Outputs", @outputs);
end

function setInputPortDimensions(block, port, dimensions)
block.InputPort(port).Dimensions = dimensions;
end

function setInputPortSamplingMode(block, port, mode)
block.InputPort(port).SamplingMode = mode;
for k = 1:3
    block.OutputPort(k).SamplingMode = mode;
end
end

function outputs(block)
u = block.InputPort(1).Data;
block.OutputPort(1).Data = min(u);
block.OutputPort(2).Data = mean(u);
block.OutputPort(3).Data = max(u);
end
