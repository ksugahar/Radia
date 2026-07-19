function power_W = evaluateIHPowerLUT(lut, inputs)
%EVALUATEIHPOWERLUT Evaluate a validated IH power LUT at row-wise inputs.
%   INPUTS is N-by-D in the order stored in LUT.input_names. With the
%   default clip policy, values outside the Radia training range are clipped
%   instead of silently extrapolated.

arguments
    lut (1,1) struct
    inputs double {mustBeFinite}
end

if ~isfield(lut, "schema") || lut.schema ~= "radia.ih.simulink.power_lut.v1"
    error("radia:simulink:InvalidLUT", "unsupported IH power LUT schema.");
end
if size(inputs, 2) ~= numel(lut.breakpoints)
    error("radia:simulink:LUTInputSize", ...
        "inputs must have one column per LUT breakpoint.");
end

query = inputs;
for k = 1:size(query, 2)
    grid = lut.breakpoints{k};
    if lut.extrapolation == "clip"
        query(:, k) = min(max(query(:, k), grid(1)), grid(end));
    elseif any(query(:, k) < grid(1) | query(:, k) > grid(end))
        error("radia:simulink:LUTExtrapolation", ...
            "input column %d is outside its Radia training range.", k);
    end
end

F = griddedInterpolant(lut.breakpoints, lut.table_power_W, "linear", "nearest");
args = num2cell(query, 1);
power_W = F(args{:});
power_W = power_W(:);
end
