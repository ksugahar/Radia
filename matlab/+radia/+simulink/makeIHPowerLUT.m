function lut = makeIHPowerLUT(breakpoints, power_W, options)
%MAKEIHPOWERLUT Validate a Radia-to-Simulink power lookup table.
%   LUT = radia.simulink.makeIHPowerLUT(BREAKPOINTS,POWER_W) returns the
%   portable contract used by an n-D Lookup Table block. BREAKPOINTS is a
%   cell array of monotonically increasing physical grids, for example
%   {position_rad, drive_A, temperature_K}. POWER_W has the corresponding
%   ndgrid shape. The first dimensions may represent motion position and
%   speed; the table can therefore encode moving-workpiece data without
%   changing the thermal plant.

arguments
    breakpoints (1,:) cell
    power_W double
    options.InputNames (1,:) string = string.empty(1, 0)
    options.Extrapolation (1,1) string = "clip"
    options.Source (1,1) string = "Radia VIM/FEM/SIBC/ESIM"
end

nDim = numel(breakpoints);
if nDim < 1
    error("radia:simulink:EmptyLUT", "at least one breakpoint is required.");
end
if options.InputNames.isempty()
    options.InputNames = "input" + (1:nDim);
end
if ~ismember(options.Extrapolation, ["clip", "error"])
    error("radia:simulink:LUTExtrapolation", ...
        "Extrapolation must be clip or error.");
end
if numel(options.InputNames) ~= nDim
    error("radia:simulink:LUTNames", "InputNames must match breakpoint count.");
end

shape = zeros(1, nDim);
for k = 1:nDim
    grid = breakpoints{k};
    if ~isvector(grid) || isempty(grid) || any(~isfinite(grid)) || ...
            any(diff(grid) <= 0)
        error("radia:simulink:LUTGrid", ...
            "each breakpoint must be finite, nonempty, and strictly increasing.");
    end
    breakpoints{k} = grid(:).';
    shape(k) = numel(grid);
end
actualShape = size(power_W);
actualShape(end + 1:nDim) = 1;
if (nDim == 1 && numel(power_W) ~= shape(1)) || ...
        (nDim > 1 && ~isequal(actualShape(1:nDim), shape))
    error("radia:simulink:LUTShape", ...
        "power_W must have ndgrid shape [%s].", join(string(shape), ", "));
end
if any(~isfinite(power_W(:))) || any(power_W(:) < 0)
    error("radia:simulink:LUTPower", "power_W must be finite and nonnegative.");
end
power_W = reshape(power_W, shape);

lut = struct( ...
    "schema", "radia.ih.simulink.power_lut.v1", ...
    "breakpoints", {breakpoints}, ...
    "table_power_W", power_W, ...
    "input_names", options.InputNames, ...
    "output_name", "power_W", ...
    "extrapolation", options.Extrapolation, ...
    "source", options.Source, ...
    "table_dimensions", shape);
end
