function winding = makeWindingSpec(options)
%MAKEWINDINGSPEC Construct one winding dictionary value.

arguments
    options.Regions (:,1) string
    options.RegionPolarity (:,1) double = zeros(0,1)
    options.Turns (1,1) double {mustBeInteger,mustBePositive}
    options.Polarity (1,1) double {mustBeMember(options.Polarity,[-1 1])} = 1
    options.ParallelPaths (1,1) double {mustBeInteger,mustBePositive} = 1
    options.Resistance_ohm (1,1) double {mustBeFinite,mustBeNonnegative} = 0
    options.PositiveTerminal (1,1) string = "p"
    options.NegativeTerminal (1,1) string = "n"
    options.Description (1,1) string = ""
end
regions = string(options.Regions(:));
if isempty(regions) || any(strlength(regions)==0) || numel(unique(regions))~=numel(regions)
    error("radia:simulink:WindingRegions", ...
        "Regions must be a nonempty list of unique .vol region names.");
end
regionPolarity = double(options.RegionPolarity(:));
if isempty(regionPolarity),regionPolarity=ones(numel(regions),1);end
if numel(regionPolarity)~=numel(regions) || any(~ismember(regionPolarity,[-1 1]))
    error("radia:simulink:WindingRegionPolarity", ...
        "RegionPolarity must contain one +1 or -1 value per winding region.");
end
if numel(unique(lower(regions))) ~= numel(regions)
    error("radia:simulink:WindingRegionCaseCollision", ...
        "Winding regions that differ only by case are not allowed.");
end
if strlength(options.PositiveTerminal)==0 || strlength(options.NegativeTerminal)==0 || ...
        options.PositiveTerminal==options.NegativeTerminal
    error("radia:simulink:WindingTerminals", ...
        "PositiveTerminal and NegativeTerminal must be distinct nonempty names.");
end
winding = struct("regions",regions,"region_polarity",int8(regionPolarity), ...
    "turns",uint32(options.Turns),"polarity",int8(options.Polarity), ...
    "parallel_paths",uint16(options.ParallelPaths), ...
    "resistance_ohm",double(options.Resistance_ohm), ...
    "positive_terminal",options.PositiveTerminal, ...
    "negative_terminal",options.NegativeTerminal, ...
    "description",options.Description);
end
