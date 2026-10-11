function result = propagateVariationalMap(lengthsM, APerM, options)
%PROPAGATEVARIATIONALMAP Native canonical R/T/U transfer-map propagation.
%   RESULT = radia.beam.propagateVariationalMap(LENGTHSM, APERM) propagates
%   the piecewise-constant linear generators APERM, shaped 6-by-6-by-N.
%
%   Name-value arguments F2PerM and F3PerM add the symmetric local equation
%   jets shaped 6-by-6-by-6-by-N and 6-by-6-by-6-by-6-by-N. The map uses
%   u_out = R*u + T[u,u]/2 + U[u,u,u]/6 in coordinate order
%   (x, px/p0, y, py/p0, sigma, delta).

arguments
    lengthsM double {mustBeReal,mustBeFinite,mustBeNonempty}
    APerM double {mustBeReal,mustBeFinite,mustBeNonempty}
    options.F2PerM double {mustBeReal,mustBeFinite} = []
    options.F3PerM double {mustBeReal,mustBeFinite} = []
    options.Names = strings(0,1)
    options.MaximumOrder (1,1) double {mustBeInteger,mustBePositive} = 3
    options.MaximumStepM (1,1) double {mustBeFinite,mustBePositive} = 1e-3
    options.MaximumSteps (1,1) double {mustBeInteger,mustBePositive} = 1e6
    options.MaximumRegionPairs (1,1) double {mustBeInteger,mustBePositive} = 1e5
    options.InputSymmetryTolerance (1,1) double ...
        {mustBeFinite,mustBeNonnegative} = 1e-12
end

lengthsM = double(lengthsM(:));
segmentCount = numel(lengthsM);
if any(lengthsM <= 0)
    error("radia:beam:InvalidLength", ...
        "Every segment length must be positive.");
end
if options.MaximumOrder > 3
    error("radia:beam:InvalidOrder", ...
        "MaximumOrder must be 1, 2, or 3.");
end
requireLeadingShape(APerM,2,segmentCount,"APerM");
if ~isempty(options.F2PerM)
    requireLeadingShape(options.F2PerM,3,segmentCount,"F2PerM");
end
if ~isempty(options.F3PerM)
    requireLeadingShape(options.F3PerM,4,segmentCount,"F3PerM");
end

config.schema = 'radia.beam.variational-map.v1';
config.lengths_m = lengthsM;
config.A_per_m = double(APerM);
if ~isempty(options.F2PerM)
    config.F2_per_m = double(options.F2PerM);
end
if ~isempty(options.F3PerM)
    config.F3_per_m = double(options.F3PerM);
end
if ~isempty(options.Names)
    names = string(options.Names(:));
    if numel(names) ~= segmentCount || any(ismissing(names))
        error("radia:beam:InvalidNames", ...
            "Names must contain one nonmissing entry per segment.");
    end
    config.names = cellstr(names);
end
config.maximum_order = options.MaximumOrder;
config.maximum_step_m = options.MaximumStepM;
config.maximum_steps = options.MaximumSteps;
config.maximum_region_pairs = options.MaximumRegionPairs;
config.input_symmetry_tolerance = options.InputSymmetryTolerance;

result = radia.internal.callMex( ...
    'beam.transfer.propagate_variational',config);
end

function requireLeadingShape(value, leadingRank, segmentCount, name)
expectedElements = segmentCount*6^leadingRank;
valid = numel(value) == expectedElements;
for axis = 1:leadingRank
    valid = valid && size(value,axis) == 6;
end
if ~valid
    error("radia:beam:InvalidShape", ...
        "%s must have %d leading dimensions of size 6 and %d segments.", ...
        name,leadingRank,segmentCount);
end
end
