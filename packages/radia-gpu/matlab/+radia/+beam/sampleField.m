function sample = sampleField(field,positionM,options)
%SAMPLEFIELD Evaluate a native beam field at one Cartesian point.
arguments
    field (1,1) struct
    positionM (1,3) double {mustBeReal,mustBeFinite}
    options.TimeS (1,1) double {mustBeReal,mustBeFinite} = 0
end
sample = radia.internal.callMex( ...
    'beam.field.sample',field,double(positionM),double(options.TimeS));
end
