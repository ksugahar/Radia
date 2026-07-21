function grid = cqGrid(sampleCount, timeStep, options)
%CQGRID Native Lubich convolution-quadrature Laplace/wavenumber grid.
arguments
    sampleCount (1,1) double {mustBeInteger,mustBePositive,mustBeFinite}
    timeStep (1,1) double {mustBePositive,mustBeFinite}
    options.SoundSpeed (1,1) double {mustBePositive,mustBeFinite} = 1
    options.Method (1,1) string {mustBeMember(options.Method,["BDF1","BDF2"])} = "BDF2"
end
grid = radia.internal.callMex("acoustic.cq_grid", sampleCount, timeStep, ...
    options.SoundSpeed, char(options.Method));
grid.backend = "native-mex";
end
