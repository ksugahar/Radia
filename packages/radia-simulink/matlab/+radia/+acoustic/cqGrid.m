function grid = cqGrid(sampleCount, timeStep, options)
%CQGRID Lubich convolution-quadrature Laplace/wavenumber grid.
arguments
    sampleCount (1,1) double {mustBeInteger,mustBePositive,mustBeFinite}
    timeStep (1,1) double {mustBePositive,mustBeFinite}
    options.SoundSpeed (1,1) double {mustBePositive,mustBeFinite} = 1
    options.Method (1,1) string {mustBeMember(options.Method,["BDF1","BDF2"])} = "BDF2"
end
if sampleCount > 1048576
    error("radia:acoustic:SampleCount", ...
        "SampleCount must be a positive practical integer.");
end
index = 0:sampleCount-1;
radius = eps^(1 / (2 * sampleCount));
zeta = radius * exp(-2 * 1i * pi * index / sampleCount);
nodes = radia.acoustic.bdfDelta(zeta, options.Method) / timeStep;
grid = struct( ...
    "backend", "matlab-native", ...
    "cq_radius", radius, ...
    "zeta", zeta, ...
    "cq_nodes", nodes, ...
    "cq_wavenumbers", 1i * nodes / options.SoundSpeed);
end
