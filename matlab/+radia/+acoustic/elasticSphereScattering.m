function result = elasticSphereScattering(wavenumber, radius, points, options)
%ELASTICSPHERESCATTERING Native Faran elastic-sphere scattering reference.
arguments
    wavenumber (1,1) double {mustBePositive,mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.LongitudinalSpeed (1,1) double {mustBePositive,mustBeFinite} = 2
    options.ShearSpeed (1,1) double {mustBeNonnegative,mustBeFinite} = 1
    options.DensityRatio (1,1) double {mustBePositive,mustBeFinite} = 1.5
    options.Terms (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
end
result = radia.internal.callMex("acoustic.elastic_sphere", ...
    wavenumber, radius, points, options.LongitudinalSpeed, ...
    options.ShearSpeed, options.DensityRatio, options.Terms);
result.backend = "native-mex";
end
