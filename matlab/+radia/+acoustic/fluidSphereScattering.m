function result = fluidSphereScattering(wavenumber, radius, points, options)
%FLUIDSPHERESCATTERING Native penetrable-fluid sphere transmission reference.
arguments
    wavenumber (1,1) double {mustBePositive,mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.InteriorWavenumber (1,1) double {mustBePositive,mustBeFinite} = wavenumber
    options.DensityRatio (1,1) double {mustBePositive,mustBeFinite} = 1
    options.Terms (1,1) double {mustBeInteger,mustBeFinite} = -1
end
if options.Terms < -1
    error("radia:acoustic:Terms", "Terms must be -1 (automatic) or nonnegative.");
end
result = radia.internal.callMex("acoustic.fluid_sphere", ...
    wavenumber, radius, points, options.InteriorWavenumber, ...
    options.DensityRatio, options.Terms);
result.backend = "native-mex";
end
