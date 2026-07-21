function result = softSphereScattering(wavenumber, radius, points, options)
%SOFTSPHERESCATTERING Native sound-soft sphere scattering reference.
arguments
    wavenumber (1,1) double {mustBePositive,mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.Terms (1,1) double {mustBeInteger,mustBeFinite} = -1
end
if options.Terms < -1
    error("radia:acoustic:Terms", "Terms must be -1 (automatic) or nonnegative.");
end
result = radia.internal.callMex("acoustic.soft_sphere", ...
    wavenumber, radius, points, options.Terms);
result.backend = "native-mex";
end
