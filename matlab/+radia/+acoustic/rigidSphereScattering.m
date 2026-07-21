function result = rigidSphereScattering(wavenumber, radius, points, options)
%RIGIDSPHERESCATTERING Native sound-hard sphere scattering reference.
arguments
    wavenumber (1,1) double {mustBePositive,mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.Terms (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 1
end
result = radia.internal.callMex("acoustic.rigid_sphere", ...
    wavenumber, radius, points, options.Terms);
result.backend = "native-mex";
end
