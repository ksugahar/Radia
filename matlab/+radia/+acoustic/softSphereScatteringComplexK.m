function scattered = softSphereScatteringComplexK(wavenumber, radius, points, options)
%SOFTSPHERESCATTERINGCOMPLEXK Native complex-wavenumber sphere reference.
arguments
    wavenumber (1,1) double {mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.Terms (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 28
end
if wavenumber == 0
    error("radia:acoustic:Wavenumber", "Wavenumber must be nonzero.");
end
scattered = radia.internal.callMex("acoustic.soft_sphere_complex_k", ...
    wavenumber, radius, points, options.Terms);
end
