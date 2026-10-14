function scattered = softSphereScatteringComplexK(wavenumber, radius, points, options)
%SOFTSPHERESCATTERINGCOMPLEXK Readable complex-wavenumber sphere reference.
arguments
    wavenumber (1,1) double {mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.Terms (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 28
end
if wavenumber == 0
    error("radia:acoustic:Wavenumber", "Wavenumber must be nonzero.");
end
if options.Terms > 512
    error("radia:acoustic:Terms", "Terms must be between 0 and 512.");
end
if isempty(points)
    error("radia:acoustic:Points", "Points must be a nonempty N-by-3 matrix.");
end
distance = sqrt(sum(points.^2, 2));
if any(distance < radius * (1 - 1e-9))
    error("radia:acoustic:Exterior", ...
        "Evaluation points must lie on or outside the sphere r >= R.");
end
cosine = points(:,3) ./ distance;
scattered = zeros(size(distance));
for order = 0:options.Terms
    coefficient = -(1i^order) * (2 * order + 1) * ...
        sphJ(order, wavenumber * radius) / ...
        sphH(order, wavenumber * radius);
    scattered = scattered + coefficient * ...
        sphH(order, wavenumber * distance) .* legendreP(order, cosine);
end
end

function value = sphJ(order, argument)
value = sqrt(pi ./ (2 * argument)) .* besselj(order + 0.5, argument);
end

function value = sphH(order, argument)
value = sqrt(pi ./ (2 * argument)) .* ...
    (besselj(order + 0.5, argument) + 1i * bessely(order + 0.5, argument));
end
