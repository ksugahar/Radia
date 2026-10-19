function result = softSphereScattering(wavenumber, radius, points, options)
%SOFTSPHERESCATTERING Readable sound-soft sphere validation reference.
arguments
    wavenumber (1,1) double {mustBePositive,mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.Terms (1,1) double {mustBeInteger,mustBeFinite} = -1
end
validateTerms(options.Terms, true);
[distance, cosine] = exteriorPoints(points, radius);
terms = options.Terms;
if terms < 0
    terms = ceil(wavenumber * radius) + 12;
end
validateTerms(terms, false);

scattered = zeros(size(distance));
lastMode = zeros(size(distance));
for order = 0:terms
    coefficient = -(1i^order) * (2 * order + 1) * ...
        sphJ(order, wavenumber * radius) / ...
        sphH(order, wavenumber * radius);
    lastMode = coefficient * sphH(order, wavenumber * distance) .* ...
        legendreP(order, cosine);
    scattered = scattered + lastMode;
end
incident = exp(1i * wavenumber .* points(:,3));
result = struct( ...
    "backend", "matlab-reference", ...
    "kind", "soft_sphere_plane_wave_scattering_series", ...
    "wavenumber", wavenumber, ...
    "radius", radius, ...
    "terms", terms, ...
    "truncation_tail", max(abs(lastMode)), ...
    "scattered", scattered, ...
    "incident", incident, ...
    "total", incident + scattered);
end

function [distance, cosine] = exteriorPoints(points, radius)
if isempty(points)
    error("radia:acoustic:Points", "Points must be a nonempty N-by-3 matrix.");
end
distance = sqrt(sum(points.^2, 2));
if any(distance < radius * (1 - 1e-9))
    error("radia:acoustic:Exterior", ...
        "Evaluation points must lie on or outside the sphere r >= R.");
end
cosine = points(:,3) ./ distance;
end

function validateTerms(terms, allowAutomatic)
if (allowAutomatic && terms == -1) || (terms >= 0 && terms <= 512)
    return
end
if allowAutomatic
    error("radia:acoustic:Terms", ...
        "Terms must be -1 or between 0 and 512.");
end
error("radia:acoustic:Terms", "Terms must be between 0 and 512.");
end

function value = sphJ(order, argument)
value = sqrt(pi ./ (2 * argument)) .* besselj(order + 0.5, argument);
end

function value = sphH(order, argument)
value = sqrt(pi ./ (2 * argument)) .* ...
    (besselj(order + 0.5, argument) + 1i * bessely(order + 0.5, argument));
end
