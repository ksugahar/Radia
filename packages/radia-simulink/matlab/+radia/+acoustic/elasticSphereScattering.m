function result = elasticSphereScattering(wavenumber, radius, points, options)
%ELASTICSPHERESCATTERING Readable Faran elastic-sphere validation reference.
arguments
    wavenumber (1,1) double {mustBePositive,mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.LongitudinalSpeed (1,1) double {mustBePositive,mustBeFinite} = 2
    options.ShearSpeed (1,1) double {mustBeNonnegative,mustBeFinite} = 1
    options.DensityRatio (1,1) double {mustBePositive,mustBeFinite} = 1.5
    options.Terms (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
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

terms = options.Terms;
if terms == 0
    terms = ceil(wavenumber * radius) + 10;
end
if terms > 512
    error("radia:acoustic:Terms", "Automatic Terms exceeds 512.");
end
cosine = points(:,3) ./ distance;
coefficients = zeros(terms + 1, 1);
for order = 0:terms
    coefficients(order + 1) = elasticCoefficient( ...
        order, wavenumber, radius, options.LongitudinalSpeed, ...
        options.ShearSpeed, options.DensityRatio);
end

scattered = zeros(size(distance));
lastMode = zeros(size(distance));
for order = 0:terms
    lastMode = (1i^order) * (2 * order + 1) * ...
        coefficients(order + 1) * ...
        sphH(order, wavenumber * distance) .* legendreP(order, cosine);
    scattered = scattered + lastMode;
end
incident = exp(1i * wavenumber .* points(:,3));
result = struct( ...
    "backend", "matlab-reference", ...
    "kind", "elastic_solid_sphere_faran_scattering_series", ...
    "wavenumber", wavenumber, ...
    "radius", radius, ...
    "longitudinal_speed", options.LongitudinalSpeed, ...
    "shear_speed", options.ShearSpeed, ...
    "density_ratio", options.DensityRatio, ...
    "terms", terms, ...
    "truncation_tail", max(abs(lastMode)), ...
    "incident", incident, ...
    "scattered", scattered, ...
    "total", incident + scattered);
end

function coefficient = elasticCoefficient( ...
    order, wavenumber, radius, longitudinalSpeed, shearSpeed, densityRatio)
omega = wavenumber;
kLongitudinal = omega / longitudinalSpeed;
x = wavenumber * radius;
xl = kLongitudinal * radius;
fluidFactor = wavenumber / omega^2;
mu = densityRatio * shearSpeed^2;
lameLambda = densityRatio * (longitudinalSpeed^2 - 2 * shearSpeed^2);

if shearSpeed == 0
    matrix = [ ...
        fluidFactor * sphHDerivative(order, x), ...
        -kLongitudinal * sphJDerivative(order, xl); ...
        sphH(order, x), ...
        -lameLambda * kLongitudinal^2 * sphJ(order, xl)];
    rhs = [-fluidFactor * sphJDerivative(order, x); -sphJ(order, x)];
    solution = scaledSolve(matrix, rhs);
    coefficient = solution(1);
    return
end

kTransverse = omega / shearSpeed;
xt = kTransverse * radius;
angular = order * (order + 1);
urA = kLongitudinal * sphJDerivative(order, xl);
urB = angular / radius * sphJ(order, xt);
durA = kLongitudinal^2 * sphJSecondDerivative(order, xl);
durB = angular * ( ...
    kTransverse * sphJDerivative(order, xt) / radius - ...
    sphJ(order, xt) / radius^2);
srrA = -lameLambda * kLongitudinal^2 * sphJ(order, xl) + 2 * mu * durA;
srrB = 2 * mu * durB;
vaA = sphJ(order, xl) / radius;
vaB = sphJ(order, xt) / radius + kTransverse * sphJDerivative(order, xt);
vpA = -sphJ(order, xl) / radius^2 + ...
    kLongitudinal * sphJDerivative(order, xl) / radius;
vpB = -sphJ(order, xt) / radius^2 + ...
    kTransverse * sphJDerivative(order, xt) / radius + ...
    kTransverse^2 * sphJSecondDerivative(order, xt);
srtA = mu * (urA / radius + vpA - vaA / radius);
srtB = mu * (urB / radius + vpB - vaB / radius);
matrix = [ ...
    fluidFactor * sphHDerivative(order, x), -urA, -urB; ...
    sphH(order, x), srrA, srrB; ...
    0, srtA, srtB];
rhs = [-fluidFactor * sphJDerivative(order, x); -sphJ(order, x); 0];
solution = scaledSolve(matrix, rhs);
coefficient = solution(1);
end

function solution = scaledSolve(matrix, rhs)
% Scale the modal traction system without changing its solution. High-order
% spherical Bessel terms span many decades and otherwise trigger misleading
% near-singular warnings in MATLAB's generic backslash diagnostics.
columnScale = max(abs(matrix), [], 1);
columnScale(columnScale == 0) = 1;
scaledMatrix = matrix ./ columnScale;
rowScale = max(abs(scaledMatrix), [], 2);
rowScale(rowScale == 0) = 1;
scaledSolution = (scaledMatrix ./ rowScale) \ (rhs ./ rowScale);
solution = scaledSolution ./ columnScale.';
end

function value = sphJ(order, argument)
value = sqrt(pi ./ (2 * argument)) .* besselj(order + 0.5, argument);
end

function value = sphH(order, argument)
value = sqrt(pi ./ (2 * argument)) .* ...
    (besselj(order + 0.5, argument) + 1i * bessely(order + 0.5, argument));
end

function value = sphJDerivative(order, argument)
value = sphJ(order - 1, argument) - ...
    (order + 1) ./ argument .* sphJ(order, argument);
end

function value = sphHDerivative(order, argument)
value = sphH(order - 1, argument) - ...
    (order + 1) ./ argument .* sphH(order, argument);
end

function value = sphJSecondDerivative(order, argument)
value = -(2 ./ argument) .* sphJDerivative(order, argument) - ...
    (1 - order * (order + 1) ./ argument.^2) .* sphJ(order, argument);
end
