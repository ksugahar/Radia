function result = fluidSphereScattering(wavenumber, radius, points, options)
%FLUIDSPHERESCATTERING Readable Anderson fluid-sphere validation reference.
arguments
    wavenumber (1,1) double {mustBePositive,mustBeFinite}
    radius (1,1) double {mustBePositive,mustBeFinite}
    points (:,3) double {mustBeFinite}
    options.InteriorWavenumber (1,1) double {mustBePositive,mustBeFinite} = wavenumber
    options.DensityRatio (1,1) double {mustBePositive,mustBeFinite} = 1
    options.Terms (1,1) double {mustBeInteger,mustBeFinite} = -1
end
validateTerms(options.Terms, true);
if isempty(points)
    error("radia:acoustic:Points", "Points must be a nonempty N-by-3 matrix.");
end

k0 = wavenumber;
k1 = options.InteriorWavenumber;
densityRatio = options.DensityRatio;
distance = sqrt(sum(points.^2, 2));
safeDistance = max(distance, 1e-30);
cosine = points(:,3) ./ safeDistance;
inside = distance <= radius * (1 + 1e-12);
requested = max(options.Terms, 0);
terms = max(requested, ...
    ceil(max(k0 * max([radius; distance]), k1 * radius)) + 12);
validateTerms(terms, false);

x0 = k0 * radius;
x1 = k1 * radius;
total = zeros(size(distance));
lastMode = zeros(size(distance));
for order = 0:terms
    polynomial = legendreP(order, cosine);
    incidentCoefficient = (1i^order) * (2 * order + 1);
    j0 = sphJ(order, x0);
    h0 = sphH(order, x0);
    j1 = sphJ(order, x1);
    beta = (k1 / densityRatio) * sphJDerivative(order, x1) / j1;
    scatteredCoefficient = -incidentCoefficient * ...
        (k0 * sphJDerivative(order, x0) - beta * j0) / ...
        (k0 * sphHDerivative(order, x0) - beta * h0);
    interiorCoefficient = ...
        (incidentCoefficient * j0 + scatteredCoefficient * h0) / j1;

    mode = zeros(size(distance));
    mode(inside) = interiorCoefficient * ...
        sphJ(order, k1 * safeDistance(inside)) .* polynomial(inside);
    mode(~inside) = (incidentCoefficient * ...
        sphJ(order, k0 * safeDistance(~inside)) + ...
        scatteredCoefficient * sphH(order, k0 * safeDistance(~inside))) .* ...
        polynomial(~inside);
    total = total + mode;
    lastMode = mode;
end
result = struct( ...
    "backend", "matlab-reference", ...
    "kind", "fluid_sphere_transmission_scattering_series", ...
    "wavenumber", k0, ...
    "interior_wavenumber", k1, ...
    "density_ratio", densityRatio, ...
    "radius", radius, ...
    "terms", terms, ...
    "truncation_tail", max(abs(lastMode)), ...
    "incident", exp(1i * k0 .* points(:,3)), ...
    "total", total, ...
    "inside_mask", inside);
end

function validateTerms(terms, allowAutomatic)
if (allowAutomatic && terms == -1) || (terms >= 0 && terms <= 512)
    return
end
error("radia:acoustic:Terms", "Terms are outside the supported range.");
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
