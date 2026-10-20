function force_N = virtualWorkForce(positions_m, energy_J, energyKind)
%VIRTUALWORKFORCE Differentiate an energy table into force samples.
%   Fixed-current coenergy uses +dW'/dx. Fixed-flux stored energy uses
%   -dW/dx. Every sample, endpoints included, comes from the quadratic
%   through three consecutive rows, so the whole table is second-order
%   accurate and the spacing need not be uniform. Matches
%   radia.force.virtual_work_force_from_displacement_samples.

if nargin < 3 || isempty(energyKind)
    energyKind = "coenergy";
end
positions = localTable(positions_m, "positions_m");
energy = localTable(energy_J, "energy_J");
if numel(positions) ~= numel(energy)
    error("radia:force:SampleCount", "positions_m and energy_J must have the same length");
end
if any(diff(positions) <= 0)
    error("radia:force:Positions", "positions_m must be strictly increasing");
end
derivative = localThreePointDerivative(positions, energy);
key = lower(replace(replace(strtrim(string(energyKind)), "-", "_"), " ", "_"));
if any(key == ["coenergy", "magnetic_coenergy", "constant_current", "w_prime"])
    signValue = 1.0;
elseif any(key == ["stored_energy", "field_energy", "magnetic_energy", "constant_flux"])
    signValue = -1.0;
else
    error("radia:force:EnergyKind", "energyKind must be coenergy/constant_current or stored_energy/constant_flux");
end
force_N = signValue * derivative;
end

function derivative = localThreePointDerivative(nodes, values)
%LOCALTHREEPOINTDERIVATIVE Second-order derivative on an arbitrary 1D grid.
%   Every sample uses the quadratic through three consecutive nodes, so the
%   endpoints are second-order too and the spacing need not be uniform. On a
%   uniform grid this reduces to the classic (f(i+1)-f(i-1))/(2h) interior and
%   (-3f0+4f1-f2)/(2h) / (f(n-2)-4f(n-1)+3f(n))/(2h) endpoint stencils.
step = diff(nodes);
back = step(1:end-1);
forward = step(2:end);
derivative = zeros(size(values));
derivative(2:end-1) = ...
    -forward ./ (back .* (back + forward)) .* values(1:end-2) ...
    + (forward - back) ./ (back .* forward) .* values(2:end-1) ...
    + back ./ (forward .* (back + forward)) .* values(3:end);
first = step(1);
second = step(2);
derivative(1) = ...
    -(2*first + second) / (first * (first + second)) * values(1) ...
    + (first + second) / (first * second) * values(2) ...
    - first / (second * (first + second)) * values(3);
last = step(end);
prior = step(end-1);
derivative(end) = ...
    last / (prior * (last + prior)) * values(end-2) ...
    - (last + prior) / (last * prior) * values(end-1) ...
    + (2*last + prior) / (last * (last + prior)) * values(end);
end

function table = localTable(value, name)
if ~isnumeric(value) || ~isreal(value) || any(~isfinite(value), "all") || ~isvector(value) || numel(value) < 3
    error("radia:force:Table", "%s must be a finite real vector with at least three samples", name);
end
table = reshape(double(value), [], 1);
end
