function force_N = virtualWorkForce(positions_m, energy_J, energyKind)
%VIRTUALWORKFORCE Differentiate an energy table into force samples.
%   Fixed-current coenergy uses +dW'/dx. Fixed-flux stored energy uses
%   -dW/dx. Endpoints use one-sided differences; interiors use centred
%   differences, matching radia.force.

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
derivative = zeros(size(energy));
derivative(1) = (energy(2) - energy(1)) / (positions(2) - positions(1));
derivative(end) = (energy(end) - energy(end-1)) / (positions(end) - positions(end-1));
derivative(2:end-1) = (energy(3:end) - energy(1:end-2)) ./ ...
    (positions(3:end) - positions(1:end-2));
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

function table = localTable(value, name)
if ~isnumeric(value) || ~isreal(value) || any(~isfinite(value), "all") || ~isvector(value) || numel(value) < 3
    error("radia:force:Table", "%s must be a finite real vector with at least three samples", name);
end
table = reshape(double(value), [], 1);
end
