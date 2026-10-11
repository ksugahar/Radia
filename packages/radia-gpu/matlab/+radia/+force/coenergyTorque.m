function torque_Nm = coenergyTorque(angles_rad, coenergy_J, periodic, period_rad)
%COENERGYTORQUE Differentiate fixed-current coenergy by mechanical angle.
%   Second-order accurate at every sample, endpoints of an open sweep
%   included, and valid on a non-uniform angle grid. A periodic sweep omits
%   the duplicate endpoint and wraps across period_rad. Matches
%   radia.force.coenergy_torque_from_angle_samples.

if nargin < 3 || isempty(periodic)
    periodic = false;
end
if nargin < 4 || isempty(period_rad)
    period_rad = 2*pi;
end
angles = localTable(angles_rad, "angles_rad");
coenergy = localTable(coenergy_J, "coenergy_J");
if numel(angles) ~= numel(coenergy)
    error("radia:force:SampleCount", "angles_rad and coenergy_J must have the same length");
end
if any(diff(angles) <= 0)
    error("radia:force:Angles", "angles_rad must be strictly increasing");
end
if ~isscalar(periodic)
    error("radia:force:Periodic", "periodic must be scalar logical");
end
if logical(periodic)
    if ~isscalar(period_rad) || ~isreal(period_rad) || ~isfinite(period_rad) || period_rad <= 0
        error("radia:force:Period", "period_rad must be finite and positive");
    end
    if angles(end) - angles(1) >= period_rad
        error("radia:force:Period", ...
            "periodic angles must omit the duplicate endpoint and span less than period_rad");
    end
    count = numel(angles);
    extendedAngles = [angles(end) - period_rad; angles; angles(1) + period_rad];
    extendedValues = [coenergy(end); coenergy; coenergy(1)];
    torque_Nm = localThreePointDerivative(extendedAngles, extendedValues);
    torque_Nm = torque_Nm(2:count+1);
else
    torque_Nm = localThreePointDerivative(angles, coenergy);
end
torque_Nm = reshape(torque_Nm, size(coenergy));
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
