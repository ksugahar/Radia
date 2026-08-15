function torque_Nm = coenergyTorque(angles_rad, coenergy_J, periodic, period_rad)
%COENERGYTORQUE Differentiate fixed-current coenergy by mechanical angle.

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
torque_Nm = zeros(size(coenergy));
if logical(periodic)
    if ~isscalar(period_rad) || ~isreal(period_rad) || ~isfinite(period_rad) || period_rad <= 0
        error("radia:force:Period", "period_rad must be finite and positive");
    end
    count = numel(angles);
    for index = 1:count
        minusIndex = mod(index - 2, count) + 1;
        plusIndex = mod(index, count) + 1;
        angleMinus = angles(minusIndex);
        anglePlus = angles(plusIndex);
        if minusIndex > index
            angleMinus = angleMinus - period_rad;
        end
        if plusIndex < index
            anglePlus = anglePlus + period_rad;
        end
        torque_Nm(index) = (coenergy(plusIndex) - coenergy(minusIndex)) / ...
            (anglePlus - angleMinus);
    end
else
    torque_Nm(1) = (coenergy(2) - coenergy(1)) / (angles(2) - angles(1));
    torque_Nm(end) = (coenergy(end) - coenergy(end-1)) / (angles(end) - angles(end-1));
    torque_Nm(2:end-1) = (coenergy(3:end) - coenergy(1:end-2)) ./ ...
        (angles(3:end) - angles(1:end-2));
end
end

function table = localTable(value, name)
if ~isnumeric(value) || ~isreal(value) || any(~isfinite(value), "all") || ~isvector(value) || numel(value) < 3
    error("radia:force:Table", "%s must be a finite real vector with at least three samples", name);
end
table = reshape(double(value), [], 1);
end
