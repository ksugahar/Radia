function derivative = relativisticLorentzRhs(time,state,chargeC,massKg, ...
    electricField,magneticField,relativistic)
%RELATIVISTICLORENTZRHS Lorentz equation with momentum as the state.
%   STATE is [x;y;z;wx;wy;wz], where w=p/(m*c)=gamma*v/c. The returned
%   derivative is [v;Q*(E+cross(v,B))/(m*c)]. The normalized momentum keeps
%   the ODE state well scaled for electrons, ions, and relativistic beams.
arguments
    time (1,1) double {mustBeFinite} %#ok<INUSD>
    state (6,1) double {mustBeFinite}
    chargeC (1,1) double {mustBeFinite,mustBeNonzero}
    massKg (1,1) double {mustBeFinite,mustBePositive}
    electricField (1,1) function_handle
    magneticField (1,1) function_handle
    relativistic (1,1) logical = true
end

speedOfLight = 299792458.0;
position = state(1:3);
normalizedMomentum = state(4:6);
if relativistic
    gamma = sqrt(1.0+dot(normalizedMomentum,normalizedMomentum));
else
    gamma = 1.0;
end
velocity = speedOfLight*normalizedMomentum/gamma;
electric = localFieldVector(electricField,position,"electric");
magnetic = localFieldVector(magneticField,position,"magnetic");
derivative = [velocity;chargeC*(electric+cross(velocity,magnetic))/ ...
    (massKg*speedOfLight)];
end

function value = localFieldVector(field,position,name)
value = double(field(position(1),position(2),position(3)));
value = value(:);
if ~isequal(size(value),[3,1]) || any(~isfinite(value))
    error("radia:tracking:Field", ...
        "%s field must return three finite components",name);
end
end
