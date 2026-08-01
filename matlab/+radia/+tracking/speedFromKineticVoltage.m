function speed = speedFromKineticVoltage(chargeC,massKg,voltageV,options)
%SPEEDFROMKINETICVOLTAGE Convert accelerating voltage magnitude to speed.
arguments
    chargeC (1,1) double {mustBeFinite,mustBeNonzero}
    massKg (1,1) double {mustBeFinite,mustBePositive}
    voltageV (1,1) double {mustBeFinite,mustBeNonnegative}
    options.Relativistic (1,1) logical = true
end
energyJ = abs(chargeC)*voltageV;
if options.Relativistic
    speedOfLight = 299792458.0;
    gamma = 1.0+energyJ/(massKg*speedOfLight^2);
    speed = speedOfLight*sqrt(1.0-1.0/gamma^2);
else
    speed = sqrt(2.0*energyJ/massKg);
end
end
