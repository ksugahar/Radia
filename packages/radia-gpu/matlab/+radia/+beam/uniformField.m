function field = uniformField(magneticT,options)
%UNIFORMFIELD Create a uniform E/B field description in SI units.
arguments
    magneticT (1,3) double {mustBeReal,mustBeFinite}
    options.ElectricVM (1,3) double {mustBeReal,mustBeFinite} = [0 0 0]
end
field = struct(type='uniform',magnetic_t=double(magneticT), ...
    electric_v_m=double(options.ElectricVM));
end
