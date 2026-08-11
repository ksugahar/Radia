function state = cartesianState(positionM,kineticMomentumKgMS,options)
%CARTESIANSTATE Create an explicit SI Cartesian beam state struct.
arguments
    positionM (1,3) double {mustBeReal,mustBeFinite}
    kineticMomentumKgMS (1,3) double {mustBeReal,mustBeFinite}
    options.TimeS (1,1) double {mustBeReal,mustBeFinite} = 0
    options.PathLengthM (1,1) double {mustBeReal,mustBeFinite} = 0
end
state = struct( ...
    position_m=double(positionM), ...
    kinetic_momentum_kg_m_s=double(kineticMomentumKgMS), ...
    time_s=double(options.TimeS), ...
    path_length_m=double(options.PathLengthM));
end
