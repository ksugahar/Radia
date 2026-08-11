function config = trackingConfig(species,state,field,independent)
%TRACKINGCONFIG Checked common input for native beam tracking commands.
arguments
    species (1,1) struct
    state (1,1) struct
    field (1,1) struct
    independent (1,1) string {mustBeMember(independent, ...
        ["time","path_length","azimuth"])}
end
requiredSpecies = ["charge_c","rest_mass_kg"];
requiredState = ["position_m","kinetic_momentum_kg_m_s"];
if ~all(isfield(species,requiredSpecies))
    error("radia:beam:InvalidSpecies", ...
        "Species must define charge_c and rest_mass_kg.");
end
if ~all(isfield(state,requiredState))
    error("radia:beam:InvalidState", ...
        "State must define position_m and kinetic_momentum_kg_m_s.");
end
config = struct(schema='radia.beam.tracking.v1',species=species, ...
    state=state,field=field,independent=char(independent));
end
