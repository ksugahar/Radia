function reference = referenceParticle(species,kineticEnergyEV)
%REFERENCEPARTICLE Construct native relativistic reference-particle data.
arguments
    species (1,1) struct
    kineticEnergyEV (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
end
reference = radia.internal.callMex( ...
    'beam.reference_particle.from_kinetic_energy_ev', ...
    species,double(kineticEnergyEV));
end
