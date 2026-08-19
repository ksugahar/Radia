function result = acceleratorSectionOptics(functionName, positional, options)
%ACCELERATORSECTIONOPTICS Call section-optics inverse design operations.
%   Explicit batch boundary over radia.accelerator_section_optics: an
%   optics specification written on one section's transfer block becomes
%   the field difference that delivers it, resting on the exact ordered
%   composition M = M_after . M_S . M_before. Provides the
%   representation-independent whitened minimum-norm step under a
%   physical field metric, the inert-knob filter, the horizontal/vertical
%   coupling diagnostic, and both the multipole and CanonicalHCurl chain
%   representations. Intended for design steps and artifact generation,
%   not Simulink step-time use.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython( ...
    "radia.accelerator_section_optics",functionName,positional, ...
    Keywords=options.Keywords);
end
