function family = makeMagLevSmokeFamily(options)
%MAKEMAGLEVSMOKEFAMILY Create a small passive MagLev Foster test family.
%   This family exists so the packaged block can compile and simulate
%   without external artifacts. It is diagnostic data, not a validated
%   levitation design. Replace it with a shared-mode family exported by
%   radia.vim.ExportHCurlEddyFosterFamilyJSON for engineering work.

arguments
    options.SampleTime_s (1,1) double {mustBeFinite, mustBePositive} = 1.0e-4
end

positions = [0.0; 0.01];
decayRates = [160.0; 950.0];
modalPort = cat(3, [0.30; 0.12], [0.24; 0.10]);
forceOperator = zeros(3, 2, 1, 2);
forceOperator(3, :, 1, 1) = [-0.18, -0.08];
forceOperator(3, :, 1, 2) = 0.75 * [-0.18, -0.08];

family = radia.simulink.makeHCurlEddyFosterFamily(positions, decayRates, ...
    modalPort, forceOperator, SampleTime_s=options.SampleTime_s);
family.metadata = struct( ...
    "purpose", "diagnostic-smoke-only", ...
    "production_replacement", "radia.vim.ExportHCurlEddyFosterFamilyJSON");
end
