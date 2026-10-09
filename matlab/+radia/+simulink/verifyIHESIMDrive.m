function verifyIHESIMDrive(config, current_A)
%VERIFYIHESIMDRIVE Reject use outside a frozen ESIM constitutive operating band.
arguments
    config (1,1) struct
    current_A (1,1) double
end
if ~isfield(config,"surface_impedance") || ...
        string(config.surface_impedance.mode) ~= "per-panel-esim", return; end
state = config.surface_impedance;
reference = double(state.reference_current_A);
band = double(state.drive_relative_band);
if ~(isfinite(reference) && reference > 0 && isfinite(band) && band >= 0 && band <= .1)
    error("radia:simulink:IHESIMDriveContract", "Invalid frozen ESIM drive contract.");
end
if ~isfinite(current_A) || abs(abs(current_A)/reference-1) > band+32*eps
    error("radia:simulink:IHESIMDriveBand", ...
        "Peak current amplitude (envelope) is outside frozen ESIM reference %.9g A (relative band %.3g). Instantaneous carrier samples and zero-start/soft-start ramps are unsupported. Rebuild at a new reference current.", reference, band);
end
end
