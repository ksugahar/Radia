function result = analyzeRLCResonance(resistance_Ohm, inductance_H, capacitance_F)
%ANALYZERLCRESONANCE Return the characteristic frequencies of a series RLC.
arguments
    resistance_Ohm (1,1) double {mustBeNonnegative,mustBeFinite}
    inductance_H (1,1) double {mustBePositive,mustBeFinite}
    capacitance_F (1,1) double {mustBePositive,mustBeFinite}
end

omega0 = 1 / sqrt(inductance_H * capacitance_F);
alpha = resistance_Ohm / (2 * inductance_H);
omegaRingSquared = omega0^2 - alpha^2;
if omegaRingSquared <= 0
    error("radia:simulink:RLCOverdamped", ...
        "The selected RLC circuit is not underdamped.");
end

omegaCapacitorPeakSquared = omega0^2 - 2 * alpha^2;
capacitorPeak = NaN;
if omegaCapacitorPeakSquared > 0
    capacitorPeak = sqrt(omegaCapacitorPeakSquared) / (2*pi);
end

result = struct( ...
    "natural_frequency_Hz", omega0 / (2*pi), ...
    "ringdown_frequency_Hz", sqrt(omegaRingSquared) / (2*pi), ...
    "series_current_peak_frequency_Hz", omega0 / (2*pi), ...
    "capacitor_voltage_peak_frequency_Hz", capacitorPeak, ...
    "decay_rate_per_s", alpha, ...
    "quality_factor", omega0 / (2*alpha));
if resistance_Ohm == 0
    result.quality_factor = Inf;
end
end
