function [frequency_Hz, details] = measureSeriesResonance(source, options)
%MEASURESERIESRESONANCE Measure AC resonance from Im(V/I) = 0.
arguments
    source
    options.VoltageTrace (1,1) string = "V(in)"
    options.CurrentTrace (1,1) string = "I(Vdrive)"
end

raw = resolveRaw(source);
if ~strcmpi(string(raw.Data.analysis), "ac")
    error("radia:ltspice:ACRequired", ...
        "Series resonance measurement requires an LTspice AC analysis.");
end
frequency = real(raw.getTrace("frequency"));
voltage = raw.getTrace(options.VoltageTrace);
sourceCurrent = raw.getTrace(options.CurrentTrace);
valid = isfinite(frequency) & frequency > 0 & ...
    isfinite(voltage) & isfinite(sourceCurrent) & sourceCurrent ~= 0;
frequency = frequency(valid);
voltage = voltage(valid);
sourceCurrent = sourceCurrent(valid);
if numel(frequency) < 3 || any(diff(frequency) <= 0)
    error("radia:ltspice:ACFrequencyAxis", ...
        "AC frequency samples must be finite and strictly increasing.");
end

% LTspice reports voltage-source current into the positive source terminal.
impedance = -voltage ./ sourceCurrent;
if median(real(impedance)) < 0
    impedance = -impedance;
end
reactance = imag(impedance);
indices = find((reactance(1:end-1) <= 0 & reactance(2:end) >= 0) | ...
    (reactance(1:end-1) >= 0 & reactance(2:end) <= 0));
if isempty(indices)
    error("radia:ltspice:ResonanceNotBracketed", ...
        "The AC sweep does not bracket an input-reactance zero crossing.");
end
[~, choice] = min(max(abs([reactance(indices), reactance(indices+1)]), [], 2));
index = indices(choice);
f0 = frequency(index);
f1 = frequency(index+1);
x0 = reactance(index);
x1 = reactance(index+1);
if x1 == x0
    frequency_Hz = 0.5 * (f0 + f1);
else
    frequency_Hz = f0 - x0 * (f1-f0) / (x1-x0);
end

currentMagnitude = abs(sourceCurrent);
[~, peakIndex] = max(currentMagnitude);
peakFrequency_Hz = frequency(peakIndex);
if peakIndex > 1 && peakIndex < numel(frequency)
    fitFrequency = frequency(peakIndex-1:peakIndex+1);
    fitMagnitude = log(currentMagnitude(peakIndex-1:peakIndex+1));
    coefficients = polyfit(fitFrequency, fitMagnitude, 2);
    if coefficients(1) < 0
        candidate = -coefficients(2) / (2*coefficients(1));
        if candidate >= fitFrequency(1) && candidate <= fitFrequency(end)
            peakFrequency_Hz = candidate;
        end
    end
end

weight = (frequency_Hz-f0) / (f1-f0);
resistanceAtResonance = real(impedance(index)) + weight * ...
    (real(impedance(index+1))-real(impedance(index)));
details = struct( ...
    "analysis", "ac", ...
    "method", "linear interpolation of input-reactance zero", ...
    "point_count", numel(frequency), ...
    "bracket_frequency_Hz", [f0 f1], ...
    "bracket_reactance_Ohm", [x0 x1], ...
    "input_resistance_at_resonance_Ohm", resistanceAtResonance, ...
    "current_peak_frequency_Hz", peakFrequency_Hz, ...
    "zero_vs_peak_relative_difference", ...
        abs(frequency_Hz-peakFrequency_Hz)/frequency_Hz, ...
    "peak_current_A", max(currentMagnitude), ...
    "voltage_trace", options.VoltageTrace, ...
    "current_trace", options.CurrentTrace);
end

function raw = resolveRaw(source)
if isa(source, "radia.ltspice.RawRead")
    raw = source;
elseif isstruct(source) && isscalar(source) && isfield(source, "raw_file")
    raw = radia.ltspice.RawRead(string(source.raw_file));
elseif (isstring(source) && isscalar(source)) || ...
        (ischar(source) && isrow(source))
    raw = radia.ltspice.RawRead(string(source));
else
    error("radia:ltspice:RawSource", ...
        "source must be a run result, RawRead object, or RAW-file path.");
end
end
