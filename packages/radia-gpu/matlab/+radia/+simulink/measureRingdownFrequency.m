function [frequency_Hz, details] = measureRingdownFrequency(signal, options)
%MEASURERINGDOWNFREQUENCY Estimate frequency from interpolated zero crossings.
arguments
    signal
    options.IgnoreBefore_s (1,1) double {mustBeNonnegative,mustBeFinite} = 0
    options.MinimumCrossings (1,1) double ...
        {mustBeInteger,mustBeGreaterThanOrEqual(options.MinimumCrossings,4)} = 8
end

[time_s, values] = unpackSignal(signal);
valid = isfinite(time_s) & isfinite(values) & time_s >= options.IgnoreBefore_s;
time_s = time_s(valid);
values = values(valid);
if numel(time_s) < 3 || any(diff(time_s) <= 0)
    error("radia:simulink:RingdownTime", ...
        "Ring-down samples require a finite, strictly increasing time vector.");
end

indices = find((values(1:end-1) > 0 & values(2:end) <= 0) | ...
    (values(1:end-1) < 0 & values(2:end) >= 0));
if numel(indices) < options.MinimumCrossings
    error("radia:simulink:RingdownCrossings", ...
        "Need at least %d zero crossings; found %d.", ...
        options.MinimumCrossings, numel(indices));
end

t0 = time_s(indices);
t1 = time_s(indices + 1);
y0 = values(indices);
y1 = values(indices + 1);
crossingTimes_s = t0 - y0 .* (t1 - t0) ./ (y1 - y0);
halfPeriods_s = diff(crossingTimes_s);
medianHalfPeriod_s = median(halfPeriods_s);
frequency_Hz = 1 / (2 * medianHalfPeriod_s);
absoluteDeviation_s = median(abs(halfPeriods_s - medianHalfPeriod_s));

details = struct( ...
    "crossing_count", numel(crossingTimes_s), ...
    "median_half_period_s", medianHalfPeriod_s, ...
    "relative_half_period_mad", absoluteDeviation_s / medianHalfPeriod_s, ...
    "first_crossing_s", crossingTimes_s(1), ...
    "last_crossing_s", crossingTimes_s(end));
end

function [time_s, values] = unpackSignal(signal)
if isa(signal, "timeseries")
    time_s = double(signal.Time(:));
    values = double(signal.Data(:));
elseif isnumeric(signal) && ismatrix(signal) && size(signal,2) == 2
    time_s = double(signal(:,1));
    values = double(signal(:,2));
else
    error("radia:simulink:RingdownSignal", ...
        "signal must be a timeseries or a numeric [time, value] matrix.");
end
if numel(time_s) ~= numel(values)
    error("radia:simulink:RingdownSignal", ...
        "Ring-down time and value arrays must have equal lengths.");
end
end
