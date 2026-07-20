function info = writePwl(file, time_s, values, options)
%WRITEPWL Write a Simulink/MATLAB waveform for an LTspice PWL source.
%   Use the generated file from a voltage or current source with
%   PWL FILE="path" in the LTspice netlist.
arguments
    file (1,1) string
    time_s (:,1) double {mustBeFinite, mustBeNonnegative}
    values (:,1) double {mustBeFinite}
    options.TimeUnit (1,1) string = "s"
end
if numel(time_s) < 2 || numel(values) ~= numel(time_s)
    error("radia:ltspice:PwlSize", ...
        "time_s and values must have equal length with at least two samples.");
end
if any(diff(time_s) <= 0)
    error("radia:ltspice:PwlTime", "time_s must be strictly increasing.");
end
folder = fileparts(file);
if strlength(folder) > 0 && ~isfolder(folder)
    mkdir(folder);
end
handle = fopen(file, 'w');
if handle < 0
    error("radia:ltspice:PwlWrite", "Could not write %s.", file);
end
cleanup = onCleanup(@() fclose(handle));
for k = 1:numel(time_s)
    fprintf(handle, '%.17g%s\t%.17g\n', ...
        time_s(k), options.TimeUnit, values(k));
end
clear cleanup
info = struct( ...
    "schema", "radia.ltspice.pwl.v1", ...
    "path", file, ...
    "sample_count", numel(time_s), ...
    "start_time_s", time_s(1), ...
    "stop_time_s", time_s(end), ...
    "minimum", min(values), ...
    "maximum", max(values));
end
