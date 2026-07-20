function result = runLTspice(netlistFile, options)
%RUNLTSPICE Run an LTspice circuit from a Simulink-oriented workflow.
%   RESULT = radia.simulink.runLTspice(FILE, InputSignals=SIGNALS) writes
%   each logged Simulink signal to an LTspice PWL file before running the
%   netlist. SIGNALS is a scalar struct whose fields contain Nx2 [time,value]
%   arrays. The netlist can reference <field>.pwl with PWL FILE="...".
arguments
    netlistFile (1,1) string {mustBeFile}
    options.InputSignals (1,1) struct = struct()
    options.Parameters (1,1) struct = struct()
    options.Executable (1,1) string = ""
    options.OutputDirectory (1,1) string = ""
    options.Timeout_s (1,1) double {mustBePositive} = 300
end

outputDirectory = options.OutputDirectory;
if strlength(outputDirectory) == 0
    outputDirectory = string(tempname("C:\temp"));
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end

names = fieldnames(options.InputSignals);
pwl = struct();
for k = 1:numel(names)
    samples = options.InputSignals.(names{k});
    if ~isnumeric(samples) || size(samples, 2) ~= 2
        error("radia:simulink:LTspiceInput", ...
            "InputSignals.%s must be an Nx2 [time_s,value] array.", names{k});
    end
    if size(samples,1)==1
        samples=[samples; samples(1,1)+eps(max(1,samples(1,1))), samples(1,2)];
    end
    destination = fullfile(outputDirectory, string(names{k}) + ".pwl");
    pwl.(names{k}) = radia.ltspice.writePwl( ...
        destination, samples(:, 1), samples(:, 2));
end

sourceNetlist = string(fileread(netlistFile));
[~, stem, extension] = fileparts(netlistFile);
stagedNetlist = fullfile(outputDirectory, stem + extension);
for k = 1:numel(names)
    sourceNetlist = replace(sourceNetlist, ...
        "${" + string(names{k}) + "_PWL}", ...
        replace(string(pwl.(names{k}).path), "\", "/"));
end
writeText(stagedNetlist, sourceNetlist);
result = radia.ltspice.run(stagedNetlist, ...
    Parameters=options.Parameters, Executable=options.Executable, ...
    OutputDirectory=outputDirectory, Timeout_s=options.Timeout_s);
result.schema = "radia.simulink.ltspice.run.v1";
result.input_pwl = pwl;
end

function writeText(path, text)
handle = fopen(path, 'w');
if handle < 0
    error("radia:simulink:LTspiceWrite", "Could not write %s.", path);
end
cleanup = onCleanup(@() fclose(handle));
fprintf(handle, '%s', text);
clear cleanup
end
