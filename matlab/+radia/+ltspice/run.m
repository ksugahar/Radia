function result = run(netlistFile, options)
%RUN Execute an LTspice netlist and return parsed waveform data.
%   RESULT = radia.ltspice.run(FILE, Parameters=STRUCT) copies FILE to a
%   run directory, replaces matching .param assignments, runs LTspice in
%   batch mode, and reads the ASCII RAW output without Simscape.
arguments
    netlistFile (1,1) string {mustBeFile}
    options.Parameters (1,1) struct = struct()
    options.Executable (1,1) string = ""
    options.OutputDirectory (1,1) string = ""
    options.Timeout_s (1,1) double {mustBePositive} = 300
    options.RawFormat (1,1) string {mustBeMember(options.RawFormat,["ascii","binary"])} = "binary"
end

executable = radia.ltspice.findExecutable(Executable=options.Executable);
sourceInput=netlistFile; [sourceFolder, stem, extension] = fileparts(sourceInput);
if ~ismember(lower(string(extension)), [".asc", ".cir", ".net", ".sp", ".spi"])
    error("radia:ltspice:NetlistRequired", ...
        "run currently accepts SPICE netlists (.cir/.net/.sp/.spi), not %s.", extension);
end
if strlength(options.OutputDirectory) == 0
    root = "C:\temp";
    if ~ispc
        root = string(tempdir);
    end
    if ~isfolder(root)
        mkdir(root);
    end
    outputDirectory = string(tempname(root));
else
    outputDirectory = options.OutputDirectory;
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end

schematic=struct.empty;
if lower(string(extension))==".asc"
    schematic=radia.ltspice.schematicToNetlist(sourceInput,Executable=executable,OutputDirectory=outputDirectory);
    netlistFile=schematic.netlist; extension=".net";
end

dependencies=radia.ltspice.collectDependencies(netlistFile);
stageDependencies(dependencies,fileparts(netlistFile),outputDirectory);

runNetlist = fullfile(outputDirectory, stem + extension);
source = string(fileread(netlistFile));
source = applyParameters(source, options.Parameters);
writeText(runNetlist, source);

rawFile = fullfile(outputDirectory, stem + ".raw");
logFile = fullfile(outputDirectory, stem + ".log");
if isfile(rawFile),delete(rawFile);end
if isfile(logFile),delete(logFile);end
started = tic;
[status, commandOutput] = executeLTspice( ...
    executable, runNetlist, rawFile, logFile, options.Timeout_s, options.RawFormat);
elapsed = toc(started);
logText = "";
if isfile(logFile)
    logText = string(fileread(logFile));
end
if status ~= 0 || ~isfile(rawFile)
    error("radia:ltspice:SimulationFailed", ...
        "LTspice failed (status %d).\n%s\n%s", status, commandOutput, logText);
end

waveform = radia.ltspice.readRaw(rawFile);
result = struct( ...
    "schema", "radia.ltspice.run.v1", ...
    "executable", executable, ...
    "source_netlist", sourceInput, ...
    "schematic_conversion", schematic, ...
    "dependencies", dependencies, ...
    "run_netlist", string(runNetlist), ...
    "raw_file", string(rawFile), ...
    "log_file", string(logFile), ...
    "output_directory", outputDirectory, ...
    "elapsed_s", elapsed, ...
    "parameters", options.Parameters, ...
    "log", logText, ...
    "waveform", waveform);
end

function stageDependencies(manifest,rootFolder,outputDirectory)
rootFolder=string(java.io.File(char(rootFolder)).getCanonicalPath());
for path=manifest.local_files(:)'
    if path==manifest.root,continue,end
    relative=extractAfter(path,strlength(rootFolder)+1); if startsWith(relative,["\","/"]),relative=extractAfter(relative,1);end
    if strlength(relative)==0||startsWith(relative,".."),continue,end
    destination=fullfile(outputDirectory,relative); folder=fileparts(destination); if ~isfolder(folder),mkdir(folder);end
    if string(java.io.File(char(path)).getCanonicalPath())~=string(java.io.File(char(destination)).getCanonicalPath()),copyfile(path,destination);end
end
end

function [status, output] = executeLTspice( ...
        executable, runNetlist, rawFile, logFile, timeout_s, rawFormat)
if ~ispc
    asciiFlag=""; if rawFormat=="ascii", asciiFlag=" -ascii"; end
    command = sprintf('"%s" "%s"%s -b -run', executable, runNetlist, asciiFlag);
    [status, output] = system(command);
    return
end

asciiFlag=""; if rawFormat=="ascii", asciiFlag=" -ascii"; end
ltArgs={'-Run','-b'}; if rawFormat=="ascii",ltArgs{end+1}='-ascii';end; ltArgs{end+1}=char(runNetlist);
quoted=cellfun(@(x)"'"+replace(string(x),"'","''")+"'",ltArgs);
script="$p=Start-Process -FilePath '"+replace(executable,"'","''")+"' -ArgumentList @("+join(quoted,",")+") -Wait -PassThru -WindowStyle Hidden; exit $p.ExitCode";
encoded=matlab.net.base64encode(unicode2native(char(script),'UTF-16LE'));
info = System.Diagnostics.ProcessStartInfo();
info.FileName = 'pwsh'; info.Arguments='-NoLogo -NoProfile -NonInteractive -EncodedCommand '+string(encoded);
info.UseShellExecute=false; info.CreateNoWindow=true;
process = System.Diagnostics.Process();
process.StartInfo = info;
if ~process.Start()
    error("radia:ltspice:ProcessStart", "Could not start LTspice.");
end
cleanup = onCleanup(@() stopOwnedProcess(process));
started = tic;
complete = false;
while toc(started) <= timeout_s
    if process.HasExited
        complete = isfile(rawFile);
        break
    end
    if isfile(rawFile) && isfile(logFile)
        logText = string(fileread(logFile));
        if contains(logText, "Total elapsed time:")
            pause(2.0);
            complete = true;
            break
        end
    end
    pause(0.05);
end
if ~complete
    if toc(started) > timeout_s
        error("radia:ltspice:Timeout", ...
            "LTspice exceeded Timeout_s=%g.", timeout_s);
    end
    status = double(process.ExitCode);
else
    status = 0;
end
output = "";
clear cleanup
stopOwnedProcess(process);
end

function stopOwnedProcess(process)
try
    if ~process.HasExited
        process.Kill();
        process.WaitForExit(5000);
    end
catch
    % The owned process may have exited between HasExited and Kill.
end
end

function text = applyParameters(text, parameters)
names = fieldnames(parameters);
for k = 1:numel(names)
    name = string(names{k});
    value = parameters.(names{k});
    if isnumeric(value) && isscalar(value) && isfinite(value)
        replacementValue = string(sprintf('%.17g', value));
    elseif isstring(value) && isscalar(value)
        replacementValue = value;
    elseif ischar(value) && isrow(value)
        replacementValue = string(value);
    else
        error("radia:ltspice:ParameterValue", ...
            "Parameter %s must be a finite scalar number or text.", name);
    end
    escaped = regexptranslate('escape', char(name));
    pattern = "(?im)(^\s*\.param\s+" + escaped + "\s*=\s*)([^\s;]+)";
    if isempty(regexp(text, pattern, 'once'))
        error("radia:ltspice:ParameterNotFound", ...
            "No .param assignment named %s exists in the netlist.", name);
    end
    text = regexprep(text, pattern, "$1" + replacementValue);
end
end

function writeText(path, text)
file = fopen(path, 'w');
if file < 0
    error("radia:ltspice:Write", "Could not write %s.", path);
end
cleanup = onCleanup(@() fclose(file));
fprintf(file, '%s', text);
clear cleanup
end
