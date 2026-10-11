function result = runLTspiceTrial(netlistFile, parameters, executable, runDirectory)
%RUNLTSPICETRIAL Worker-safe LTspice execution entry point.
result = radia.ltspice.run(netlistFile, Parameters=parameters, ...
    Executable=executable, OutputDirectory=runDirectory);
end
