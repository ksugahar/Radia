function target = openElectromagnet(options)
%OPENELECTROMAGNET Open the standalone electromagnet topology model.
arguments
    options.ModelName (1,1) string = "radia_electromagnet"
    options.OutputDirectory (1,1) string = ""
    options.Open (1,1) logical = true
end

modelFile = "";
if options.ModelName == "radia_electromagnet"
    modelFile = string(which("radia_electromagnet.slx"));
end
if strlength(modelFile) == 0
    outputDirectory = options.OutputDirectory;
    if strlength(outputDirectory) == 0
        outputDirectory = fullfile("C:\temp","radia_electromagnet");
    end
    modelFile = radia.simulink.buildElectromagnetOptimizationModel( ...
        ModelName=options.ModelName,OutputDirectory=outputDirectory, ...
        Open=false);
end
load_system(modelFile);
target = options.ModelName;
if options.Open
    open_system(target);
end
end
