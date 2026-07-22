function target = openIH(options)
%OPENIH Open the production IH block or a configured IH control model.
%   radia.simulink.openIH() opens Applications/Induction Heating in the
%   installed Radia Simulink library.
%
%   radia.simulink.openIH(Plant=PLANT,EddyLUT=LUT) builds and opens the
%   separated Eddy Current / Thermal model with its IH Parameters block.
%   Current, rotation angle, and ambient temperature remain external
%   Simulink signals.

arguments
    options.ModelName (1,1) string = "radia_ih"
    options.Plant (1,1) struct = struct()
    options.EddyLUT (1,1) struct = struct()
    options.PlantBlock (1,1) string = "standard"
    options.OutputDirectory (1,1) string = ""
    options.Save (1,1) logical = true
    options.Open (1,1) logical = true
end

hasPlant = ~isempty(fieldnames(options.Plant));
hasEddyLUT = ~isempty(fieldnames(options.EddyLUT));
if xor(hasPlant, hasEddyLUT)
    error("radia:simulink:IHLaunchInputs", ...
        "Plant and EddyLUT must be supplied together.");
end

if hasPlant
    if bdIsLoaded(char(options.ModelName))
        target = options.ModelName;
        if options.Open
            open_system(target);
        end
        return
    end
    target = radia.simulink.buildIHControlModel( ...
        options.ModelName, options.Plant, options.EddyLUT, ...
        PlantBlock=options.PlantBlock, ...
        OutputDirectory=options.OutputDirectory, ...
        Save=options.Save, Open=options.Open);
    if ~options.Save
        target = options.ModelName;
    end
    return
end

libraryFile = string(which("radia_simulink_library.slx"));
if strlength(libraryFile) == 0
    if ispc
        outputDirectory = fullfile("C:\temp", "radia_simulink_library");
    else
        outputDirectory = fullfile(tempdir, "radia_simulink_library");
    end
    libraryFile = radia.simulink.buildLibrary( ...
        OutputDirectory=string(outputDirectory));
end
load_system(libraryFile);
target = "radia_simulink_library/Applications/Induction Heating";
if options.Open
    open_system(target);
end
end
