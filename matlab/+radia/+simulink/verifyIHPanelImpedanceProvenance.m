function verifyIHPanelImpedanceProvenance(configFile, sourceFile, expectedHash)
%VERIFYIHPANELIMPEDANCEPROVENANCE Require assembly to use the watched Zs bytes.
arguments
    configFile (1,1) string
    sourceFile (1,1) string
    expectedHash (1,1) string
end
config = jsondecode(fileread(configFile));
if ~isfield(config, "surface_impedance")
    error("radia:simulink:IHPanelZsProvenance", "Assembler omitted surface impedance provenance.");
end
actual = config.surface_impedance;
if strlength(sourceFile) == 0
    valid = string(actual.mode) == "uniform-linear";
else
    valid = string(actual.mode) == "specified-panel" && ...
        isfield(actual, "file_sha256") && string(actual.file_sha256) == expectedHash && ...
        isfield(actual, "source_file") && string(actual.source_file) == sourceFile;
end
if ~valid
    error("radia:simulink:IHPanelZsProvenance", ...
        "Assembled impedance does not match the watched input; rebuild from unchanged inputs.");
end
end
