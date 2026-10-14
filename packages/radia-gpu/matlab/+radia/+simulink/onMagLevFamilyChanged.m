function onMagLevFamilyChanged(blockPath)
%ONMAGLEVFAMILYCHANGED Apply the standalone MagLev parameter mask.

arguments
    blockPath (1,1) string
end

modelName = string(bdroot(blockPath));
familyFile = string(get_param(blockPath, "family_file"));
interpolation = string(get_param(blockPath, "interpolation"));
extrapolation = string(get_param(blockPath, "extrapolation"));
sampleTime = str2double(string(get_param(blockPath, "sample_time_s")));
if ~isfinite(sampleTime) || sampleTime <= 0
    error("radia:simulink:MagLevSampleTime", ...
        "Sample time must be a positive finite scalar.");
end

workspace = get_param(modelName, "ModelWorkspace");
if strlength(strtrim(familyFile)) > 0
    family = radia.simulink.loadHCurlEddyFosterFamily( ...
        familyFile, SampleTime_s=sampleTime, ...
        Interpolation=interpolation, Extrapolation=extrapolation);
else
    try
        family = workspace.getVariable("radia_maglev_family");
    catch
        family = radia.simulink.makeMagLevSmokeFamily( ...
            SampleTime_s=sampleTime);
    end
    family = resampleFamily(family, sampleTime);
    family.interpolation = interpolation;
    family.extrapolation = extrapolation;
end

workspace.assignin("radia_maglev_family", family);
set_param(modelName, "FixedStep", char(compose("%.17g", sampleTime)));
set_param(modelName, "SimulationCommand", "update");
end

function family = resampleFamily(family, sampleTime)
% Only the zero-order hold depends on the sample time; the shared Foster
% modes and the height-indexed operators are unchanged.
updated = radia.simulink.makeHCurlEddyFosterFamily(family.positions_m, ...
    family.decay_rates, family.modal_port_rhs, family.force_operator, ...
    SampleTime_s=sampleTime, Interpolation=family.interpolation, ...
    Extrapolation=family.extrapolation, InitialState=family.x0);
updated.source_file = family.source_file;
updated.metadata = family.metadata;
family = updated;
end
