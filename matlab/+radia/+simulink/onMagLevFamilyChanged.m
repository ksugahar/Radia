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
    family = radia.simulink.loadHCurlEddyCLNFamily( ...
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
models = family.models;
for index = 1:numel(models)
    old = models{index};
    updated = radia.simulink.makeHCurlEddyCLNModel( ...
        old.resistance, old.inductance, old.port_rhs, ...
        SampleTime_s=sampleTime, InitialState=old.x0);
    preserved = ["height_m", "exchange_schema", "force_operator", "metadata"];
    for fieldName = preserved
        if isfield(old, fieldName)
            updated.(fieldName) = old.(fieldName);
        end
    end
    models{index} = updated;
end
family.models = models;
family.sample_time_s = sampleTime;
end
