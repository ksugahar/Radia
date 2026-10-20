function paths = requireIHNativeRuntime()
%REQUIREIHNATIVERUNTIME Require readable wrappers and the shared native ABI.
names = ["radia_mex","radia_ih_eddy_sfun","radia_ih_thermal_sfun", ...
    "radia_ih_monitor_sfun"];
paths = strings(size(names));
for index = 1:numel(names)
    expectedKind = 2;
    if names(index) == "radia_mex", expectedKind = 3; end
    if exist(names(index), "file") ~= expectedKind
        error("radia:simulink:IHRuntimeMissing", ...
            "Required IH runtime entry point is missing: %s", names(index));
    end
    paths(index) = string(which(names(index)));
end
end
