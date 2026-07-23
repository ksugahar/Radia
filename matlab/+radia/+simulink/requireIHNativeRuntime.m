function paths = requireIHNativeRuntime()
%REQUIREIHNATIVERUNTIME Fail before model construction when IH MEX is absent.
names = ["radia_ih_eddy_sfun","radia_ih_thermal_sfun"];
paths = strings(size(names));
for index = 1:numel(names)
    if exist(names(index), "file") ~= 3
        error("radia:simulink:IHMexMissing", ...
            "Required native IH MEX S-Function is missing: %s.%s", ...
            names(index), mexext);
    end
    paths(index) = string(which(names(index)));
end
end
