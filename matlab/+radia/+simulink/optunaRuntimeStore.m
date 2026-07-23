function value = optunaRuntimeStore(action, key, runtime)
%OPTUNARUNTIMESTORE Keep handle-valued study state outside numeric DWork.
persistent runtimes
if isempty(runtimes)
    runtimes = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
key = char(string(key));
switch string(action)
    case "set"
        runtimes(key) = runtime;
        value = runtime;
    case "get"
        if isKey(runtimes, key)
            value = runtimes(key);
        else
            value = [];
        end
    case "remove"
        if isKey(runtimes, key)
            remove(runtimes, key);
        end
        value = [];
    otherwise
        error('radia:simulink:OptunaRuntimeAction', ...
            'Unknown Optuna runtime-store action: %s', action);
end
end
