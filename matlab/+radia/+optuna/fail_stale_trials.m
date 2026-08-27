function fail_stale_trials(study)
%FAIL_STALE_TRIALS Fail trials reported stale by a heartbeat backend.
arguments
    study
end
warning("radia:optuna:ExperimentalWarning", ...
    "fail_stale_trials is experimental in Optuna 4.9.0.");
if isa(study,"radia.optuna.Study")
    % Native table/MAT studies do not implement heartbeat, matching
    % upstream's no-op for storage backends outside BaseHeartbeat.
    return
end
try
    module=py.importlib.import_module("optuna.storages");
    operation=py.builtins.getattr(module,"fail_stale_trials");
    operation(study);
catch cause
    error("radia:optuna:HeartbeatStudy", ...
        "study must be a MATLAB Study or upstream optuna.Study: %s", ...
        cause.message);
end
end
