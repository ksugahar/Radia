function value = radia_optuna_mesh_failure(trial)
%RADIA_OPTUNA_MESH_FAILURE Deterministic CAE failure fixture.
trial.suggestFloat("mesh_scale", 0.5, 1.5);
error("radia:solver:meshInvalid", ...
    "Mesh Jacobian is invalid for this frozen negative control.");
value = NaN;
end
