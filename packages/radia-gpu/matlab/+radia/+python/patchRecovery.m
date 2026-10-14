function result = patchRecovery(functionName, positional, options)
%PATCHRECOVERY Optional batch gradient recovery; does not modify the FE solve.
% Use recover_vertex_patches or recover_mixed_omega_p1 in radia.patch_recovery.
% NGSolve objects remain in result.value. Not a per-step Simulink backend.
arguments
    functionName (1,1) string
    positional (1,:) cell = {}
    options.Keywords (1,1) struct = struct()
end
result = radia.internal.callPython("radia.patch_recovery", ...
    functionName, positional, Keywords=options.Keywords);
end
