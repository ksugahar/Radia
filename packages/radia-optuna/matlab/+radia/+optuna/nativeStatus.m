function info = nativeStatus()
%NATIVESTATUS Report optional MEX acceleration used by radia.optuna.
info = radia.optuna.internal.NativeKernels.status();
end
