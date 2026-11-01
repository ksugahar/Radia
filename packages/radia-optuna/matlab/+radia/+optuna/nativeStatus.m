function info = nativeStatus()
%NATIVESTATUS Report the required plain MEX/shared-kernel gateway used by radia.optuna.
info = radia.optuna.internal.NativeKernels.status();
end
