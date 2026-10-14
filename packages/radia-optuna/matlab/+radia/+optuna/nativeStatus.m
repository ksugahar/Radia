function info = nativeStatus()
%NATIVESTATUS Report the required mexcuda gateway used by radia.optuna.
info = radia.optuna.internal.NativeKernels.status();
end
