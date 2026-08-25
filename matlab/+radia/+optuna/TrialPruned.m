function exception = TrialPruned(message)
%TRIALPRUNED Construct the exception used to prune an objective trial.
%   throw(radia.optuna.TrialPruned()) is the MATLAB equivalent of
%   raising optuna.TrialPruned from an Optuna objective. Study.optimize
%   catches this identifier, records a PRUNED trial, preserves the last
%   reported intermediate value, and continues with callbacks/next trials.
arguments
    message (1,1) string = "Trial was pruned."
end
exception = MException("radia:optuna:TrialPruned","%s",message);
end
