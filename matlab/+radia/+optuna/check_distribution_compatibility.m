function compatible=check_distribution_compatibility(left,right)
%CHECK_DISTRIBUTION_COMPATIBILITY Apply Optuna 4.9 dynamic-space rules.
arguments
    left (1,1) struct
    right (1,1) struct
end
if left.kind~=right.kind
    error("radia:optuna:IncompatibleDistribution", ...
        "Cannot set different distribution kind to the same parameter name.");
end
if ismember(left.kind,["float","integer"]) && left.log~=right.log
    error("radia:optuna:IncompatibleDistribution", ...
        "Cannot set different log configuration to the same parameter name.");
end
if left.kind=="categorical" && ...
        ~isequal(radia.optuna.internal.DistributionCodec.choiceTokens( ...
        left.choices),radia.optuna.internal.DistributionCodec.choiceTokens( ...
        right.choices))
    error("radia:optuna:IncompatibleDistribution", ...
        "CategoricalDistribution does not support dynamic value space.");
end
compatible=true;
end
