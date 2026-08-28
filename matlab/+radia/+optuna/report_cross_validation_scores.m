function report_cross_validation_scores(trial,scores)
%REPORT_CROSS_VALIDATION_SCORES Attach fold scores to a running trial.

scores=reshape(double(scores),1,[]);
if numel(scores)<=1
    error("radia:optuna:CrossValidationScores", ...
        "The length of scores must be greater than one.");
end
trial.set_system_attr("terminator:cv_scores",scores);
end
