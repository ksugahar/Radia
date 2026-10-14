function build_mexcuda(repoRoot)
% Compatibility entry point sharing the canonical package build recipe.
addpath(fullfile(repoRoot,"packages","radia-optuna"));
buildOptunaMex;
end
