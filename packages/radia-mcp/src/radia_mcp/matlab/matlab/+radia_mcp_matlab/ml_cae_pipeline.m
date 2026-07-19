function report = ml_cae_pipeline(stage, seed)
%ML_CAE_PIPELINE Reproducible, solver-gated CAE/ML teaching pipeline.
arguments
    stage (1,1) string
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 17
end

rng(seed, "twister");
stage = lower(stage);
allowed = ["schema","pod","gp","active_learning","forward_verify", ...
    "generative_preflight","rl_preflight","all"];
if ~ismember(stage, allowed)
    error("radia_mcp_matlab:UnknownMlCaeStage", ...
        "stage must be one of: %s", strjoin(allowed, ", "));
end

report = struct("schema", "radia-mcp.matlab.ml-cae-pipeline.v1", ...
    "stage", stage, "seed", seed, "units", struct("x", "1", "qoi", "Pa"), ...
    "provenance", struct("source", "analytic_forward_model", ...
    "generated_candidate_is_ground_truth", false), "checks", struct());

% A cheap analytic forward model keeps the test deterministic. Production
% callers replace this with a FEM/BEM/NGSolve/attached LiveLink callback.
x = linspace(0, 1, 32)';
snapshots = [sin(pi*x), sin(2*pi*x), 0.5*cos(pi*x)] + 0.02*x*(1:3);
qoi = forwardModel(x);

[u,s,~] = svd(snapshots - mean(snapshots, 1), "econ");
energy = cumsum(diag(s).^2) / sum(diag(s).^2);
rank99 = find(energy >= 0.99, 1, "first");
pod = struct("rank99", rank99, "energy99", energy(rank99), ...
    "orthogonality_error", norm(u(:,1:rank99)'*u(:,1:rank99)-eye(rank99), "fro"));

trainIdx = (1:4:numel(x))';
xt = x(trainIdx); yt = qoi(trainIdx);
[yp, ys] = explicitGpPredict(xt, yt, x);
gpBackend = "explicit_cholesky";
relRmse = sqrt(mean((yp-qoi).^2)) / max(rms(qoi), eps);
[~, nextIdx] = max(ys);
active = struct("next_x", x(nextIdx), "max_predictive_sigma", ys(nextIdx));

tol = 0.08;
verification = struct("relative_rmse", relRmse, "tolerance", tol, ...
    "forward_solver", "analytic_forward_model", "passed", relRmse <= tol, ...
    "candidate_promoted", relRmse <= tol);

versionInfo = ver;
installed = string({versionInfo.Name});
toolboxes = struct( ...
    "statistics_machine_learning", any(contains(installed, "Statistics and Machine Learning")), ...
    "deep_learning", any(contains(installed, "Deep Learning")), ...
    "optimization", any(contains(installed, "Optimization")), ...
    "reinforcement_learning", any(contains(installed, "Reinforcement Learning")));

report.toolboxes = toolboxes;
report.experiment = struct("split", "deterministic_stride_4", ...
    "num_samples", numel(x), "num_train", numel(trainIdx));
report.pod = pod;
report.gp = struct("relative_rmse", relRmse, ...
    "mean_predictive_sigma", mean(ys), "kernel", "squaredexponential", ...
    "backend", gpBackend);
report.active_learning = active;
report.forward_verification = verification;
report.generative_preflight = struct("ready", toolboxes.deep_learning, ...
    "requires_forward_verification", true, "required_contracts", ...
    ["seed","units","provenance","latent_schema","solver_callback"]);
report.rl_preflight = struct("ready", toolboxes.reinforcement_learning, ...
    "requires_forward_verification", true, "required_contracts", ...
    ["observation_schema","action_bounds","reward_units","termination","solver_callback"]);

report.checks.schema_complete = seed >= 0 && report.units.qoi ~= "";
report.checks.pod_orthogonal = pod.orthogonality_error < 1e-10;
report.checks.gp_uncertainty_finite = all(isfinite(ys)) && all(ys >= 0);
report.checks.active_learning_in_bounds = active.next_x >= 0 && active.next_x <= 1;
report.checks.forward_verified = verification.passed;

switch stage
    case "schema"
        required = report.checks.schema_complete;
    case "pod"
        required = report.checks.schema_complete && report.checks.pod_orthogonal;
    case "gp"
        required = report.checks.gp_uncertainty_finite;
    case "active_learning"
        required = report.checks.active_learning_in_bounds;
    case "forward_verify"
        required = report.checks.forward_verified;
    case "generative_preflight"
        required = report.generative_preflight.ready;
    case "rl_preflight"
        required = report.rl_preflight.ready;
    otherwise
        required = all(structfun(@(v) logical(v), report.checks));
end
report.ok = required;
end

function y = forwardModel(x)
y = 1.2 + sin(2*pi*x) + 0.15*cos(6*pi*x) + 0.1*x;
end

function [mu, sigma] = explicitGpPredict(xTrain, yTrain, xQuery)
% Small transparent GP fallback for sessions without fitrgp on the path.
ell = 0.18;
noise = 1e-8;
ktt = exp(-0.5*((xTrain-xTrain')./ell).^2) + noise*eye(numel(xTrain));
ktq = exp(-0.5*((xTrain-xQuery')./ell).^2);
l = chol(ktt, "lower");
alpha = l' \ (l \ yTrain);
mu = ktq' * alpha;
v = l \ ktq;
variance = max(1 - sum(v.^2, 1)', 0);
sigma = sqrt(variance);
end
