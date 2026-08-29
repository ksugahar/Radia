function sampler = samplerFromName(name, seed, options)
%SAMPLERFROMNAME Construct a sampler from its short name.
%   SAMPLER = SAMPLERFROMNAME(NAME, SEED) returns the sampler the short
%   name selects, seeded with SEED. An empty SEED leaves the sampler's own
%   default, which upstream defines as fresh private entropy.
%
%   This is the single place the short names map onto sampler objects. The
%   mapping used to live inside optunaSFunction as a local function, where
%   nothing outside Simulink could reach it and its seed was pinned to 0;
%   both the Simulink block and radia.optuna.optimize now call this.
%
%   "auto" is rejected here on purpose. Choosing a sampler needs the study's
%   objective count, trial budget and search-space shape, which this function
%   does not have; radia.optuna.internal.AutoSamplerPolicy.choose owns that
%   decision, and the caller resolves "auto" before calling.
%
%   See also radia.optuna.internal.AutoSamplerPolicy, radia.optuna.optimize.

arguments
    name (1,1) string
    seed = []
    options.NStartupTrials double = []
    options.PopulationSize double = []
end

positionalSeed = {};
namedSeed = {};
if ~isempty(seed)
    positionalSeed = {seed};
    namedSeed = {"Seed", seed};
end

startupTrials = options.NStartupTrials;
populationSize = options.PopulationSize;

switch name
    case "random"
        sampler = radia.optuna.RandomSampler(positionalSeed{:});
    case "tpe"
        sampler = radia.optuna.TPESampler(namedSeed{:}, ...
            "NStartupTrials", defaulted(startupTrials, 10));
    case "cmaes"
        sampler = radia.optuna.CmaEsSampler(namedSeed{:}, ...
            "NStartupTrials", defaulted(startupTrials, 1));
    case "gp"
        sampler = radia.optuna.GPSampler(namedSeed{:}, ...
            "NStartupTrials", defaulted(startupTrials, 10), ...
            "DeterministicObjective", true);
    case "motpe"
        sampler = radia.optuna.MOTPESampler(namedSeed{:}, ...
            "NStartupTrials", defaulted(startupTrials, 20));
    case "nsgaii"
        sampler = radia.optuna.NSGAIISampler(namedSeed{:}, ...
            "PopulationSize", defaulted(populationSize, 24));
    case "nsgaiii"
        sampler = radia.optuna.NSGAIIISampler(namedSeed{:}, ...
            "PopulationSize", defaulted(populationSize, 24));
    case "bruteforce"
        sampler = radia.optuna.BruteForceSampler(namedSeed{:});
    case "qmc"
        sampler = radia.optuna.QMCSampler(namedSeed{:}, ...
            "QMCType", "sobol", "Scramble", true);
    case "auto"
        error("radia:optuna:SamplerName", ...
            "Resolve 'auto' through " + ...
            "radia.optuna.internal.AutoSamplerPolicy.choose before " + ...
            "calling samplerFromName; the choice needs the study's " + ...
            "objective count, trial budget and search-space shape.");
    otherwise
        error("radia:optuna:SamplerName", ...
            "Unknown sampler '%s'. Available: random, tpe, cmaes, gp, " + ...
            "motpe, nsgaii, nsgaiii, bruteforce, qmc.", name);
end
end

function value = defaulted(candidate, fallback)
if isempty(candidate)
    value = fallback;
else
    value = candidate;
end
end
