function result = solveSeriesResonance(coil, capacitance_F, options)
%SOLVESERIESRESONANCE Find the self-consistent PEEC series resonance.
%
% The root is defined by
%   imag(Z_coil(f) + 1/(j*2*pi*f*C)) = 0.
% Therefore L_eff is evaluated at the resonance being sought, rather than
% frozen at DC or at an unrelated frequency.
arguments
    coil (1,1) struct
    capacitance_F (1,1) double {mustBePositive,mustBeFinite}
    options.DCProperties (1,1) struct = struct()
    options.FrequencyBounds_Hz double = zeros(0,2)
    options.RelativeFrequencyTolerance (1,1) double ...
        {mustBePositive,mustBeFinite} = 1e-10
    options.RelativeReactanceTolerance (1,1) double ...
        {mustBePositive,mustBeFinite} = 1e-10
    options.MaxIterations (1,1) double ...
        {mustBeInteger,mustBePositive} = 64
    options.MaxBracketExpansions (1,1) double ...
        {mustBeInteger,mustBeNonnegative} = 20
end

if isempty(fieldnames(options.DCProperties))
    dc = radia.peec.seriesCoilProperties(coil);
else
    dc = options.DCProperties;
end
if ~all(isfield(dc,["inductance_H","resistance_Ohm"]))
    error("radia:peec:DCProperties", ...
        "DCProperties must contain inductance_H and resistance_Ohm.");
end
dcEstimate = 1/(2*pi*sqrt(double(dc.inductance_H)*capacitance_F));
if isempty(options.FrequencyBounds_Hz)
    low = dcEstimate/2;
    high = dcEstimate*2;
else
    if numel(options.FrequencyBounds_Hz) ~= 2 || ...
            any(~isfinite(options.FrequencyBounds_Hz)) || ...
            any(options.FrequencyBounds_Hz <= 0)
        error("radia:peec:ResonanceBounds", ...
            "FrequencyBounds_Hz must contain two finite positive values.");
    end
    bounds = sort(double(options.FrequencyBounds_Hz(:)));
    low = bounds(1);
    high = bounds(2);
    if low == high
        error("radia:peec:ResonanceBounds", ...
            "FrequencyBounds_Hz must be strictly increasing.");
    end
end

historyFrequency = zeros(options.MaxIterations+2* ...
    options.MaxBracketExpansions+4,1);
historyResidual = zeros(size(historyFrequency));
historyResistance = zeros(size(historyFrequency));
historyInductance = zeros(size(historyFrequency));
evaluationCount = 0;
[~,lowResidual] = evaluate(low);
[~,highResidual] = evaluate(high);

expansionCount = 0;
while sameNonzeroSign(lowResidual,highResidual) && ...
        expansionCount < options.MaxBracketExpansions
    expansionCount = expansionCount+1;
    if lowResidual > 0
        high = low;
        highResidual = lowResidual;
        low = low/2;
        [~,lowResidual] = evaluate(low);
    else
        low = high;
        lowResidual = highResidual;
        high = high*2;
        [~,highResidual] = evaluate(high);
    end
end
if sameNonzeroSign(lowResidual,highResidual)
    error("radia:peec:ResonanceNotBracketed", ...
        ["Could not bracket the PEEC-capacitor reactance root after %d " ...
         "frequency expansions."],expansionCount);
end

converged = false;
rootFrequency = NaN;
rootEvaluation = struct();
rootResidual = NaN;
iterationCount = 0;
for iteration = 1:options.MaxIterations
    iterationCount = iteration;
    if highResidual == lowResidual
        candidate = 0.5*(low+high);
    else
        candidate = (low*highResidual-high*lowResidual)/ ...
            (highResidual-lowResidual);
    end
    margin = 0.1*(high-low);
    if ~isfinite(candidate) || candidate <= low+margin || ...
            candidate >= high-margin
        candidate = 0.5*(low+high);
    end
    [candidateEvaluation,candidateResidual] = evaluate(candidate);
    reactanceScale = max([abs(candidateEvaluation.reactance_Ohm), ...
        1/(2*pi*candidate*capacitance_F),eps]);
    residualRelative = abs(candidateResidual)/reactanceScale;
    bracketRelative = (high-low)/candidate;
    if residualRelative <= options.RelativeReactanceTolerance || ...
            bracketRelative <= options.RelativeFrequencyTolerance
        converged = true;
        rootFrequency = candidate;
        rootEvaluation = candidateEvaluation;
        rootResidual = candidateResidual;
        break
    end
    if sameNonzeroSign(lowResidual,candidateResidual)
        low = candidate;
        lowResidual = candidateResidual;
    else
        high = candidate;
        highResidual = candidateResidual;
    end
end
if ~converged
    error("radia:peec:ResonanceConvergence", ...
        "The self-consistent PEEC resonance did not converge in %d iterations.", ...
        options.MaxIterations);
end

capacitorReactance = -1/(2*pi*rootFrequency*capacitance_F);
reactanceScale = max([abs(rootEvaluation.reactance_Ohm), ...
    abs(capacitorReactance),eps]);
result = struct( ...
    "frequency_Hz",rootFrequency, ...
    "capacitance_F",capacitance_F, ...
    "impedance_Ohm",rootEvaluation.impedance_Ohm, ...
    "resistance_Ohm",rootEvaluation.resistance_Ohm, ...
    "inductance_H",rootEvaluation.inductance_H, ...
    "coil_reactance_Ohm",rootEvaluation.reactance_Ohm, ...
    "capacitor_reactance_Ohm",capacitorReactance, ...
    "net_reactance_Ohm",rootResidual, ...
    "relative_reactance_residual",abs(rootResidual)/reactanceScale, ...
    "dc_frequency_estimate_Hz",dcEstimate, ...
    "frequency_shift_from_dc_relative", ...
        (rootFrequency-dcEstimate)/dcEstimate, ...
    "iteration_count",iterationCount, ...
    "bracket_expansion_count",expansionCount, ...
    "evaluation_count",evaluationCount, ...
    "converged",converged, ...
    "frequency_dependent_properties",rootEvaluation, ...
    "history",struct( ...
        "frequency_Hz",historyFrequency(1:evaluationCount), ...
        "net_reactance_Ohm",historyResidual(1:evaluationCount), ...
        "resistance_Ohm",historyResistance(1:evaluationCount), ...
        "inductance_H",historyInductance(1:evaluationCount)));

    function [evaluation,residual] = evaluate(frequency)
        evaluation = radia.peec.seriesCoilImpedance( ...
            coil,frequency,DCProperties=dc);
        residual = evaluation.reactance_Ohm - ...
            1/(2*pi*frequency*capacitance_F);
        evaluationCount = evaluationCount+1;
        historyFrequency(evaluationCount) = frequency;
        historyResidual(evaluationCount) = residual;
        historyResistance(evaluationCount) = evaluation.resistance_Ohm;
        historyInductance(evaluationCount) = evaluation.inductance_H;
    end
end

function value = sameNonzeroSign(left,right)
value = left ~= 0 && right ~= 0 && sign(left) == sign(right);
end
