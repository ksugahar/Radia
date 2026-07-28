function model = makeCircuitFieldStateSpace(K, sourceMatrix, resistance_ohm, options)
%MAKECIRCUITFIELDSTATESPACE Reduce a linear 2-D field/circuit model for MEX.
%   K is the constrained planar or axisymmetric magnetic field matrix and
%   sourceMatrix contains one signed coil source vector per branch. The
%   reduced inductance is L = sourceMatrix' * (K \ sourceMatrix), giving
%       L * di/dt + R * i = inputMap * voltage(t).
%   The exact-ZOH model runs through radia_state_space_mex_sfunction, so
%   Simulink performs no Python or finite-element factorization per step.

arguments
    K double
    sourceMatrix double
    resistance_ohm double
    options.SampleTime_s (1,1) double {mustBeFinite, mustBePositive} = 1.0e-5
    options.VoltageInputMode (1,1) string {mustBeMember(options.VoltageInputMode, ...
        ["common", "per_branch"])} = "common"
    options.InitialCurrent_A double = []
end

K = double(K);
sourceMatrix = double(sourceMatrix);
if ndims(K) ~= 2 || size(K, 1) ~= size(K, 2) || isempty(K) || any(~isfinite(K), "all")
    error("radia:simulink:CircuitFieldMatrix", ...
        "K must be a non-empty finite square matrix.");
end
if ndims(sourceMatrix) ~= 2 || size(sourceMatrix, 1) ~= size(K, 1) || ...
        isempty(sourceMatrix) || any(~isfinite(sourceMatrix), "all")
    error("radia:simulink:CircuitSourceMatrix", ...
        "sourceMatrix must have one row per field degree of freedom.");
end

branchCount = size(sourceMatrix, 2);
resistance = double(resistance_ohm(:));
if numel(resistance) ~= branchCount || any(~isfinite(resistance)) || any(resistance <= 0)
    error("radia:simulink:CircuitResistance", ...
        "resistance_ohm must contain one positive finite value per branch.");
end

K = 0.5 * (K + K.');
if min(eig(K)) <= 0 || rcond(K) <= eps
    error("radia:simulink:CircuitFieldMatrix", ...
        "The constrained field matrix must be positive definite.");
end
inductance = sourceMatrix.' * (K \ sourceMatrix);
inductance = 0.5 * (inductance + inductance.');
if min(eig(inductance)) <= 0 || rcond(inductance) <= eps
    error("radia:simulink:CircuitInductance", ...
        "The reduced branch inductance matrix must be positive definite.");
end

if options.VoltageInputMode == "common"
    inputMap = ones(branchCount, 1);
else
    inputMap = eye(branchCount);
end
A = -(inductance \ diag(resistance));
B = inductance \ inputMap;
C = [eye(branchCount); inductance];
D = zeros(2 * branchCount, size(inputMap, 2));

if isempty(options.InitialCurrent_A)
    x0 = zeros(branchCount, 1);
else
    x0 = double(options.InitialCurrent_A(:));
    if numel(x0) ~= branchCount || any(~isfinite(x0))
        error("radia:simulink:CircuitInitialCurrent", ...
            "InitialCurrent_A must contain one finite value per branch.");
    end
end

inputCount = size(B, 2);
augmented = [A, B; zeros(inputCount, branchCount + inputCount)];
discrete = expm(augmented * options.SampleTime_s);
Ad = discrete(1:branchCount, 1:branchCount);
Bd = discrete(1:branchCount, branchCount + (1:inputCount));

model = struct( ...
    "schema", "radia.circuit-field.state-space.v1", ...
    "backend", "native-mex-sfunction", ...
    "field_matrix", K, "source_matrix", sourceMatrix, ...
    "resistance_ohm", resistance, "inductance_H", inductance, ...
    "A", A, "B", B, "C", C, "D", D, ...
    "Ad", Ad, "Bd", Bd, "Cd", C, "Dd", D, ...
    "x0", x0, "sample_time_s", options.SampleTime_s, ...
    "voltage_input_mode", options.VoltageInputMode, ...
    "state_order", branchCount, "input_count", inputCount, ...
    "output_count", 2 * branchCount, ...
    "output_convention", "branch currents followed by branch flux linkages", ...
    "mex_s_function", "radia_state_space_mex_sfunction", ...
    "python_per_step", false, "field_factorization_per_step", false);
end
