function result = solveAbeElementFillPlan( ...
        specificationResponse, requestedDifference, materialActive, ...
        elementVolumes, options)
%SOLVEABEELEMENTFILLPLAN Native bounded Abe/DUCAS material-fill solve.
%   The design coordinate is one signed fill fraction per element. Existing
%   iron has capacity [-1,0] and addable air has capacity [0,1]. The native
%   kernel reuses one ACA-QR-TSVD factor throughout the clipping iteration.
arguments
    specificationResponse (:,:) double
    requestedDifference (:,1) double
    materialActive (:,1) logical
    elementVolumes (:,1) double {mustBePositive}
    options.ElementIds (:,1) double = []
    options.CurrentFill (:,1) double = []
    options.FieldResponse (:,:) double = []
    options.ResidualPeakToPeak (1,1) double = -1
    options.ResidualRms (1,1) double = 1e-10
    options.MaxIterations (1,1) double {mustBeInteger,mustBePositive} = 64
    options.Relaxation (1,1) double {mustBePositive} = 1
    options.StagnationTolerance (1,1) double {mustBeNonnegative} = 1e-12
    options.RelativeSingularThreshold (1,1) double {mustBeNonnegative} = 1e-12
    options.Modes (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.MaxRank (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.AcaTolerance (1,1) double {mustBePositive} = 1e-8
    options.Method (1,1) string = "aca_qr_tsvd"
end

[rowCount, elementCount] = size(specificationResponse);
if rowCount == 0 || elementCount == 0
    error("radia:topopt:EmptyAbeResponse", ...
        "specificationResponse must be nonempty.");
end
if numel(requestedDifference) ~= rowCount || ...
        numel(materialActive) ~= elementCount || ...
        numel(elementVolumes) ~= elementCount
    error("radia:topopt:AbeDimensionMismatch", ...
        "Target, activity, and volume dimensions must match the response.");
end
if isempty(options.CurrentFill)
    currentFill = zeros(elementCount,1);
else
    currentFill = options.CurrentFill;
end
if numel(currentFill) ~= elementCount
    error("radia:topopt:AbeCurrentFillSize", ...
        "CurrentFill must have one value per response column.");
end
if isempty(options.FieldResponse)
    fieldResponse = zeros(0,elementCount);
else
    fieldResponse = options.FieldResponse;
end
if size(fieldResponse,2) ~= elementCount
    error("radia:topopt:AbeFieldResponseSize", ...
        "FieldResponse must have one column per response column.");
end
if options.ResidualPeakToPeak < 0 && options.ResidualRms < 0
    error("radia:topopt:AbeResidualTolerance", ...
        "At least one residual tolerance must be nonnegative.");
end
if options.Relaxation > 1
    error("radia:topopt:AbeRelaxation", ...
        "Relaxation must not exceed one.");
end
method = lower(options.Method);
if ~ismember(method,["aca_qr_tsvd","dense"])
    error("radia:topopt:AbeMethod", ...
        "Method must be 'aca_qr_tsvd' or 'dense'.");
end

result = radia.internal.callMex( ...
    'topopt.abe_element_fill_plan', double(specificationResponse), ...
    double(requestedDifference), logical(materialActive), ...
    double(elementVolumes), double(currentFill), double(fieldResponse), ...
    options.ResidualPeakToPeak, options.ResidualRms, ...
    options.MaxIterations, options.Relaxation, ...
    options.StagnationTolerance, options.RelativeSingularThreshold, ...
    options.Modes, options.MaxRank, options.AcaTolerance, char(method));
if isempty(options.ElementIds)
    result.element_ids = (1:elementCount).';
elseif numel(options.ElementIds) == elementCount
    result.element_ids = options.ElementIds;
else
    error("radia:topopt:AbeElementIdsSize", ...
        "ElementIds must have one value per response column.");
end
end
