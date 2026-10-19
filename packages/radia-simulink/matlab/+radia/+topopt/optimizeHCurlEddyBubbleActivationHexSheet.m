function result=optimizeHCurlEddyBubbleActivationHexSheet(initialState,buildStepInputs,deformationFactory,evaluateObjective,rebuildHMatrix,cubitBackend,elementSizes,driverOptions)
%OPTIMIZEHCURLEDDYBUBBLEACTIVATIONHEXSHEET Activation-adjoint two-level driver.
arguments
    initialState (1,1) struct
    buildStepInputs (1,1) function_handle
    deformationFactory (1,1) function_handle
    evaluateObjective (1,1) function_handle
    rebuildHMatrix (1,1) function_handle
    cubitBackend
    elementSizes (:,1) double {mustBePositive,mustBeFinite}
    driverOptions (1,1) struct=struct
end
linearizeStep=@(state)localLinearize(state,buildStepInputs);
args=namedargs2cell(driverOptions);
result=radia.topopt.optimizeHexSheetTopology(initialState,linearizeStep, ...
    deformationFactory,evaluateObjective,rebuildHMatrix,cubitBackend, ...
    elementSizes,args{:});
end

function step=localLinearize(state,builder)
inputs=builder(state);
required=["operator","cellCurlGrams","frequenciesHz","rhs", ...
    "conductivity","area","lpOptions"];
if ~isstruct(inputs)||any(~isfield(inputs,required))
    error("radia:topopt:HCurlStepInputs", ...
        "buildStepInputs did not return the complete HCurl activation-LP contract.");
end
args=namedargs2cell(inputs.lpOptions);
step=radia.topopt.linearizeAndSolveHCurlActivationSheetJouleLP( ...
    state,inputs.operator,inputs.cellCurlGrams,inputs.frequenciesHz, ...
    inputs.rhs,inputs.conductivity,inputs.area,args{:});
end
