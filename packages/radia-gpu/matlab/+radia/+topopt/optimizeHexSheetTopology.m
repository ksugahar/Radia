function result=optimizeHexSheetTopology(initialState,linearizeStep,deformationFactory,evaluateObjective,rebuildHMatrix,cubitBackend,elementSizes,options)
%OPTIMIZEHEXSHEETTOPOLOGY Two-level GetTrafo/Cubit topology optimization.
%   Five to twenty inner iterations keep activation continuous and update an
%   NGSolve VectorH1 deformation. Activation hysteresis suppresses topology
%   chatter, and sparse topology changes are accumulated until their fraction
%   or age reaches a Cubit batching threshold. H-matrix reconstruction occurs
%   exactly once after each successful Cubit rebuild and never in the inner loop.
arguments
    initialState (1,1) struct
    linearizeStep (1,1) function_handle
    deformationFactory (1,1) function_handle
    evaluateObjective (1,1) function_handle
    rebuildHMatrix (1,1) function_handle
    cubitBackend
    elementSizes (:,1) double {mustBePositive,mustBeFinite}
    options.InnerIterations (1,1) double {mustBeInteger,mustBePositive}=10
    options.MinimumInnerIterations (1,1) double {mustBeInteger,mustBePositive}=5
    options.MaxOuterIterations (1,1) double {mustBeInteger,mustBePositive}=10
    options.ObjectiveTolerance (1,1) double {mustBeNonnegative}=1e-4
    options.DesignTolerance (1,1) double {mustBeNonnegative}=1e-3
    options.ActivationThreshold (1,1) double {mustBeBetween(options.ActivationThreshold,0,1)}=0.5
    options.ActivationRemoveThreshold (1,1) double ...
        {mustBeBetween(options.ActivationRemoveThreshold,0,1)}=0.35
    options.ActivationRestoreThreshold (1,1) double ...
        {mustBeBetween(options.ActivationRestoreThreshold,0,1)}=0.65
    options.CubitBatchInterval (1,1) double {mustBeInteger,mustBePositive}=5
    options.CubitBatchFraction (1,1) double ...
        {mustBePositive,mustBeLessThanOrEqual(options.CubitBatchFraction,1)}=0.05
    options.MinimumScale (1,1) double {mustBePositive}=1/64
    options.Contraction (1,1) double {mustBePositive}=0.5
    options.MinimumJacobianRatio (1,1) double {mustBePositive}=0.2
    options.MaximumCondition (1,1) double {mustBePositive}=20
    options.RefineThreshold (1,1) double {mustBePositive}=0.25
    options.RebuildThreshold (1,1) double {mustBePositive}=0.5
    options.IntegrationOrder (1,1) double {mustBeInteger,mustBePositive}=2
    options.WorkDirectory (1,1) string="C:\temp\radia_hex_topopt"
    options.FinalizeTopology (1,1) logical=true
    options.ProgressFcn=[]
end
if options.InnerIterations<5||options.InnerIterations>20
    error("radia:topopt:InnerIterations", ...
        "InnerIterations must be between 5 and 20.");
end
if options.MinimumInnerIterations>options.InnerIterations
    error("radia:topopt:InnerIterations", ...
        "MinimumInnerIterations must not exceed InnerIterations.");
end
if options.ActivationRemoveThreshold>=options.ActivationRestoreThreshold
    error("radia:topopt:ActivationHysteresis", ...
        "Activation thresholds must satisfy remove < restore.");
end
required=["mesh","model","normal_displacement","thickness","activation","objective"];
missing=required(~isfield(initialState,required));
if ~isempty(missing)
    error("radia:topopt:State","Initial state is missing: %s",strjoin(missing,", "));
end
if ~isa(initialState.mesh,"radia.ngsolve.Mesh")
    error("radia:topopt:State","Initial state mesh must be radia.ngsolve.Mesh.");
end
state=normalizeState(initialState,elementSizes);
if ~isfolder(options.WorkDirectory), mkdir(options.WorkDirectory); end

history=table('Size',[0,14], ...
    'VariableTypes',{'double','double','double','double','double','double', ...
    'string','double','double','logical','double','logical','logical','string'}, ...
    'VariableNames',{'OuterIteration','InnerIteration','GlobalIteration', ...
    'ObjectiveBefore','ObjectiveAfter','AcceptedScale','Route', ...
    'MinimumJacobian','MaximumCondition','TopologyChanged', ...
    'PendingTopologyChanges','CubitRebuilt','HMatrixRebuilt','CubitReason'});
committedNormal=state.normal_displacement;
committedTopology=state.activation>=options.ActivationThreshold;
pendingTopology=false(size(committedTopology));
desiredTopology=committedTopology;
pendingCount=0; lastCubitIteration=0;
globalIteration=0; cubitCount=0; hmatrixCount=0; converged=false;

for outerIteration=1:options.MaxOuterIterations
    cubitRequired=false; cubitReason=""; batchConverged=false;
    for innerIteration=1:options.InnerIterations
        globalIteration=globalIteration+1;
        step=linearizeStep(state);
        if isfield(step,"update"), update=step.update; else, update=step; end
        [targetNormal,targetThickness,targetActivation]=readUpdate(update,elementSizes);
        oldNormal=state.normal_displacement;
        oldThickness=state.thickness;
        oldActivation=state.activation;
        relative=abs(targetNormal-oldNormal)./elementSizes;
        desiredTopology=committedTopology;
        desiredTopology(committedTopology & ...
            targetActivation<=options.ActivationRemoveThreshold)=false;
        desiredTopology(~committedTopology & ...
            targetActivation>=options.ActivationRestoreThreshold)=true;
        pendingTopology=desiredTopology~=committedTopology;
        pendingCount=nnz(pendingTopology);
        topologyCommitDue=pendingCount>0 && ( ...
            pendingCount/max(1,numel(pendingTopology))>=options.CubitBatchFraction || ...
            globalIteration-lastCubitIteration>=options.CubitBatchInterval);
        acceptance=radia.topopt.backtrackTrafoDeformation( ...
            state.mesh,deformationFactory,oldNormal,targetNormal,relative, ...
            MinimumScale=options.MinimumScale,Contraction=options.Contraction, ...
            MinimumJacobianRatio=options.MinimumJacobianRatio, ...
            MaximumCondition=options.MaximumCondition, ...
            RefineThreshold=options.RefineThreshold, ...
            RebuildThreshold=options.RebuildThreshold, ...
            IntegrationOrder=options.IntegrationOrder, ...
            TopologyChanged=topologyCommitDue);
        if acceptance.accepted
            scale=acceptance.scale;
        else
            % Cubit can safely realize the full target even when an in-place
            % VectorH1 deformation cannot pass the Jacobian gate.
            scale=1;
        end
        normal=oldNormal+scale*(targetNormal-oldNormal);
        thickness=oldThickness+scale*(targetThickness-oldThickness);
        activation=oldActivation+scale*(targetActivation-oldActivation);
        candidate=state;
        candidate.normal_displacement=normal;
        candidate.thickness=thickness;
        candidate.activation=activation;
        objectiveBefore=state.objective;
        candidate.objective=double(evaluateObjective(candidate));
        if ~isscalar(candidate.objective)||~isfinite(candidate.objective)
            error("radia:topopt:Objective","Topology objective must be finite and scalar.");
        end
        change=max([abs(normal-oldNormal);abs(thickness-oldThickness); ...
            abs(activation-oldActivation)]);
        relativeObjective=abs(candidate.objective-objectiveBefore)/ ...
            max(1,abs(objectiveBefore));
        route=string(acceptance.decision.route);
        explicitCubit=isfield(step,"requires_cubit")&&logical(step.requires_cubit);
        if ~acceptance.accepted
            cubitRequired=true; cubitReason="Trafo backtracking exhausted";
        elseif topologyCommitDue
            cubitRequired=true; cubitReason="batched activation hysteresis commit";
        elseif route~="ngsolve_deform"
            cubitRequired=true; cubitReason="Trafo quality requires mesh rebuild";
        elseif explicitCubit
            cubitRequired=true; cubitReason="linearization requested topology rebuild";
        end
        batchConverged=change<=options.DesignTolerance&& ...
            relativeObjective<=options.ObjectiveTolerance;
        state=candidate;
        history(end+1,:)={outerIteration,innerIteration,globalIteration, ...
            objectiveBefore,state.objective,scale,route, ...
            acceptance.decision.minimum_jacobian, ...
            acceptance.decision.maximum_condition,topologyCommitDue, ...
            pendingCount,false,false,""}; %#ok<AGROW>
        if ~isempty(options.ProgressFcn)
            options.ProgressFcn(struct("outer_iteration",outerIteration, ...
                "inner_iteration",innerIteration,"state",state, ...
                "acceptance",acceptance, ...
                "topology_changed",topologyCommitDue, ...
                "pending_topology_changes",pendingCount));
        end
        if cubitRequired
            break
        end
        if batchConverged&&innerIteration>=options.MinimumInnerIterations
            break
        end
    end

    if ~cubitRequired&&batchConverged&&options.FinalizeTopology&& ...
            any(abs(state.normal_displacement-committedNormal)>options.DesignTolerance)
        cubitRequired=true;
        cubitReason="converged deformation requires CAD commit";
    end

    if cubitRequired
        state.mesh.unsetDeformation();
        request=struct("schema","radia.topopt.cubit-hex-remesh/v1", ...
            "outer_iteration",outerIteration, ...
            "normal_displacement",state.normal_displacement, ...
            "thickness",state.thickness,"activation",state.activation, ...
            "desired_topology",desiredTopology, ...
            "pending_topology_changes",pendingCount, ...
            "journal_path",fullfile(options.WorkDirectory, ...
                sprintf("hex_topopt_%04d.jou",outerIteration-1)), ...
            "mesh_path",fullfile(options.WorkDirectory, ...
                sprintf("hex_topopt_%04d.vol",outerIteration-1)), ...
            "reason",cubitReason);
        state.mesh=invokeCubit(cubitBackend,request);
        state.model=rebuildHMatrix(state.mesh,state.normal_displacement, ...
            state.thickness,state.activation,"cubit_rebuild");
        state.objective=double(evaluateObjective(state));
        if ~isscalar(state.objective)||~isfinite(state.objective)
            error("radia:topopt:Objective", ...
                "Post-Cubit topology objective must be finite and scalar.");
        end
        cubitCount=cubitCount+1; hmatrixCount=hmatrixCount+1;
        committedNormal=state.normal_displacement;
        committedTopology=desiredTopology;
        pendingTopology(:)=false;
        lastCubitIteration=globalIteration;
        history.CubitRebuilt(end)=true;
        history.HMatrixRebuilt(end)=true;
        history.CubitReason(end)=cubitReason;
        history.Route(end)="cubit_rebuild";
        history.ObjectiveAfter(end)=state.objective;
    end
    if batchConverged
        converged=true;
        break
    end
end

result=struct("schema","radia.topopt.hex-sheet-two-level/v2", ...
    "state",state,"history",history,"converged",converged, ...
    "inner_iteration_count",globalIteration, ...
    "outer_iteration_count",max([0;history.OuterIteration]), ...
    "cubit_rebuild_count",cubitCount, ...
    "hmatrix_rebuild_count",hmatrixCount, ...
    "inner_iterations_per_batch",options.InnerIterations, ...
    "committed_topology",committedTopology, ...
    "pending_topology_changes",nnz(pendingTopology));
end

function state=normalizeState(state,elementSizes)
state.normal_displacement=state.normal_displacement(:);
state.thickness=state.thickness(:);
state.activation=state.activation(:);
if ~isequal(size(state.normal_displacement),size(elementSizes))|| ...
        ~isequal(size(state.thickness),size(elementSizes))|| ...
        ~isequal(size(state.activation),size(elementSizes))
    error("radia:topopt:Shape", ...
        "HEX sheet design arrays must match ElementSizes.");
end
if any(~isfinite(state.normal_displacement))||any(~isfinite(state.thickness))|| ...
        any(~isfinite(state.activation))||any(state.thickness<=0)|| ...
        any(state.activation<0|state.activation>1)|| ...
        ~isscalar(state.objective)||~isfinite(state.objective)
    error("radia:topopt:State","Initial topology state is not physically valid.");
end
state.objective=double(state.objective);
end

function [normal,thickness,activation]=readUpdate(update,elementSizes)
required=["normal_displacement","thickness","activation"];
if ~isstruct(update)||any(~isfield(update,required))
    error("radia:topopt:Update", ...
        "LinearizeStep must return an update with displacement, thickness, and activation.");
end
normal=update.normal_displacement(:);
thickness=update.thickness(:);
activation=update.activation(:);
if ~isequal(size(normal),size(elementSizes))|| ...
        ~isequal(size(thickness),size(elementSizes))|| ...
        ~isequal(size(activation),size(elementSizes))
    error("radia:topopt:Shape","Topology update arrays must match ElementSizes.");
end
if any(~isfinite(normal))||any(~isfinite(thickness))|| ...
        any(~isfinite(activation))||any(thickness<=0)|| ...
        any(activation<0|activation>1)
    error("radia:topopt:Update","Topology update is not physically valid.");
end
end

function mesh=invokeCubit(backend,request)
if isa(backend,"function_handle")
    mesh=backend(request);
elseif isa(backend,"radia.topopt.CubitHexRemeshBackend")
    mesh=backend.rebuild(request);
else
    error("radia:topopt:CubitBackend", ...
        "Cubit backend must be a function handle or CubitHexRemeshBackend.");
end
if ~isa(mesh,"radia.ngsolve.Mesh")
    error("radia:topopt:CubitMesh", ...
        "Cubit backend must return a radia.ngsolve.Mesh.");
end
end
