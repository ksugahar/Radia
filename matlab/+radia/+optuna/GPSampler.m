classdef GPSampler < radia.optuna.BaseSampler
    %GPSAMPLER Matérn-5/2 ARD Bayesian optimization for CAE objectives.
    %   Supports mixed stable search spaces, single-objective expected
    %   improvement, multi-objective expected hypervolume improvement,
    %   soft constraints c<=0, and pending-trial repulsion. The hot loop is
    %   MATLAB-native and does not require Statistics or Optimization Toolbox.

    properties (SetAccess=private)
        Seed (1,1) double = 0
        NStartupTrials (1,1) double = 10
        DeterministicObjective (1,1) logical = false
        Backend (1,1) string = "upstream-python"
        CandidateCount (1,1) double = 2048
        LocalSearchCount (1,1) double = 10
        MonteCarloSamples (1,1) double = 128
        ConstraintsFcn = []
        Stream
    end

    properties (Access=private)
        IndependentSampler
        AttachedStudy = []
        Restored (1,1) logical = false
        ObjectiveTheta cell = cell(0,1)
        ConstraintTheta cell = cell(0,1)
    end

    properties (Transient, Access=private)
        PythonOptuna = []
        PythonStudy = []
        PythonTrial = []
        PythonTrialNumber (1,1) double = NaN
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.gp-sampler-state.v1"
        SamplerName = "gp"
    end

    methods
        function obj=GPSampler(options)
            arguments
                options.Seed double = double.empty(1,0)
                options.NStartupTrials (1,1) double ...
                    {mustBeInteger,mustBeNonnegative} = 10
                options.DeterministicObjective (1,1) logical = false
                options.Backend (1,1) string = "upstream-python"
                options.CandidateCount (1,1) double ...
                    {mustBeInteger,mustBePositive} = 2048
                options.LocalSearchCount (1,1) double ...
                    {mustBeInteger,mustBeNonnegative} = 10
                options.MonteCarloSamples (1,1) double ...
                    {mustBeInteger,mustBePositive} = 128
                options.ConstraintsFcn = []
            end
            if options.CandidateCount<16
                error("radia:optuna:GPCandidates", ...
                    "CandidateCount must be at least 16.");
            end
            if ~isempty(options.ConstraintsFcn) && ...
                    ~isa(options.ConstraintsFcn,"function_handle")
                error("radia:optuna:GPConstraints", ...
                    "ConstraintsFcn must be a function handle.");
            end
            if ~ismember(options.Backend,["upstream-python","matlab-native"])
                error("radia:optuna:GPBackend", ...
                    "Backend must be 'upstream-python' or 'matlab-native'.");
            end
            obj.Seed=radia.optuna.internal.resolveSeed(options.Seed);
            obj.NStartupTrials=double(options.NStartupTrials);
            obj.DeterministicObjective=options.DeterministicObjective;
            obj.Backend=options.Backend;
            obj.CandidateCount=double(options.CandidateCount);
            obj.LocalSearchCount=double(options.LocalSearchCount);
            obj.MonteCarloSamples=double(options.MonteCarloSamples);
            obj.ConstraintsFcn=options.ConstraintsFcn;
            obj.Stream=radia.optuna.internal.NumpyRandomState(obj.Seed);
            obj.IndependentSampler=radia.optuna.RandomSampler(options.Seed);
        end

        function searchSpace=inferRelativeSearchSpace(~,study,trial) %#ok<INUSD>
            searchSpace=radia.optuna.internal.IntersectionSearchSpace. ...
                calculate(study,IncludePruned=false);
        end

        function searchSpace=infer_relative_search_space(obj,study,trial)
            if nargin<3, trial=[]; end
            searchSpace=obj.inferRelativeSearchSpace(study,trial);
        end

        function before_trial(obj,study,trial)
            % Upstream's public hook delegates only to its independent
            % sampler. Python-study mirroring belongs to MATLAB's internal
            % ask lifecycle in beforeTrial, not to this public hook.
            obj.IndependentSampler.before_trial(study,trial);
        end

        function after_trial(obj,study,trial,state,values)
            if ~isempty(obj.ConstraintsFcn) && ...
                    state==radia.optuna.TrialState.COMPLETE
                study.recordConstraints(trial,obj.ConstraintsFcn(trial));
            end
            obj.IndependentSampler.after_trial( ...
                study,trial,state,values);
        end

        function beforeTrial(obj,study,trial)
            if obj.Backend=="upstream-python" && ...
                    ~isempty(obj.PythonTrial) && ...
                    obj.PythonTrialNumber==trial.Number
                obj.attach(study);
                obj.preparePythonTrial(study,trial);
                trial.setSystemAttr("gp_sampling_mode","upstream_optuna_4_9_0");
                trial.setSystemAttr("gp_backend","upstream-python");
                return
            end
            obj.IndependentSampler.beforeTrial(study,trial);
            obj.attach(study);
            completed=sum(study.TrialTable.State=="COMPLETE");
            if completed<obj.NStartupTrials
                trial.setSystemAttr("gp_sampling_mode","startup_random");
                obj.recordState(study,trial.Number);
                return
            end
            searchSpace=obj.inferRelativeSearchSpace(study,trial);
            if isempty(searchSpace)
                trial.setSystemAttr("gp_sampling_mode", ...
                    "independent_dynamic_space");
                obj.recordState(study,trial.Number);
                return
            end
            [observations,objectives,numbers,pending,constraints, ...
                constraintPresent]=obj.observations(study,searchSpace,trial.Number);
            finished=~pending;
            if sum(finished)<obj.NStartupTrials || ...
                    size(objectives,1)<obj.NStartupTrials
                trial.setSystemAttr("gp_sampling_mode","startup_random");
                obj.recordState(study,trial.Number);
                return
            end
            values=obj.sampleRelative(study,searchSpace,observations, ...
                objectives,numbers,pending,constraints,constraintPresent);
            trial.setRelativeParameters(searchSpace,values,"gp");
            if isscalar(study.Directions)
                acquisition="expected_improvement";
            else
                acquisition="expected_hypervolume_improvement";
            end
            trial.setSystemAttr("gp_sampling_mode","matern52_ard");
            trial.setSystemAttr("gp_acquisition",acquisition);
            trial.setSystemAttr("gp_pending_count",sum(pending));
            obj.recordState(study,trial.Number);
        end

        function value=sampleFloat(obj,study,trial,name,low,high,options)
            if obj.Backend=="upstream-python" && ...
                    ~isempty(obj.PythonTrial) && ...
                    obj.PythonTrialNumber==trial.Number
                obj.ensurePythonTrial(trial);
                arguments={"log",logical(options.Log)};
                if isfinite(options.Step)
                    arguments=[arguments,{"step",double(options.Step)}];
                end
                value=double(obj.PythonTrial.suggest_float( ...
                    char(name),double(low),double(high),pyargs(arguments{:})));
                return
            end
            value=obj.IndependentSampler.sampleFloat( ...
                study,trial,name,low,high,options);
        end

        function value=sampleInteger(obj,study,trial,name,low,high)
            if obj.Backend=="upstream-python"
                value=obj.sampleIntegerDetailed( ...
                    study,trial,name,low,high,1,false);
                return
            end
            value=obj.IndependentSampler.sampleInteger( ...
                study,trial,name,low,high);
        end

        function value=sampleIntegerDetailed(obj,study,trial,name,low,high,step,logScale)
            if obj.Backend~="upstream-python" || isempty(obj.PythonTrial) || ...
                    obj.PythonTrialNumber~=trial.Number
                if ~logScale && step==1
                    value=obj.IndependentSampler.sampleInteger( ...
                        study,trial,name,low,high);
                else
                    value=obj.IndependentSampler.sampleFloat( ...
                        study,trial,name,low,high, ...
                        struct("Log",logical(logScale),"Step",double(step)));
                    value=low+round((value-low)/step)*step;
                end
                return
            end
            obj.ensurePythonTrial(trial);
            value=double(obj.PythonTrial.suggest_int( ...
                char(name),int64(low),int64(high), ...
                pyargs("step",int64(step),"log",logical(logScale))));
        end

        function value=sampleCategorical(obj,study,trial,name,choices)
            if obj.Backend=="upstream-python" && ...
                    ~isempty(obj.PythonTrial) && ...
                    obj.PythonTrialNumber==trial.Number
                obj.ensurePythonTrial(trial);
                raw=obj.PythonTrial.suggest_categorical( ...
                    char(name),obj.pythonChoices(choices));
                value=obj.matlabChoice(raw);
                return
            end
            value=obj.IndependentSampler.sampleCategorical( ...
                study,trial,name,choices);
        end

        function values=sampleJoint(obj,study,trial,names,lows,highs,options)
            values=zeros(1,numel(names));
            for index=1:numel(names)
                if obj.Backend=="upstream-python"
                    values(index)=obj.sampleFloat( ...
                        study,trial,names(index),lows(index),highs(index), ...
                        struct("Log",options.Log(index),"Step",NaN));
                else
                    values(index)=obj.IndependentSampler.sampleFloat( ...
                        study,trial,names(index),lows(index),highs(index), ...
                        struct("Log",options.Log(index),"Step",NaN));
                end
            end
        end

        function afterTrial(obj,study,trial)
            if obj.Backend=="upstream-python"
                constraintPresent=false;
                constraints=zeros(1,0);
                if trial.State=="COMPLETE" && ~isempty(obj.ConstraintsFcn)
                    [constraintPresent,constraints]= ...
                        study.constraintRecord(trial.Number);
                    if ~constraintPresent
                        constraints=reshape(double(obj.ConstraintsFcn(trial)),1,[]);
                        study.recordConstraints(trial,constraints);
                        constraintPresent=true;
                    end
                end
                obj.finishPythonTrial(trial,constraintPresent,constraints);
                return
            end
            obj.IndependentSampler.afterTrial(study,trial);
            if trial.State=="COMPLETE" && ~isempty(obj.ConstraintsFcn)
                study.recordConstraints(trial,obj.ConstraintsFcn(trial));
            end
            obj.recordState(study,trial.Number);
        end
    end

    methods (Access=private)
        function preparePythonTrial(obj,study,trial)
            if isempty(obj.PythonStudy)
                obj.initializePythonStudy(study);
                obj.replayPythonHistory(study,trial.Number);
            end
            if ~isempty(obj.PythonTrial)
                error("radia:optuna:GPPythonTrial", ...
                    "The previous upstream Python GP trial is still running.");
            end
            obj.PythonTrial=obj.PythonStudy.ask();
            obj.PythonTrialNumber=trial.Number;
        end

        function initializePythonStudy(obj,study)
            environment=pyenv;
            if environment.Status=="NotLoaded"
                try
                    environment=pyenv(ExecutionMode="InProcess");
                catch exception
                    cause=MException("radia:optuna:GPPython", ...
                        "GPSampler Backend='upstream-python' requires the " + ...
                         "configured in-process Python 3.12 environment.");
                    throw(addCause(cause,exception));
                end
            end
            if environment.ExecutionMode~="InProcess" || ...
                    ~startsWith(string(environment.Version),"3.12")
                error("radia:optuna:GPPython", ...
                    "GPSampler Backend='upstream-python' requires " + ...
                     "in-process Python 3.12; pyenv reports %s (%s).", ...
                    string(environment.Version),string(environment.ExecutionMode));
            end
            try
                obj.PythonOptuna=py.importlib.import_module("optuna");
                version=string(py.builtins.getattr( ...
                    obj.PythonOptuna,"__version__"));
                if version~="4.9.0"
                    error("radia:optuna:GPPythonVersion", ...
                        "Expected optuna==4.9.0, found %s.",version);
                end
                samplerArguments=obj.pythonSamplerArguments();
                sampler=obj.PythonOptuna.samplers.GPSampler( ...
                    pyargs(samplerArguments{:}));
                if isscalar(study.Directions)
                    obj.PythonStudy=obj.PythonOptuna.create_study(pyargs( ...
                        "direction",char(study.Directions),"sampler",sampler));
                else
                    directions=py.list(cellstr(study.Directions));
                    obj.PythonStudy=obj.PythonOptuna.create_study(pyargs( ...
                        "directions",directions,"sampler",sampler));
                end
            catch exception
                if exception.identifier=="radia:optuna:GPPythonVersion"
                    rethrow(exception)
                end
                cause=MException("radia:optuna:GPPython", ...
                    "Could not initialize the pinned Optuna 4.9.0 GP " + ...
                     "backend with its NumPy, SciPy, and PyTorch runtime.");
                throw(addCause(cause,exception));
            end
        end

        function ensurePythonTrial(obj,trial)
            if isempty(obj.PythonTrial) || obj.PythonTrialNumber~=trial.Number
                error("radia:optuna:GPPythonTrial", ...
                    "No upstream Python GP trial is active for trial %d.", ...
                    trial.Number);
            end
        end

        function finishPythonTrial(obj,trial,constraintPresent,constraints)
            obj.ensurePythonTrial(trial);
            if nargin<3, constraintPresent=false; end
            if nargin<4, constraints=zeros(1,0); end
            if trial.State=="COMPLETE" && constraintPresent
                obj.setPythonConstraints(obj.PythonTrial,constraints);
            end
            switch trial.State
                case "COMPLETE"
                    values=reshape(double(trial.Values),1,[]);
                    if isscalar(values)
                        obj.PythonStudy.tell(obj.PythonTrial,values);
                    else
                        obj.PythonStudy.tell(obj.PythonTrial, ...
                            py.list(num2cell(values)));
                    end
                case "PRUNED"
                    obj.PythonStudy.tell(obj.PythonTrial,pyargs( ...
                        "state",obj.PythonOptuna.trial.TrialState.PRUNED));
                case "FAIL"
                    obj.PythonStudy.tell(obj.PythonTrial,pyargs( ...
                        "state",obj.PythonOptuna.trial.TrialState.FAIL));
                otherwise
                    error("radia:optuna:GPPythonTrial", ...
                        "Unsupported upstream GP trial state '%s'.",trial.State);
            end
            obj.PythonTrial=[];
            obj.PythonTrialNumber=NaN;
        end

        function pairs=pythonSamplerArguments(obj)
            pairs={ ...
                "seed",int64(obj.Seed), ...
                "n_startup_trials",int64(obj.NStartupTrials), ...
                "deterministic_objective",obj.DeterministicObjective};
            if ~isempty(obj.ConstraintsFcn)
                code=join([ ...
                    "def radia_matlab_constraints(trial):"; ...
                    "    return trial.user_attrs['__radia_matlab_constraints']"; ...
                    "result = radia_matlab_constraints"],newline);
                namespace=py.dict;
                py.builtins.exec(char(code),namespace);
                callback=namespace{char("result")};
                pairs=[pairs,{"constraints_func",callback}];
            end
        end

        function replayPythonHistory(obj,study,beforeTrialNumber)
            try
                obj.replayPythonHistoryExactly(study,beforeTrialNumber);
            catch exception
                if exception.identifier~="radia:optuna:GPPythonReplayMismatch"
                    rethrow(exception)
                end
                % A study created by another sampler has valid observations
                % but cannot reproduce this GP sampler's random draws. Start
                % from a fresh upstream sampler and import those exact frozen
                % trials instead of rejecting a public Optuna study history.
                obj.initializePythonStudy(study);
                obj.importPythonHistory(study,beforeTrialNumber);
            end
        end

        function replayPythonHistoryExactly(obj,study,beforeTrialNumber)
            rows=study.TrialTable.TrialNumber<beforeTrialNumber;
            prior=sortrows(study.TrialTable(rows,:),"TrialNumber");
            for index=1:height(prior)
                number=prior.TrialNumber(index);
                state=prior.State(index);
                if ~ismember(state,["COMPLETE","PRUNED","FAIL"])
                    error("radia:optuna:GPPythonReplay", ...
                        "Cannot replay prior trial %d in state %s.",number,state);
                end
                pythonTrial=obj.PythonStudy.ask();
                if double(pythonTrial.number)~=number
                    error("radia:optuna:GPPythonReplay", ...
                        "Python replay trial number %d does not match MATLAB trial %d.", ...
                        double(pythonTrial.number),number);
                end
                parameterRows=find(study.ParamTable.TrialNumber==number)';
                for row=parameterRows
                    distribution= ...
                        radia.optuna.internal.DistributionCodec.decode( ...
                        study.ParamTable.Kind(row), ...
                        study.ParamTable.Distribution(row));
                    proposed=obj.pythonSuggest( ...
                        pythonTrial,study.ParamTable.Name(row),distribution);
                    recorded=obj.parameterValue(study.ParamTable(row,:));
                    if ~obj.sameParameterValue(proposed,recorded)
                        error("radia:optuna:GPPythonReplayMismatch", ...
                            "Replayed parameter '%s' in trial %d differs from the stored value.", ...
                            study.ParamTable.Name(row),number);
                    end
                end
                intermediate=sortrows(study.IntermediateTable( ...
                    study.IntermediateTable.TrialNumber==number,:),"Step");
                for reportIndex=1:height(intermediate)
                    pythonTrial.report(intermediate.Value(reportIndex), ...
                        int64(intermediate.Step(reportIndex)));
                end
                [constraintPresent,constraints]=study.constraintRecord(number);
                if state=="COMPLETE" && constraintPresent && ...
                        ~isempty(obj.ConstraintsFcn)
                    obj.setPythonConstraints(pythonTrial,constraints);
                end
                obj.tellPythonReplay(study,pythonTrial,number,state);
            end
        end

        function importPythonHistory(obj,study,beforeTrialNumber)
            rows=study.TrialTable.TrialNumber<beforeTrialNumber;
            prior=sortrows(study.TrialTable(rows,:),"TrialNumber");
            for index=1:height(prior)
                number=prior.TrialNumber(index);
                state=prior.State(index);
                if ~ismember(state,["COMPLETE","PRUNED","FAIL"])
                    error("radia:optuna:GPPythonReplay", ...
                        "Cannot replay prior trial %d in state %s.",number,state);
                end
                params=py.dict;
                distributions=py.dict;
                parameterRows=find(study.ParamTable.TrialNumber==number)';
                for row=parameterRows
                    name=char(study.ParamTable.Name(row));
                    distribution= ...
                        radia.optuna.internal.DistributionCodec.decode( ...
                        study.ParamTable.Kind(row), ...
                        study.ParamTable.Distribution(row));
                    recorded=obj.parameterValue(study.ParamTable(row,:));
                    params{name}=obj.pythonValue(recorded);
                    distributions{name}=obj.pythonDistribution(distribution);
                end
                intermediate=sortrows(study.IntermediateTable( ...
                    study.IntermediateTable.TrialNumber==number,:),"Step");
                intermediateValues=py.dict;
                for reportIndex=1:height(intermediate)
                    intermediateValues{int64(intermediate.Step(reportIndex))}= ...
                        intermediate.Value(reportIndex);
                end
                userAttrs=obj.pythonAttributes( ...
                    study.UserAttrTable,number);
                systemAttrs=obj.pythonAttributes( ...
                    study.SystemAttrTable,number);
                [constraintPresent,constraints]=study.constraintRecord(number);
                if state=="COMPLETE" && constraintPresent && ...
                        ~isempty(obj.ConstraintsFcn)
                    userAttrs{char("__radia_matlab_constraints")}= ...
                        py.list(num2cell(reshape(double(constraints),1,[])));
                end
                pairs={"state",obj.pythonTrialState(state), ...
                    "params",params,"distributions",distributions, ...
                    "user_attrs",userAttrs,"system_attrs",systemAttrs, ...
                    "intermediate_values",intermediateValues};
                objectiveRows=sortrows(study.ObjectiveTable( ...
                    study.ObjectiveTable.TrialNumber==number,:), ...
                    "ObjectiveIndex");
                if ~isempty(objectiveRows)
                    values=reshape(double(objectiveRows.Value),1,[]);
                    if isscalar(values)
                        pairs=[pairs,{"value",values}];
                    else
                        pairs=[pairs,{"values",py.list(num2cell(values))}];
                    end
                end
                frozen=obj.PythonOptuna.trial.create_trial( ...
                    pyargs(pairs{:}));
                obj.PythonStudy.add_trial(frozen);
                imported=obj.PythonStudy.trials{int64(number+1)};
                if double(imported.number)~=number
                    error("radia:optuna:GPPythonReplay", ...
                        "Imported Python trial number %d does not match MATLAB trial %d.", ...
                        double(imported.number),number);
                end
            end
        end

        function result=pythonDistribution(obj,distribution)
            switch distribution.kind
                case "float"
                    pairs={"log",logical(distribution.log)};
                    if isfinite(distribution.step)
                        pairs=[pairs,{"step",distribution.step}];
                    end
                    result=obj.PythonOptuna.distributions.FloatDistribution( ...
                        distribution.low,distribution.high,pyargs(pairs{:}));
                case "integer"
                    result=obj.PythonOptuna.distributions.IntDistribution( ...
                        int64(distribution.low),int64(distribution.high), ...
                        pyargs("log",logical(distribution.log), ...
                        "step",int64(distribution.step)));
                case "categorical"
                    result=obj.PythonOptuna.distributions. ...
                        CategoricalDistribution( ...
                        obj.pythonChoices(distribution.choices));
                otherwise
                    error("radia:optuna:GPPythonReplay", ...
                        "Unsupported stored distribution kind '%s'.", ...
                        distribution.kind);
            end
        end

        function result=pythonAttributes(~,source,number)
            result=py.dict;
            json=py.importlib.import_module("json");
            rows=find(source.TrialNumber==number)';
            for row=rows
                result{char(source.Name(row))}= ...
                    json.loads(char(source.ValueJSON(row)));
            end
        end

        function result=pythonTrialState(obj,state)
            switch state
                case "COMPLETE"
                    result=obj.PythonOptuna.trial.TrialState.COMPLETE;
                case "PRUNED"
                    result=obj.PythonOptuna.trial.TrialState.PRUNED;
                case "FAIL"
                    result=obj.PythonOptuna.trial.TrialState.FAIL;
                otherwise
                    error("radia:optuna:GPPythonReplay", ...
                        "Cannot import trial state '%s'.",state);
            end
        end

        function result=pythonValue(~,value)
            if isstring(value) || ischar(value)
                result=py.builtins.str(char(string(value)));
            elseif islogical(value) && isscalar(value)
                result=py.builtins.bool(value);
            elseif isnumeric(value) && isscalar(value)
                result=py.builtins.float(double(value));
            else
                error("radia:optuna:GPPythonReplay", ...
                    "Only scalar string, logical, and numeric parameters are supported.");
            end
        end

        function value=pythonSuggest(obj,pythonTrial,name,distribution)
            switch distribution.kind
                case "float"
                    args={"log",logical(distribution.log)};
                    if isfinite(distribution.step)
                        args=[args,{"step",double(distribution.step)}];
                    end
                    value=double(pythonTrial.suggest_float( ...
                        char(name),distribution.low,distribution.high, ...
                        pyargs(args{:})));
                case "integer"
                    value=double(pythonTrial.suggest_int( ...
                        char(name),int64(distribution.low), ...
                        int64(distribution.high),pyargs( ...
                        "step",int64(distribution.step), ...
                        "log",logical(distribution.log))));
                case "categorical"
                    value=obj.matlabChoice(pythonTrial.suggest_categorical( ...
                        char(name),obj.pythonChoices(distribution.choices)));
                otherwise
                    error("radia:optuna:GPPythonReplay", ...
                        "Unsupported stored distribution kind '%s'.", ...
                        distribution.kind);
            end
        end

        function value=parameterValue(~,row)
            if ~isnan(row.ValueNumeric)
                value=row.ValueNumeric;
            else
                value=jsondecode(row.ValueText);
            end
        end

        function result=sameParameterValue(~,left,right)
            if isnumeric(left) && isnumeric(right) && ...
                    isscalar(left) && isscalar(right)
                result=isequaln(double(left),double(right));
            else
                result=radia.optuna.internal.DistributionCodec. ...
                    choiceToken(left)== ...
                    radia.optuna.internal.DistributionCodec.choiceToken(right);
            end
        end

        function setPythonConstraints(~,pythonTrial,constraints)
            pythonTrial.set_user_attr(char("__radia_matlab_constraints"), ...
                py.list(num2cell(reshape(double(constraints),1,[]))));
        end

        function tellPythonReplay(obj,study,pythonTrial,number,state)
            switch state
                case "COMPLETE"
                    rows=study.ObjectiveTable.TrialNumber==number;
                    values=sortrows(study.ObjectiveTable(rows,:), ...
                        "ObjectiveIndex").Value;
                    values=reshape(double(values),1,[]);
                    if isscalar(values)
                        obj.PythonStudy.tell(pythonTrial,values);
                    else
                        obj.PythonStudy.tell(pythonTrial,py.list(num2cell(values)));
                    end
                case "PRUNED"
                    obj.PythonStudy.tell(pythonTrial,pyargs( ...
                        "state",obj.PythonOptuna.trial.TrialState.PRUNED));
                case "FAIL"
                    obj.PythonStudy.tell(pythonTrial,pyargs( ...
                        "state",obj.PythonOptuna.trial.TrialState.FAIL));
            end
        end

        function result=pythonChoices(~,choices)
            count=numel(choices);
            items=cell(1,count);
            for index=1:count
                if iscell(choices)
                    item=choices{index};
                else
                    item=choices(index);
                end
                if isstring(item) || ischar(item)
                    items{index}=py.str(char(string(item)));
                elseif islogical(item)
                    items{index}=py.bool(logical(item));
                elseif isnumeric(item) && isscalar(item)
                    items{index}=py.float(double(item));
                else
                    error("radia:optuna:GPChoices", ...
                        "The upstream GP backend supports scalar string, logical, and numeric choices.");
                end
            end
            result=py.list(items);
        end

        function value=matlabChoice(~,raw)
            if isa(raw,"py.str")
                value=string(raw);
            elseif isa(raw,"py.bool")
                value=logical(raw);
            elseif isa(raw,"py.int") || isa(raw,"py.float")
                value=double(raw);
            else
                error("radia:optuna:GPChoices", ...
                    "Optuna returned an unsupported categorical choice type '%s'.", ...
                    string(class(raw)));
            end
        end

        function values=sampleRelative(obj,study,searchSpace,observations, ...
                objectives,trialNumbers,pending,constraints,constraintPresent)
            finished=~pending;
            x=observations(finished,:);
            y=objectives(finished,:);
            finishedNumbers=trialNumbers(finished);
            categorical=obj.categoricalMask(searchSpace);
            candidates=obj.candidatePool(x,categorical,searchSpace);
            objectiveModels=cell(1,size(y,2));
            if numel(obj.ObjectiveTheta)~=size(y,2)
                obj.ObjectiveTheta=cell(size(y,2),1);
            end
            for objective=1:size(y,2)
                [objectiveModels{objective},obj.ObjectiveTheta{objective}]= ...
                    radia.optuna.internal.GaussianProcess.fit( ...
                    x,y(:,objective),categorical, ...
                    obj.DeterministicObjective,obj.ObjectiveTheta{objective});
            end
            [probability,constraintModels]=obj.probabilityFeasible( ...
                candidates,x,constraints(finished,:), ...
                constraintPresent(finished),categorical);
            obj.ConstraintTheta=cell(numel(constraintModels),1);
            for index=1:numel(constraintModels)
                obj.ConstraintTheta{index}=constraintModels{index}.theta;
            end
            feasible=obj.feasibleMask(study,finishedNumbers);
            if size(y,2)==1
                acquisition=obj.expectedImprovement( ...
                    objectiveModels{1},candidates,y(:,1), ...
                    study.Directions(1),feasible,probability);
                running=observations(pending,:);
                acquisition=acquisition.*obj.pendingPenalty( ...
                    candidates,running,categorical);
                [refined,refinedAcquisition]=obj.refineSingleAcquisition( ...
                    candidates,acquisition,objectiveModels{1},y(:,1), ...
                    study.Directions(1),feasible,constraintModels, ...
                    running,categorical,searchSpace);
                candidates=[candidates;refined];
                acquisition=[acquisition;refinedAcquisition];
            else
                acquisition=obj.expectedHypervolumeImprovement( ...
                    objectiveModels,candidates,y,study.Directions, ...
                    feasible,probability);
                running=observations(pending,:);
                acquisition=acquisition.*obj.pendingPenalty( ...
                    candidates,running,categorical);
            end
            if all(~isfinite(acquisition)) || all(acquisition<=0)
                uncertainty=zeros(size(candidates,1),1);
                for objective=1:numel(objectiveModels)
                    [~,standardDeviation]=radia.optuna.internal. ...
                        GaussianProcess.predict( ...
                        objectiveModels{objective},candidates);
                    uncertainty=uncertainty+standardDeviation;
                end
                acquisition=probability.*uncertainty.* ...
                    obj.pendingPenalty(candidates,running,categorical);
            end
            [~,best]=max(acquisition);
            values=obj.decodePoint(candidates(best,:),searchSpace);
        end

        function candidates=candidatePool( ...
                obj,observations,categorical,searchSpace)
            count=obj.CandidateCount;
            dimensions=size(observations,2);
            % Match Optuna's 2048-point QMC preliminary search where the
            % native Sobol table covers the requested dimensionality.
            if dimensions<=32
                qmcSeed=randi(obj.Stream,2^30)-1;
                qmc=radia.optuna.QMCSampler(QMCType="sobol", ...
                    Scramble=true,Seed=qmcSeed);
                candidates=qmc.unitPoints(dimensions,count);
            else
                % Explicit bounded fallback above the native direction
                % table: randomized Latin hypercube, not iid sampling.
                candidates=zeros(count,dimensions);
                for dimension=1:dimensions
                    order=randperm(obj.Stream,count)';
                    candidates(:,dimension)=((order-1)+ ...
                        rand(obj.Stream,count,1))/count;
                end
            end
            candidates=obj.quantizeCategorical( ...
                candidates,categorical,searchSpace);
            if ~isempty(observations)
                candidates=[candidates;observations];
            end
            candidates=unique(candidates,"rows","stable");
        end

        function [refined,acquisition]=refineSingleAcquisition( ...
                obj,candidates,initialAcquisition,model,values,direction, ...
                feasible,constraintModels,pending,categorical,searchSpace)
            refined=zeros(0,size(candidates,2));
            acquisition=zeros(0,1);
            continuous=find(~categorical);
            if obj.LocalSearchCount==0 || isempty(continuous)
                return
            end
            finite=initialAcquisition;
            finite(~isfinite(finite))=-Inf;
            [~,order]=sort(finite,"descend");
            count=min(obj.LocalSearchCount,numel(order));
            anchors=candidates(order(1:count),:);
            anchorAcquisition=initialAcquisition(order(1:count));
            steps=[0.2 0.08 0.032 0.0128 0.00512];
            for step=steps
                probes=zeros(count*(1+2*numel(continuous)),size(anchors,2));
                owners=zeros(size(probes,1),1);
                cursor=1;
                for anchorIndex=1:count
                    probes(cursor,:)=anchors(anchorIndex,:);
                    owners(cursor)=anchorIndex;
                    cursor=cursor+1;
                    for dimension=reshape(continuous,1,[])
                        lower=anchors(anchorIndex,:);
                        upper=anchors(anchorIndex,:);
                        lower(dimension)=max(0,lower(dimension)-step);
                        upper(dimension)=min(1,upper(dimension)+step);
                        probes(cursor,:)=lower;
                        owners(cursor)=anchorIndex;
                        probes(cursor+1,:)=upper;
                        owners(cursor+1)=anchorIndex;
                        cursor=cursor+2;
                    end
                end
                probes=obj.quantizeCategorical(probes,categorical,searchSpace);
                probabilities=obj.probabilityFromModels( ...
                    constraintModels,probes);
                probeAcquisition=obj.expectedImprovement( ...
                    model,probes,values,direction,feasible,probabilities).* ...
                    obj.pendingPenalty(probes,pending,categorical);
                for anchorIndex=1:count
                    local=find(owners==anchorIndex);
                    [bestValue,bestOffset]=max(probeAcquisition(local));
                    if bestValue>anchorAcquisition(anchorIndex)
                        anchors(anchorIndex,:)=probes(local(bestOffset),:);
                        anchorAcquisition(anchorIndex)=bestValue;
                    end
                end
            end
            refined=anchors;
            acquisition=anchorAcquisition;
        end

        function [probability,models]=probabilityFeasible( ...
                obj,candidates,x,constraints,present,categorical)
            probability=ones(size(candidates,1),1);
            models=cell(0,1);
            if isempty(constraints) || ~any(present)
                return
            end
            counts=sum(isfinite(constraints),2);
            expected=max(counts(present));
            if any(counts(present)~=expected)
                error("radia:optuna:ConstraintShape", ...
                    "Trials with different numbers of constraints cannot be compared.");
            end
            models=cell(expected,1);
            if numel(obj.ConstraintTheta)~=expected
                obj.ConstraintTheta=cell(expected,1);
            end
            for constraint=1:expected
                rows=present & isfinite(constraints(:,constraint));
                [models{constraint},theta]=radia.optuna.internal. ...
                    GaussianProcess.fit(x(rows,:),constraints(rows,constraint), ...
                    categorical,obj.DeterministicObjective, ...
                    obj.ConstraintTheta{constraint});
                models{constraint}.theta=theta;
                [meanValue,stdValue]=radia.optuna.internal. ...
                    GaussianProcess.predict(models{constraint},candidates);
                z=(0-meanValue)./max(stdValue,1e-12);
                probability=probability.*(0.5*erfc(-z/sqrt(2)));
            end
            probability=max(probability,realmin("double"));
        end

        function probability=probabilityFromModels(~,models,candidates)
            probability=ones(size(candidates,1),1);
            for constraint=1:numel(models)
                [meanValue,stdValue]=radia.optuna.internal. ...
                    GaussianProcess.predict(models{constraint},candidates);
                z=(0-meanValue)./max(stdValue,1e-12);
                probability=probability.*(0.5*erfc(-z/sqrt(2)));
            end
            probability=max(probability,realmin("double"));
        end

        function acquisition=expectedImprovement(~,model,candidates,values, ...
                direction,feasible,probability)
            [meanValue,stdValue]=radia.optuna.internal. ...
                GaussianProcess.predict(model,candidates);
            if string(direction)=="minimize"
                meanValue=-meanValue;
                signedValues=-values;
            else
                signedValues=values;
            end
            if any(feasible)
                best=max(signedValues(feasible));
                improvement=meanValue-best;
                z=improvement./max(stdValue,1e-12);
                expected=improvement.*(0.5*erfc(-z/sqrt(2)))+ ...
                    stdValue.*exp(-0.5*z.^2)/sqrt(2*pi);
                expected(stdValue<1e-12)=max(improvement(stdValue<1e-12),0);
            else
                expected=ones(size(meanValue));
            end
            acquisition=max(expected,0).*probability;
        end

        function acquisition=expectedHypervolumeImprovement(obj,models, ...
                candidates,values,directions,feasible,probability)
            if ~any(feasible)
                acquisition=probability;
                return
            end
            signs=ones(1,numel(directions));
            signs(string(directions)=="maximize")=-1;
            front=values(feasible,:).*signs;
            worst=max(front,[],1);
            span=max(front,[],1)-min(front,[],1);
            reference=worst+max(0.1*max(span,abs(worst)),1e-9);
            means=zeros(size(candidates,1),numel(models));
            deviations=zeros(size(means));
            for objective=1:numel(models)
                [means(:,objective),deviations(:,objective)]= ...
                    radia.optuna.internal.GaussianProcess.predict( ...
                    models{objective},candidates);
            end
            if numel(models)==2
                uniforms=min(max(rand(obj.Stream,obj.MonteCarloSamples, ...
                    numel(models)),eps),1-eps);
                normals=sqrt(2)*erfinv(2*uniforms-1);
                first=means(:,1)+deviations(:,1).*normals(:,1)';
                second=means(:,2)+deviations(:,2).*normals(:,2)';
                samples=[first(:),second(:)].*signs;
                improvement=radia.optuna.internal.ParetoSupport. ...
                    hypervolumeImprovement2D( ...
                    samples,front,reference);
                acquisition=mean(reshape(improvement, ...
                    size(candidates,1),obj.MonteCarloSamples),2).*probability;
                return
            end
            % For three or more objectives, integrate the probability of
            % dominating objective-space Sobol nodes.  This is the same
            % EHVI integral, evaluated in a vectorized bounded domain, and
            % avoids one recursive hypervolume solve per posterior draw.
            signedMeans=means.*signs;
            lower=min([front;signedMeans-4*deviations],[],1);
            span=max(reference-lower,1e-12);
            qmcSeed=randi(obj.Stream,2^30)-1;
            qmc=radia.optuna.QMCSampler(QMCType="sobol", ...
                Scramble=true,Seed=qmcSeed);
            integration=lower+qmc.unitPoints( ...
                numel(models),obj.MonteCarloSamples).*span;
            dominated=false(size(integration,1),1);
            for point=1:size(front,1)
                dominated=dominated | all( ...
                    integration>=front(point,:),2);
            end
            integration=integration(~dominated,:);
            if isempty(integration)
                acquisition=zeros(size(candidates,1),1);
            else
                dominanceProbability=ones( ...
                    size(candidates,1),size(integration,1));
                for objective=1:numel(models)
                    z=(integration(:,objective)'- ...
                        signedMeans(:,objective))./ ...
                        max(deviations(:,objective),1e-12);
                    dominanceProbability=dominanceProbability.* ...
                        (0.5*erfc(-z/sqrt(2)));
                end
                acquisition=prod(span)*mean(dominanceProbability,2);
            end
            acquisition=acquisition.*probability;
        end

        function penalty=pendingPenalty(~,candidates,pending,categorical)
            penalty=ones(size(candidates,1),1);
            if isempty(pending), return, end
            for candidate=1:size(candidates,1)
                difference=candidates(candidate,:)-pending;
                difference(:,categorical)= ...
                    difference(:,categorical)~=0;
                distance=min(sum(difference.^2,2));
                % Kriging-believer in Optuna excludes the immediate basin
                % around each pending point.  This bounded surrogate keeps
                % the same anti-duplicate behavior without refitting the GP.
                exclusionRadiusSquared=0.01;
                if distance<=exclusionRadiusSquared
                    penalty(candidate)=1e-12;
                else
                    penalty(candidate)=max(1-exp( ...
                        -(distance-exclusionRadiusSquared)/0.04),1e-12);
                end
            end
        end

        function feasible=feasibleMask(~,study,trialNumbers)
            feasible=true(numel(trialNumbers),1);
            if ~study.hasConstraintRecords(), return, end
            feasible(:)=false;
            expected=NaN;
            for index=1:numel(trialNumbers)
                [present,values]=study.constraintRecord(trialNumbers(index));
                if ~present, continue, end
                if isnan(expected), expected=numel(values); end
                if numel(values)~=expected
                    error("radia:optuna:ConstraintShape", ...
                        "Trials with different numbers of constraints cannot be compared.");
                end
                feasible(index)=all(values<=0);
            end
        end

        function [x,y,numbers,pending,constraints,present]= ...
                observations(obj,study,searchSpace,excludedNumber)
            states=study.TrialTable.State;
            rows=find(states=="COMPLETE" | states=="RUNNING");
            rows=rows(study.TrialTable.TrialNumber(rows)~=excludedNumber);
            numbers=study.TrialTable.TrialNumber(rows);
            pending=states(rows)=="RUNNING";
            x=zeros(numel(rows),numel(searchSpace));
            y=NaN(numel(rows),numel(study.Directions));
            constraintCells=cell(numel(rows),1);
            present=false(numel(rows),1);
            maximumConstraints=0;
            for index=1:numel(rows)
                x(index,:)=obj.encodeTrial(study,numbers(index),searchSpace);
                if ~pending(index)
                    for objective=1:numel(study.Directions)
                        mask=study.ObjectiveTable.TrialNumber==numbers(index) & ...
                            study.ObjectiveTable.ObjectiveIndex==objective;
                        if sum(mask)==1
                            y(index,objective)=study.ObjectiveTable.Value(mask);
                        end
                    end
                    [present(index),constraintCells{index}]= ...
                        study.constraintRecord(numbers(index));
                    maximumConstraints=max(maximumConstraints, ...
                        numel(constraintCells{index}));
                end
            end
            constraints=NaN(numel(rows),maximumConstraints);
            for index=1:numel(rows)
                if present(index) && ~isempty(constraintCells{index})
                    constraints(index,1:numel(constraintCells{index}))= ...
                        constraintCells{index};
                end
            end
            completeRows=~pending & all(isfinite(y),2);
            keep=pending | completeRows;
            x=x(keep,:); y=y(keep,:); numbers=numbers(keep);
            pending=pending(keep); constraints=constraints(keep,:);
            present=present(keep);
        end

        function encoded=encodeTrial(obj,study,trialNumber,searchSpace)
            encoded=zeros(1,numel(searchSpace));
            for dimension=1:numel(searchSpace)
                row=study.ParamTable.TrialNumber==trialNumber & ...
                    study.ParamTable.Name==searchSpace(dimension).name;
                if sum(row)~=1
                    error("radia:optuna:GPObservations", ...
                        "A GP intersection parameter is missing from a trial.");
                end
                if isfinite(study.ParamTable.ValueNumeric(row))
                    value=study.ParamTable.ValueNumeric(row);
                else
                    value=jsondecode(study.ParamTable.ValueText(row));
                end
                encoded(dimension)=obj.encodeValue( ...
                    value,searchSpace(dimension).distribution);
            end
        end

        function encoded=encodeValue(~,value,distribution)
            if distribution.kind=="categorical"
                tokens=radia.optuna.internal.DistributionCodec. ...
                    choiceTokens(distribution.choices);
                token=radia.optuna.internal.DistributionCodec.choiceToken(value);
                index=find(tokens==token,1);
                if isempty(index)
                    error("radia:optuna:GPObservations", ...
                        "Observed categorical value is outside its distribution.");
                end
                encoded=(index-1)/max(1,numel(tokens)-1);
            elseif distribution.log
                encoded=(log(double(value))-log(distribution.low))/ ...
                    (log(distribution.high)-log(distribution.low));
            else
                encoded=(double(value)-distribution.low)/ ...
                    (distribution.high-distribution.low);
            end
            encoded=min(max(encoded,0),1);
        end

        function values=decodePoint(~,point,searchSpace)
            values=cell(1,numel(searchSpace));
            for dimension=1:numel(searchSpace)
                distribution=searchSpace(dimension).distribution;
                unit=min(max(point(dimension),0),1);
                if distribution.kind=="categorical"
                    count=numel(distribution.choices);
                    index=min(max(round(unit*max(1,count-1))+1,1),count);
                    values{dimension}=radia.optuna.internal. ...
                        DistributionCodec.choiceAt(distribution.choices,index);
                else
                    if distribution.log
                        value=exp(log(distribution.low)+unit* ...
                            (log(distribution.high)-log(distribution.low)));
                    else
                        value=distribution.low+unit* ...
                            (distribution.high-distribution.low);
                    end
                    if isfinite(distribution.step)
                        value=distribution.low+round( ...
                            (value-distribution.low)/distribution.step)* ...
                            distribution.step;
                    end
                    value=min(max(value,distribution.low),distribution.high);
                    if distribution.kind=="integer", value=round(value); end
                    values{dimension}=value;
                end
            end
        end

        function mask=categoricalMask(~,searchSpace)
            mask=false(1,numel(searchSpace));
            for index=1:numel(searchSpace)
                mask(index)=searchSpace(index).distribution.kind=="categorical";
            end
        end

        function candidates=quantizeCategorical( ...
                ~,candidates,categorical,searchSpace)
            for dimension=find(categorical)
                levelCount=numel(searchSpace(dimension).distribution.choices);
                levels=linspace(0,1,levelCount)';
                indices=min(floor(candidates(:,dimension)*numel(levels))+1, ...
                    numel(levels));
                candidates(:,dimension)=levels(indices);
            end
        end

        function value=randomValue(obj,distribution)
            if distribution.kind=="categorical"
                index=randi(obj.Stream,numel(distribution.choices));
                value=radia.optuna.internal.DistributionCodec. ...
                    choiceAt(distribution.choices,index);
                return
            end
            unit=rand(obj.Stream);
            if distribution.log
                value=exp(log(distribution.low)+unit* ...
                    (log(distribution.high)-log(distribution.low)));
            else
                value=distribution.low+unit* ...
                    (distribution.high-distribution.low);
            end
            if isfinite(distribution.step)
                value=distribution.low+round( ...
                    (value-distribution.low)/distribution.step)*distribution.step;
            end
            value=min(max(value,distribution.low),distribution.high);
            if distribution.kind=="integer", value=round(value); end
        end

        function attach(obj,study)
            changed=isempty(obj.AttachedStudy) || ~isequal(obj.AttachedStudy,study);
            if changed
                obj.AttachedStudy=study;
                obj.Stream=radia.optuna.internal.NumpyRandomState(obj.Seed);
                obj.ObjectiveTheta=cell(0,1);
                obj.ConstraintTheta=cell(0,1);
                obj.PythonOptuna=[];
                obj.PythonStudy=[];
                obj.PythonTrial=[];
                obj.PythonTrialNumber=NaN;
                obj.Restored=false;
            end
            if obj.Backend=="upstream-python"
                obj.Restored=true;
                return
            end
            if obj.Restored, return, end
            state=study.samplerState(obj.SamplerName,obj.StateSchema);
            if ~isempty(state), obj.restoreState(state); end
            obj.Restored=true;
        end

        function restoreState(obj,state)
            required=["schema","seed","random_state", ...
                "objective_theta","constraint_theta"];
            if ~isstruct(state) || ~isscalar(state) || ...
                    any(~isfield(state,required)) || ...
                    string(state.schema)~=obj.StateSchema || ...
                    double(state.seed)~=obj.Seed
                error("radia:optuna:GPState", ...
                    "Stored GP sampler state is invalid or incompatible.");
            end
            try
                obj.Stream.State=state.random_state;
            catch exception
                error("radia:optuna:GPState", ...
                    "Stored GP random state is invalid: %s",exception.message);
            end
            obj.ObjectiveTheta=state.objective_theta;
            obj.ConstraintTheta=state.constraint_theta;
        end

        function recordState(obj,study,trialNumber)
            state=struct("schema",obj.StateSchema,"seed",obj.Seed, ...
                "random_state",obj.Stream.State, ...
                "objective_theta",{obj.ObjectiveTheta}, ...
                "constraint_theta",{obj.ConstraintTheta});
            generation=sum(study.TrialTable.State=="COMPLETE");
            study.recordSamplerState(obj.SamplerName,obj.StateSchema, ...
                trialNumber,generation,state);
        end
    end
end
