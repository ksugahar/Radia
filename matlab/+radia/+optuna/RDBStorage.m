classdef RDBStorage < radia.optuna.BaseStorage
    %RDBSTORAGE Optuna 4.9 relational storage through its Python backend.
    %   SQLAlchemy schema ownership remains with pinned optuna==4.9.0. The
    %   MATLAB class converts only documented storage values and snapshots;
    %   it never substitutes MAT-file storage for a relational database.

    properties (SetAccess=immutable)
        url (1,1) string
    end

    properties (Access=protected)
        PythonStorage
        OptunaModule
        StoragesModule
        StudyModule
        TrialModule
        DistributionsModule
        JsonModule
        HeartbeatStaleTrialCallback = []
    end

    methods
        function obj=RDBStorage(url,options)
            arguments
                url (1,1) string
                options.engine_kwargs (1,1) struct = struct()
                options.skip_compatibility_check (1,1) logical = false
                options.heartbeat_interval double = NaN
                options.grace_period double = NaN
                options.heartbeat_stale_trial_callback = []
                options.failed_trial_callback = []
                options.skip_table_creation (1,1) logical = false
                options.InternalPythonStorage = []
            end
            obj.requireUpstream();
            obj.url=url;
            if ~isempty(options.InternalPythonStorage)
                obj.PythonStorage=options.InternalPythonStorage;
                return
            end
            callback=options.heartbeat_stale_trial_callback;
            if isempty(callback) && ~isempty(options.failed_trial_callback)
                warning("radia:optuna:FutureWarning", ...
                    "failed_trial_callback is deprecated; use " + ...
                    "heartbeat_stale_trial_callback instead.");
                callback=options.failed_trial_callback;
            end
            obj.HeartbeatStaleTrialCallback=callback;
            names={"skip_compatibility_check", ...
                "skip_table_creation"};
            values={options.skip_compatibility_check, ...
                options.skip_table_creation};
            if ~isempty(fieldnames(options.engine_kwargs))
                names{end+1}="engine_kwargs";
                values{end+1}=obj.toPythonValue(options.engine_kwargs);
            end
            if ~isnan(options.heartbeat_interval)
                names{end+1}="heartbeat_interval";
                values{end+1}=int64(options.heartbeat_interval);
            end
            if ~isnan(options.grace_period)
                names{end+1}="grace_period";
                values{end+1}=int64(options.grace_period);
            end
            argumentsList=cell(1,2*numel(names));
            argumentsList(1:2:end)=names;
            argumentsList(2:2:end)=values;
            obj.PythonStorage=obj.StoragesModule.RDBStorage( ...
                char(url),pyargs(argumentsList{:}));
        end

        function study_id=create_new_study(obj,directions,study_name)
            if nargin<3 || strlength(string(study_name))==0
                study_name=[];
            end
            pythonDirections=obj.pythonDirections(directions);
            try
                if isempty(study_name)
                    value=obj.PythonStorage.create_new_study(pythonDirections);
                else
                    value=obj.PythonStorage.create_new_study( ...
                        pythonDirections,pyargs("study_name",char(study_name)));
                end
            catch cause
                obj.rethrowMapped(cause);
            end
            study_id=double(value);
        end

        function trial_id=create_new_trial(obj,study_id,template_trial)
            try
                if nargin<3 || isempty(template_trial)
                    value=obj.PythonStorage.create_new_trial(int64(study_id));
                else
                    value=obj.PythonStorage.create_new_trial( ...
                        int64(study_id),obj.pythonTrial(template_trial));
                end
            catch cause
                obj.rethrowMapped(cause);
            end
            trial_id=double(value);
        end

        function delete_study(obj,study_id)
            try
                obj.PythonStorage.delete_study(int64(study_id));
            catch cause
                obj.rethrowMapped(cause);
            end
        end

        function studies=get_all_studies(obj)
            source=cell(obj.PythonStorage.get_all_studies());
            studies=radia.optuna.StudySummary.empty(0,1);
            for index=1:numel(source)
                studies(index,1)=obj.studySummary(source{index}); %#ok<AGROW>
            end
        end

        function trials=get_all_trials(obj,study_id,deepcopy,states)
            if nargin<3, deepcopy=true; end
            if nargin<4 || isempty(states)
                source=obj.PythonStorage.get_all_trials( ...
                    int64(study_id),pyargs("deepcopy",logical(deepcopy)));
            else
                source=obj.PythonStorage.get_all_trials(int64(study_id), ...
                    pyargs("deepcopy",logical(deepcopy), ...
                    "states",obj.pythonTrialStates(states)));
            end
            trials=obj.frozenTrials(source);
        end

        function trial=get_best_trial(obj,study_id)
            trial=obj.frozenTrial( ...
                obj.PythonStorage.get_best_trial(int64(study_id)));
        end

        function count=get_n_trials(obj,study_id,state)
            if nargin<3 || isempty(state)
                value=obj.PythonStorage.get_n_trials(int64(study_id));
            else
                value=obj.PythonStorage.get_n_trials(int64(study_id), ...
                    pyargs("state",obj.pythonTrialState(state)));
            end
            count=double(value);
        end

        function directions=get_study_directions(obj,study_id)
            source=cell(obj.PythonStorage.get_study_directions(int64(study_id)));
            directions(1,numel(source))=radia.optuna.StudyDirection.MINIMIZE;
            for index=1:numel(source)
                directions(index)=radia.optuna.StudyDirection.from( ...
                    obj.pythonEnumName(source{index}));
            end
        end

        function study_id=get_study_id_from_name(obj,study_name)
            try
                study_id=double(obj.PythonStorage.get_study_id_from_name( ...
                    char(study_name)));
            catch cause
                obj.rethrowMapped(cause);
            end
        end

        function study_name=get_study_name_from_id(obj,study_id)
            try
                study_name=string(obj.PythonStorage.get_study_name_from_id( ...
                    int64(study_id)));
            catch cause
                obj.rethrowMapped(cause);
            end
        end

        function attributes=get_study_system_attrs(obj,study_id)
            attributes=obj.fromPythonJson( ...
                obj.PythonStorage.get_study_system_attrs(int64(study_id)));
        end

        function attributes=get_study_user_attrs(obj,study_id)
            attributes=obj.fromPythonJson( ...
                obj.PythonStorage.get_study_user_attrs(int64(study_id)));
        end

        function trial=get_trial(obj,trial_id)
            try
                trial=obj.frozenTrial( ...
                    obj.PythonStorage.get_trial(int64(trial_id)));
            catch cause
                obj.rethrowMapped(cause);
            end
        end

        function trial_id=get_trial_id_from_study_id_trial_number( ...
                obj,study_id,trial_number)
            trial_id=double( ...
                obj.PythonStorage.get_trial_id_from_study_id_trial_number( ...
                int64(study_id),int64(trial_number)));
        end

        function trial_number=get_trial_number_from_id(obj,trial_id)
            trial_number=double( ...
                obj.PythonStorage.get_trial_number_from_id(int64(trial_id)));
        end

        function value=get_trial_param(obj,trial_id,param_name)
            value=double(obj.PythonStorage.get_trial_param( ...
                int64(trial_id),char(param_name)));
        end

        function params=get_trial_params(obj,trial_id)
            params=obj.fromPythonJson( ...
                obj.PythonStorage.get_trial_params(int64(trial_id)));
        end

        function attributes=get_trial_system_attrs(obj,trial_id)
            attributes=obj.fromPythonJson( ...
                obj.PythonStorage.get_trial_system_attrs(int64(trial_id)));
        end

        function attributes=get_trial_user_attrs(obj,trial_id)
            attributes=obj.fromPythonJson( ...
                obj.PythonStorage.get_trial_user_attrs(int64(trial_id)));
        end

        function remove_session(obj)
            obj.PythonStorage.remove_session();
        end

        function set_study_system_attr(obj,study_id,key,value)
            obj.PythonStorage.set_study_system_attr( ...
                int64(study_id),char(key),obj.toPythonValue(value));
        end

        function set_study_user_attr(obj,study_id,key,value)
            obj.PythonStorage.set_study_user_attr( ...
                int64(study_id),char(key),obj.toPythonValue(value));
        end

        function set_trial_intermediate_value(obj,trial_id,step,value)
            obj.ensureTrialUpdatable(trial_id);
            obj.PythonStorage.set_trial_intermediate_value( ...
                int64(trial_id),int64(step),double(value));
        end

        function set_trial_param( ...
                obj,trial_id,param_name,param_value_internal,distribution)
            obj.ensureTrialUpdatable(trial_id);
            pythonDistribution= ...
                obj.DistributionsModule.json_to_distribution( ...
                char(radia.optuna.distribution_to_json(distribution)));
            obj.PythonStorage.set_trial_param(int64(trial_id), ...
                char(param_name),double(param_value_internal), ...
                pythonDistribution);
        end

        function updated=set_trial_state_values(obj,trial_id,state,values)
            if nargin<4 || isempty(values)
                result=obj.PythonStorage.set_trial_state_values( ...
                    int64(trial_id),obj.pythonTrialState(state));
            else
                result=obj.PythonStorage.set_trial_state_values( ...
                    int64(trial_id),obj.pythonTrialState(state), ...
                    py.list(num2cell(reshape(double(values),1,[]))));
            end
            updated=logical(result);
        end

        function set_trial_system_attr(obj,trial_id,key,value)
            obj.ensureTrialUpdatable(trial_id);
            obj.PythonStorage.set_trial_system_attr( ...
                int64(trial_id),char(key),obj.toPythonValue(value));
        end

        function set_trial_user_attr(obj,trial_id,key,value)
            obj.ensureTrialUpdatable(trial_id);
            obj.PythonStorage.set_trial_user_attr( ...
                int64(trial_id),char(key),obj.toPythonValue(value));
        end

        function check_trial_is_updatable(~,~,trial_state)
            state=radia.optuna.TrialState.toStorage(trial_state);
            if ismember(state,["COMPLETE","PRUNED","FAIL"])
                throw(radia.optuna.UpdateFinishedTrialError( ...
                    "Trial is already finished and cannot be updated."));
            end
        end

        function record_heartbeat(obj,trial_id)
            obj.PythonStorage.record_heartbeat(int64(trial_id));
        end

        function interval=get_heartbeat_interval(obj)
            interval=obj.optionalPythonNumber( ...
                obj.PythonStorage.get_heartbeat_interval());
        end

        function callback=get_heartbeat_stale_trial_callback(obj)
            callback=obj.HeartbeatStaleTrialCallback;
        end

        function callback=get_failed_trial_callback(obj)
            warning("radia:optuna:FutureWarning", ...
                "get_failed_trial_callback has been deprecated in v4.9.0. " + ...
                "Use get_heartbeat_stale_trial_callback instead.");
            callback=obj.get_heartbeat_stale_trial_callback();
        end

        function version=get_current_version(obj)
            version=string(obj.PythonStorage.get_current_version());
        end

        function version=get_head_version(obj)
            version=string(obj.PythonStorage.get_head_version());
        end

        function versions=get_all_versions(obj)
            versions=string(cell(obj.PythonStorage.get_all_versions()));
        end

        function upgrade(obj)
            obj.PythonStorage.upgrade();
        end
    end

    methods (Hidden=true)
        function storage=pythonStorageHandle(obj)
            %PYTHONSTORAGEHANDLE Internal bridge for gRPC server ownership.
            storage=obj.PythonStorage;
        end

        function dispose(obj)
            %DISPOSE Release SQLAlchemy connections owned by this adapter.
            if ismethod(obj.PythonStorage,"remove_session")
                obj.PythonStorage.remove_session();
            end
            try
                obj.PythonStorage.engine.dispose();
            catch
                % Non-RDB subclasses such as gRPC proxies have no engine.
            end
        end
    end

    methods (Access=private)
        function ensureTrialUpdatable(obj,trial_id)
            trial=obj.get_trial(trial_id);
            obj.check_trial_is_updatable(trial_id,trial.State);
        end

        function value=pythonTrialState(obj,state)
            name=radia.optuna.TrialState.toStorage(state);
            stateClass=py.builtins.getattr(obj.TrialModule,"TrialState");
            value=py.builtins.getattr( ...
                stateClass,char(name));
        end

        function values=pythonTrialStates(obj,states)
            names=reshape(radia.optuna.TrialState.toStorage(states),1,[]);
            entries=cell(1,numel(names));
            for index=1:numel(names)
                entries{index}=obj.pythonTrialState(names(index));
            end
            values=py.tuple(entries);
        end

        function values=pythonDirections(obj,directions)
            source=reshape(radia.optuna.StudyDirection.from(directions),1,[]);
            directionClass=py.builtins.getattr( ...
                obj.StudyModule,"StudyDirection");
            entries=cell(1,numel(source));
            for index=1:numel(source)
                entries{index}=py.builtins.getattr( ...
                    directionClass, ...
                    char(string(source(index))));
            end
            values=py.list(entries);
        end

        function trial=pythonTrial(obj,source)
            if ~isa(source,"radia.optuna.FrozenTrial") || ~isscalar(source)
                error("radia:optuna:StorageTrial", ...
                    "template_trial must be one FrozenTrial.");
            end
            params=obj.toPythonValue(source.Params);
            distributions=py.dict;
            names=string(fieldnames(source.Distributions));
            for index=1:numel(names)
                spec=source.Distributions.(names(index));
                distributions{char(names(index))}= ...
                    obj.DistributionsModule.json_to_distribution( ...
                    char(radia.optuna.distribution_to_json(spec)));
            end
            userAttrs=obj.toPythonValue(source.UserAttrs);
            systemAttrs=obj.toPythonValue(source.SystemAttrs);
            intermediate=py.dict;
            for index=1:height(source.IntermediateValues)
                intermediate{int64(source.IntermediateValues.Step(index))}= ...
                    double(source.IntermediateValues.Value(index));
            end
            keyword={"state",obj.pythonTrialState(source.State), ...
                "params",params,"distributions",distributions, ...
                "intermediate_values",intermediate, ...
                "user_attrs",userAttrs,"system_attrs",systemAttrs};
            values=reshape(double(source.Values),1,[]);
            if numel(values)==1
                keyword(end+1:end+2)={"value",values(1)};
            elseif ~all(isnan(values))
                keyword(end+1:end+2)={"values",py.list(num2cell(values))};
            end
            trial=obj.TrialModule.create_trial(pyargs(keyword{:}));
        end

        function trials=frozenTrials(obj,source)
            values=cell(source);
            if isempty(values)
                trials=radia.optuna.FrozenTrial.empty(0,1);
                return
            end
            trials(numel(values),1)=obj.frozenTrial(values{end});
            for index=1:numel(values)
                trials(index,1)=obj.frozenTrial(values{index});
            end
        end

        function trial=frozenTrial(obj,source)
            params=obj.fromPythonJson(source.params);
            userAttrs=obj.fromPythonJson(source.user_attrs);
            systemAttrs=obj.fromPythonJson(source.system_attrs);
            distributions=struct();
            names=cell(py.list(source.distributions.keys()));
            for index=1:numel(names)
                name=string(names{index});
                encoded=string( ...
                    obj.DistributionsModule.distribution_to_json( ...
                    source.distributions{char(name)}));
                distribution=radia.optuna.json_to_distribution(encoded);
                distributions.(matlab.lang.makeValidName(name))= ...
                    distribution.toStruct();
            end
            intermediate=radia.optuna.Trial.emptyIntermediateTable();
            payload=obj.fromPythonJson(source.intermediate_values);
            fields=string(fieldnames(payload));
            for index=1:numel(fields)
                token=erase(fields(index),"x");
                step=str2double(token);
                intermediate(end+1,:)={step,double(payload.(fields(index))),NaT}; %#ok<AGROW>
            end
            if ~isempty(intermediate)
                intermediate=sortrows(intermediate,"Step");
            end
            values=obj.fromPythonJson(source.values);
            if isempty(values), values=NaN; end
            trial=radia.optuna.FrozenTrial( ...
                Number=double(source.number), ...
                State=obj.pythonEnumName(source.state), ...
                Values=reshape(double(values),1,[]),Params=params, ...
                Distributions=distributions,IntermediateValues=intermediate, ...
                UserAttrs=userAttrs,SystemAttrs=systemAttrs, ...
                DatetimeStart=obj.pythonDatetime(source.datetime_start), ...
                DatetimeComplete=obj.pythonDatetime(source.datetime_complete));
        end

        function summary=studySummary(obj,source)
            directions=cell(source.directions);
            names=strings(1,numel(directions));
            for index=1:numel(directions)
                names(index)=obj.pythonEnumName(directions{index});
            end
            best=[];
            try
                best=obj.frozenTrial(source.best_trial);
            catch
                % Multi-objective and empty studies have no best trial.
            end
            summary=radia.optuna.StudySummary(string(source.study_name),[], ...
                best,obj.fromPythonJson(source.user_attrs), ...
                obj.fromPythonJson(source.system_attrs),0,NaT, ...
                double(py.builtins.getattr(source,"_study_id")), ...
                directions=names);
        end

        function value=fromPythonJson(obj,source)
            dumps=py.builtins.getattr(obj.JsonModule,"dumps");
            text=string(dumps(source,pyargs("sort_keys",true)));
            if text=="null"
                value=[];
            else
                value=jsondecode(char(text));
            end
        end

        function value=toPythonValue(obj,source)
            loads=py.builtins.getattr(obj.JsonModule,"loads");
            value=loads(char(jsonencode(source)));
        end

        function value=pythonDatetime(~,source)
            value=NaT;
            if isa(source,"py.NoneType"), return, end
            text=string(source);
            try
                value=datetime(text,"InputFormat","yyyy-MM-dd HH:mm:ss.SSSSSS");
            catch
                try
                    value=datetime(text,"InputFormat","yyyy-MM-dd HH:mm:ss");
                catch
                    value=NaT;
                end
            end
        end

        function value=optionalPythonNumber(~,source)
            if isa(source,"py.NoneType")
                value=[];
            else
                value=double(source);
            end
        end

        function name=pythonEnumName(~,source)
            name=string(py.builtins.getattr(source,"name"));
        end

        function rethrowMapped(~,cause)
            message=string(cause.message);
            if contains(message,"DuplicatedStudyError")
                throw(radia.optuna.DuplicatedStudyError(message));
            elseif contains(message,"UpdateFinishedTrialError")
                throw(radia.optuna.UpdateFinishedTrialError(message));
            elseif contains(message,"KeyError") || contains(message,"not found")
                error("radia:optuna:StudyNotFound","%s",message);
            end
            rethrow(cause)
        end

        function requireUpstream(obj)
            environment=pyenv;
            if environment.Status=="NotLoaded"
                environment=pyenv(ExecutionMode="InProcess");
            end
            if environment.ExecutionMode~="InProcess"
                error("radia:optuna:RDBStoragePython", ...
                    "RDBStorage requires in-process Python. Configure pyenv first.");
            end
            try
                obj.OptunaModule=py.importlib.import_module("optuna");
                obj.StoragesModule=py.importlib.import_module( ...
                    "optuna.storages");
                obj.StudyModule=py.importlib.import_module("optuna.study");
                obj.TrialModule=py.importlib.import_module("optuna.trial");
                obj.DistributionsModule=py.importlib.import_module( ...
                    "optuna.distributions");
                obj.JsonModule=py.importlib.import_module("json");
                version=string(py.builtins.getattr( ...
                    obj.OptunaModule,"__version__"));
            catch cause
                error("radia:optuna:RDBStoragePython", ...
                    "Install radia-optuna[upstream]: %s",cause.message);
            end
            if version~="4.9.0"
                error("radia:optuna:RDBStorageVersion", ...
                    "RDBStorage requires optuna==4.9.0, found %s.",version);
            end
        end
    end
end
