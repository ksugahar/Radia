classdef InMemoryStorage < radia.optuna.BaseStorage
    %INMEMORYSTORAGE Process-local Optuna storage with integer IDs.

    properties (Access=private)
        Studies cell = cell(1,0)
        StudyNames string = strings(1,0)
        TrialStudyIds double = zeros(1,0)
        TrialNumbers double = zeros(1,0)
        TrialHandles cell = cell(1,0)
    end

    methods
        function study_id=create_new_study(obj,directions,study_name)
            if nargin<3 || isempty(study_name)
                study_name="";
            end
            study_name=string(study_name);
            if strlength(study_name)>0 && any(obj.StudyNames==study_name)
                throw(radia.optuna.DuplicatedStudyError( ...
                    "Another study with name '"+study_name+"' already exists."));
            end
            study=radia.optuna.Study(Name=study_name,Directions=directions, ...
                AutoSave=false);
            study_id=numel(obj.Studies);
            obj.Studies{study_id+1}=study;
            obj.StudyNames(study_id+1)=study.Name;
        end

        function trial_id=create_new_trial(obj,study_id,template_trial)
            if nargin<3
                template_trial=[];
            end
            study=obj.study(study_id);
            if isempty(template_trial)
                trial=study.ask();
                number=trial.Number;
            else
                if ~isa(template_trial,"radia.optuna.FrozenTrial") || ...
                        ~isscalar(template_trial)
                    error("radia:optuna:StorageTrial", ...
                        "template_trial must be one FrozenTrial.");
                end
                study.add_trial(template_trial);
                rows=study.TrialTable;
                number=rows.TrialNumber(end);
                trial=[];
                if ismember(template_trial.State,["RUNNING","WAITING"])
                    trial=radia.optuna.Trial(study,number);
                    snapshots=study.get_trials();
                    trial.restoreSnapshot(snapshots(end));
                end
            end
            trial_id=numel(obj.TrialStudyIds);
            obj.TrialStudyIds(trial_id+1)=study_id;
            obj.TrialNumbers(trial_id+1)=number;
            obj.TrialHandles{trial_id+1}=trial;
        end

        function delete_study(obj,study_id)
            obj.study(study_id);
            obj.Studies{study_id+1}=[];
            obj.StudyNames(study_id+1)=missing;
            rows=obj.TrialStudyIds==study_id;
            obj.TrialStudyIds(rows)=NaN;
            obj.TrialHandles(rows)={[]};
        end

        function study_id=get_study_id_from_name(obj,study_name)
            index=find(obj.StudyNames==string(study_name),1);
            if isempty(index)
                error("radia:optuna:StudyNotFound", ...
                    "No study with name '%s' exists.",study_name);
            end
            study_id=index-1;
        end

        function study_name=get_study_name_from_id(obj,study_id)
            study_name=obj.study(study_id).Name;
        end

        function directions=get_study_directions(obj,study_id)
            directions=radia.optuna.StudyDirection.from( ...
                obj.study(study_id).Directions);
        end

        function set_study_user_attr(obj,study_id,key,value)
            obj.study(study_id).set_user_attr(string(key),value);
        end

        function set_study_system_attr(obj,study_id,key,value)
            obj.study(study_id).set_system_attr(string(key),value);
        end

        function attributes=get_study_user_attrs(obj,study_id)
            attributes=obj.study(study_id).user_attrs();
        end

        function attributes=get_study_system_attrs(obj,study_id)
            attributes=obj.study(study_id).system_attrs();
        end

        function trial_id=get_trial_id_from_study_id_trial_number( ...
                obj,study_id,trial_number)
            obj.study(study_id);
            index=find(obj.TrialStudyIds==study_id & ...
                obj.TrialNumbers==double(trial_number),1);
            if isempty(index)
                error("radia:optuna:TrialNotFound", ...
                    "Trial %d does not exist in study %d.",trial_number,study_id);
            end
            trial_id=index-1;
        end

        function trial_number=get_trial_number_from_id(obj,trial_id)
            [~,trial_number]=obj.trialLocation(trial_id);
        end

        function trial=get_trial(obj,trial_id)
            [study,number]=obj.trialLocation(trial_id);
            trials=study.get_trials();
            index=find([trials.Number]==number,1);
            trial=trials(index);
        end

        function trials=get_all_trials(obj,study_id,deepcopy,states) %#ok<INUSD>
            if nargin<3, deepcopy=true; end %#ok<NASGU>
            study=obj.study(study_id);
            if nargin<4 || isempty(states)
                trials=study.get_trials();
            else
                trials=study.get_trials(radia.optuna.TrialState.toStorage(states));
            end
        end

        function count=get_n_trials(obj,study_id,state)
            if nargin<3 || isempty(state)
                count=numel(obj.get_all_trials(study_id));
            else
                count=numel(obj.get_all_trials(study_id,true,state));
            end
        end

        function trial=get_best_trial(obj,study_id)
            trial=obj.study(study_id).best_trial();
        end

        function studies=get_all_studies(obj)
            items=cell(1,0);
            for index=1:numel(obj.Studies)
                if isempty(obj.Studies{index}), continue; end
                study=obj.Studies{index};
                trials=study.get_trials();
                complete=trials(string({trials.State})=="COMPLETE");
                best=[];
                if isscalar(study.Directions) && ~isempty(complete)
                    best=study.best_trial();
                end
                started=NaT;
                if ~isempty(trials), started=trials(1).DatetimeStart; end
                directions=radia.optuna.StudyDirection.from(study.Directions);
                if isscalar(directions)
                    summary=radia.optuna.StudySummary(study.Name, ...
                        directions,best,study.user_attrs(), ...
                        study.system_attrs(),numel(trials),started,index-1);
                else
                    summary=radia.optuna.StudySummary(study.Name,[],best, ...
                        study.user_attrs(),study.system_attrs(), ...
                        numel(trials),started,index-1,directions=directions);
                end
                items{end+1}=summary; %#ok<AGROW>
            end
            if isempty(items)
                studies=radia.optuna.StudySummary.empty(0,1);
            else
                studies=reshape([items{:}],[],1);
            end
        end

        function set_trial_param(obj,trial_id,param_name, ...
                param_value_internal,distribution)
            trial=obj.runningTrial(trial_id);
            spec=radia.optuna.internal.DistributionCodec.normalize(distribution);
            if isa(distribution,"radia.optuna.BaseDistribution")
                value=distribution.to_external_repr(param_value_internal);
            elseif spec.kind=="categorical"
                value=radia.optuna.internal.DistributionCodec.choiceAt( ...
                    spec.choices,double(param_value_internal)+1);
            else
                value=double(param_value_internal);
            end
            trial.setStorageParameter(string(param_name),value,spec);
        end

        function value=get_trial_param(obj,trial_id,param_name)
            trial=obj.get_trial(trial_id);
            key=matlab.lang.makeValidName(string(param_name));
            if ~isfield(trial.Params,key)
                error("radia:optuna:ParameterNotFound", ...
                    "Parameter '%s' does not exist.",param_name);
            end
            external=trial.Params.(key);
            distribution=trial.Distributions.(key);
            spec=radia.optuna.internal.DistributionCodec.normalize(distribution);
            if spec.kind=="categorical"
                tokens=radia.optuna.internal.DistributionCodec.choiceTokens(spec.choices);
                value=find(tokens== ...
                    radia.optuna.internal.DistributionCodec.choiceToken(external),1)-1;
            else
                value=double(external);
            end
        end

        function params=get_trial_params(obj,trial_id)
            trial=obj.get_trial(trial_id);
            params=trial.Params;
        end

        function updated=set_trial_state_values(obj,trial_id,state,values)
            if nargin<4, values=[]; end
            [study,~,trial]=obj.trialLocation(trial_id);
            current=obj.get_trial(trial_id).State;
            target=radia.optuna.TrialState.toStorage(state);
            if current=="WAITING" && target=="RUNNING"
                trial=study.ask();
                obj.TrialHandles{trial_id+1}=trial;
                updated=true;
                return
            end
            obj.check_trial_is_updatable(trial_id,current);
            if isempty(trial)
                trial=obj.runningTrial(trial_id);
            end
            study.tell(trial,reshape(double(values),1,[]),State=target);
            obj.TrialHandles{trial_id+1}=[];
            updated=true;
        end

        function set_trial_intermediate_value(obj,trial_id,step,value)
            obj.runningTrial(trial_id).report(double(value),double(step));
        end

        function set_trial_user_attr(obj,trial_id,key,value)
            obj.runningTrial(trial_id).set_user_attr(string(key),value);
        end

        function set_trial_system_attr(obj,trial_id,key,value)
            obj.runningTrial(trial_id).set_system_attr(string(key),value);
        end

        function attributes=get_trial_user_attrs(obj,trial_id)
            trial=obj.get_trial(trial_id);
            attributes=trial.UserAttrs;
        end

        function attributes=get_trial_system_attrs(obj,trial_id)
            trial=obj.get_trial(trial_id);
            attributes=trial.SystemAttrs;
        end

        function check_trial_is_updatable(~,~,trial_state)
            state=radia.optuna.TrialState.from(trial_state);
            if state.is_finished()
                throw(radia.optuna.UpdateFinishedTrialError( ...
                    "Trial has already finished and cannot be updated."));
            end
        end

        function remove_session(~)
            % In-memory storage owns no thread-local database session.
        end
    end

    methods (Access=private)
        function study=study(obj,study_id)
            if ~isscalar(study_id) || study_id~=floor(study_id) || ...
                    study_id<0 || study_id>=numel(obj.Studies) || ...
                    isempty(obj.Studies{study_id+1})
                error("radia:optuna:StudyNotFound", ...
                    "Study %g does not exist.",study_id);
            end
            study=obj.Studies{study_id+1};
        end

        function [study,number,trial]=trialLocation(obj,trial_id)
            if ~isscalar(trial_id) || trial_id~=floor(trial_id) || ...
                    trial_id<0 || trial_id>=numel(obj.TrialStudyIds) || ...
                    ~isfinite(obj.TrialStudyIds(trial_id+1))
                error("radia:optuna:TrialNotFound", ...
                    "Trial %g does not exist.",trial_id);
            end
            study=obj.study(obj.TrialStudyIds(trial_id+1));
            number=obj.TrialNumbers(trial_id+1);
            trial=obj.TrialHandles{trial_id+1};
        end

        function trial=runningTrial(obj,trial_id)
            [~,~,trial]=obj.trialLocation(trial_id);
            current=obj.get_trial(trial_id).State;
            obj.check_trial_is_updatable(trial_id,current);
            if isempty(trial)
                error("radia:optuna:TrialState", ...
                    "Trial %d is not an active RUNNING trial.",trial_id);
            end
        end
    end
end
