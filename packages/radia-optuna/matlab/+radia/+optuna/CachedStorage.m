classdef CachedStorage < radia.optuna.BaseStorage
    %CACHEDSTORAGE MATLAB spelling of Optuna's internal _CachedStorage.
    %   The leading underscore is not a legal MATLAB identifier. This class
    %   preserves the public storage contract and delegates persistence to
    %   the supplied backend while retaining stable study metadata locally.

    properties (SetAccess=immutable)
        Backend
    end

    properties (Access=private)
        StudyNames dictionary = dictionary(double.empty,string.empty)
        StudyDirections dictionary = dictionary(double.empty,cell.empty)
    end

    methods
        function obj=CachedStorage(backend)
            arguments
                backend (1,1) radia.optuna.BaseStorage
            end
            obj.Backend=backend;
        end

        function study_id=create_new_study(obj,directions,study_name)
            if nargin<3
                study_name="";
            end
            study_id=obj.Backend.create_new_study(directions,study_name);
            obj.StudyNames(study_id)=string(study_name);
            obj.StudyDirections(study_id)={directions};
        end

        function trial_id=create_new_trial(obj,study_id,template_trial)
            if nargin<3
                trial_id=obj.Backend.create_new_trial(study_id);
            else
                trial_id=obj.Backend.create_new_trial(study_id,template_trial);
            end
        end

        function delete_study(obj,study_id)
            if isKey(obj.StudyNames,study_id)
                obj.StudyNames=remove(obj.StudyNames,study_id);
            end
            if isKey(obj.StudyDirections,study_id)
                obj.StudyDirections=remove(obj.StudyDirections,study_id);
            end
            obj.Backend.delete_study(study_id);
        end

        function studies=get_all_studies(obj)
            studies=obj.Backend.get_all_studies();
        end

        function trials=get_all_trials(obj,study_id,deepcopy,states)
            if nargin<3
                trials=obj.Backend.get_all_trials(study_id);
            elseif nargin<4
                trials=obj.Backend.get_all_trials(study_id,deepcopy);
            else
                trials=obj.Backend.get_all_trials(study_id,deepcopy,states);
            end
        end

        function trial=get_best_trial(obj,study_id)
            trial=obj.Backend.get_best_trial(study_id);
        end

        function count=get_n_trials(obj,study_id,state)
            if nargin<3
                count=obj.Backend.get_n_trials(study_id);
            else
                count=obj.Backend.get_n_trials(study_id,state);
            end
        end

        function directions=get_study_directions(obj,study_id)
            if isKey(obj.StudyDirections,study_id)
                value=obj.StudyDirections(study_id);
                directions=value{1};
                return
            end
            directions=obj.Backend.get_study_directions(study_id);
            obj.StudyDirections(study_id)={directions};
        end

        function study_id=get_study_id_from_name(obj,study_name)
            study_id=obj.Backend.get_study_id_from_name(study_name);
        end

        function study_name=get_study_name_from_id(obj,study_id)
            if isKey(obj.StudyNames,study_id)
                study_name=obj.StudyNames(study_id);
                if strlength(study_name)>0
                    return
                end
            end
            study_name=obj.Backend.get_study_name_from_id(study_id);
            obj.StudyNames(study_id)=study_name;
        end

        function attributes=get_study_system_attrs(obj,study_id)
            attributes=obj.Backend.get_study_system_attrs(study_id);
        end

        function attributes=get_study_user_attrs(obj,study_id)
            attributes=obj.Backend.get_study_user_attrs(study_id);
        end

        function trial=get_trial(obj,trial_id)
            trial=obj.Backend.get_trial(trial_id);
        end

        function trial_id=get_trial_id_from_study_id_trial_number( ...
                obj,study_id,trial_number)
            trial_id=obj.Backend.get_trial_id_from_study_id_trial_number( ...
                study_id,trial_number);
        end

        function trial_number=get_trial_number_from_id(obj,trial_id)
            trial_number=obj.Backend.get_trial_number_from_id(trial_id);
        end

        function value=get_trial_param(obj,trial_id,param_name)
            value=obj.Backend.get_trial_param(trial_id,param_name);
        end

        function params=get_trial_params(obj,trial_id)
            params=obj.Backend.get_trial_params(trial_id);
        end

        function attributes=get_trial_system_attrs(obj,trial_id)
            attributes=obj.Backend.get_trial_system_attrs(trial_id);
        end

        function attributes=get_trial_user_attrs(obj,trial_id)
            attributes=obj.Backend.get_trial_user_attrs(trial_id);
        end

        function remove_session(obj)
            obj.Backend.remove_session();
        end

        function set_study_system_attr(obj,study_id,key,value)
            obj.Backend.set_study_system_attr(study_id,key,value);
        end

        function set_study_user_attr(obj,study_id,key,value)
            obj.Backend.set_study_user_attr(study_id,key,value);
        end

        function set_trial_intermediate_value(obj,trial_id,step,value)
            obj.Backend.set_trial_intermediate_value(trial_id,step,value);
        end

        function set_trial_param( ...
                obj,trial_id,param_name,param_value_internal,distribution)
            obj.Backend.set_trial_param( ...
                trial_id,param_name,param_value_internal,distribution);
        end

        function updated=set_trial_state_values(obj,trial_id,state,values)
            if nargin<4
                updated=obj.Backend.set_trial_state_values(trial_id,state);
            else
                updated=obj.Backend.set_trial_state_values( ...
                    trial_id,state,values);
            end
        end

        function set_trial_system_attr(obj,trial_id,key,value)
            obj.Backend.set_trial_system_attr(trial_id,key,value);
        end

        function set_trial_user_attr(obj,trial_id,key,value)
            obj.Backend.set_trial_user_attr(trial_id,key,value);
        end

        function check_trial_is_updatable(obj,trial_id,trial_state)
            obj.Backend.check_trial_is_updatable(trial_id,trial_state);
        end

        function record_heartbeat(obj,trial_id)
            if ismethod(obj.Backend,"record_heartbeat")
                obj.Backend.record_heartbeat(trial_id);
            end
        end

        function interval=get_heartbeat_interval(obj)
            interval=[];
            if ismethod(obj.Backend,"get_heartbeat_interval")
                interval=obj.Backend.get_heartbeat_interval();
            end
        end

        function callback=get_heartbeat_stale_trial_callback(obj)
            callback=[];
            if ismethod(obj.Backend,"get_heartbeat_stale_trial_callback")
                callback=obj.Backend.get_heartbeat_stale_trial_callback();
            end
        end

        function callback=get_failed_trial_callback(obj)
            warning("radia:optuna:FutureWarning", ...
                "get_failed_trial_callback has been deprecated in v4.9.0. " + ...
                "Use get_heartbeat_stale_trial_callback instead.");
            callback=obj.get_heartbeat_stale_trial_callback();
        end
    end
end
