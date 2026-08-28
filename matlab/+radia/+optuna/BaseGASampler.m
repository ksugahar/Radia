classdef (Abstract) BaseGASampler < radia.optuna.BaseSampler
    %BASEGASAMPLER Optuna generation/population management contract.

    properties (Access=private)
        BasePopulationSize double = double.empty(1,0)
    end

    properties (Dependent)
        population_size
    end

    methods
        function obj=BaseGASampler(populationSize)
            if nargin<1
                populationSize=double.empty(1,0);
            end
            if ~isempty(populationSize) && (~isscalar(populationSize) || ...
                    populationSize~=floor(populationSize) || populationSize<=0)
                error("radia:optuna:GAPopulation", ...
                    "population_size must be empty or a positive integer.");
            end
            obj.BasePopulationSize=double(populationSize);
        end

        function value=get.population_size(obj)
            value=obj.BasePopulationSize;
        end

        function set.population_size(obj,value)
            if ~isscalar(value) || value~=floor(value) || value<=0
                error("radia:optuna:GAPopulation", ...
                    "population_size must be a positive integer.");
            end
            obj.BasePopulationSize=double(value);
        end

        function generation=get_trial_generation(obj,study,trial)
            key=obj.generationKey();
            if isfield(trial.SystemAttrs,key)
                generation=double(trial.SystemAttrs.(key));
                return
            end
            complete=study.get_trials("COMPLETE");
            generations=-ones(1,numel(complete));
            for index=1:numel(complete)
                if isfield(complete(index).SystemAttrs,key)
                    generations(index)=double(complete(index).SystemAttrs.(key));
                end
            end
            maxGeneration=max([0,generations]);
            count=sum(generations==maxGeneration);
            if isempty(obj.population_size)
                error("radia:optuna:GAPopulation", ...
                    "Population size must be set.");
            end
            generation=maxGeneration+(count>=obj.population_size);
            trial.set_system_attr(obj.generationAttribute(),generation);
        end

        function population=get_population(obj,study,generation)
            complete=study.get_trials("COMPLETE");
            key=obj.generationKey();
            keep=false(1,numel(complete));
            for index=1:numel(complete)
                keep(index)=isfield(complete(index).SystemAttrs,key) && ...
                    double(complete(index).SystemAttrs.(key))==generation;
            end
            population=complete(keep);
        end

        function population=get_parent_population(obj,study,generation)
            if generation==0
                population=radia.optuna.FrozenTrial.empty(0,1);
                return
            end
            key=matlab.lang.makeValidName(obj.parentAttribute(generation));
            attributes=study.system_attrs();
            if isfield(attributes,key)
                numbers=reshape(double(attributes.(key)),1,[]);
                trials=study.get_trials();
                population=trials(ismember([trials.Number],numbers));
                return
            end
            selected=obj.select_parent(study,generation);
            numbers=reshape([selected.Number],1,[]);
            study.set_system_attr(obj.parentAttribute(generation), ...
                numbers);
            trials=study.get_trials();
            population=trials(ismember([trials.Number],numbers));
        end
    end

    methods (Abstract)
        population=select_parent(obj,study,generation)
    end

    methods (Access=private)
        function value=generationAttribute(obj)
            parts=split(string(class(obj)),".");
            value=parts(end)+":generation";
        end

        function value=generationKey(obj)
            value=matlab.lang.makeValidName(obj.generationAttribute());
        end

        function value=parentAttribute(obj,generation)
            parts=split(string(class(obj)),".");
            value=parts(end)+":parent:"+string(generation);
        end
    end
end
