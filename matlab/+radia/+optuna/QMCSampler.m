classdef QMCSampler < handle
    %QMCSAMPLER Sobol or Halton sampling for a fixed numeric search space.

    properties (SetAccess=private)
        QMCType (1,1) string = "sobol"
        Scramble (1,1) logical = false
        Seed (1,1) double = 0
        IndependentSampler
    end

    properties (Access=private)
        AttachedStudy = []
        Restored (1,1) logical = false
        NextSampleId (1,1) double = 0
    end

    properties (Constant, Access=private)
        StateSchema = "radia.optuna.qmc-sampler-state.v1"
        SamplerName = "qmc"
    end

    methods
        function obj=QMCSampler(options)
            arguments
                options.QMCType (1,1) string = "sobol"
                options.Scramble (1,1) logical = false
                options.Seed (1,1) double = 0
                options.IndependentSampler = []
            end
            qmcType=lower(options.QMCType);
            if ~ismember(qmcType,["sobol","halton"])
                error("radia:optuna:QMCType", ...
                    "QMCType must be 'sobol' or 'halton'.");
            end
            obj.QMCType=qmcType;
            obj.Scramble=options.Scramble;
            obj.Seed=options.Seed;
            if isempty(options.IndependentSampler)
                obj.IndependentSampler=radia.optuna.RandomSampler(options.Seed);
            else
                obj.IndependentSampler=options.IndependentSampler;
            end
        end

        function beforeTrial(obj,study,trial)
            if ismethod(obj.IndependentSampler,"beforeTrial")
                obj.IndependentSampler.beforeTrial(study,trial);
            end
            obj.attach(study);
            searchSpace=obj.inferRelativeSearchSpace(study,trial);
            if isempty(searchSpace)
                return
            end
            point=obj.generatePoint(numel(searchSpace),obj.NextSampleId);
            values=cell(1,numel(searchSpace));
            for index=1:numel(searchSpace)
                values{index}=obj.untransform( ...
                    point(index),searchSpace(index).distribution);
            end
            trial.setRelativeParameters(searchSpace,values,"qmc");
            obj.NextSampleId=obj.NextSampleId+1;
            obj.recordState(study,trial.Number);
        end

        function searchSpace=inferRelativeSearchSpace(~,study,~)
            searchSpace=struct("name",{},"distribution",{});
            if isempty(study.ParamTable)
                return
            end
            first=min(study.ParamTable.TrialNumber);
            rows=find(study.ParamTable.TrialNumber==first)';
            rows=rows(study.ParamTable.Kind(rows)~="categorical");
            if isempty(rows)
                return
            end
            firstDistribution= ...
                radia.optuna.internal.DistributionCodec.decode( ...
                study.ParamTable.Kind(rows(1)), ...
                study.ParamTable.Distribution(rows(1)));
            searchSpace=repmat(struct("name","", ...
                "distribution",firstDistribution),1,numel(rows));
            for index=1:numel(rows)
                row=rows(index);
                distribution=radia.optuna.internal.DistributionCodec.decode( ...
                    study.ParamTable.Kind(row),study.ParamTable.Distribution(row));
                searchSpace(index)=struct( ...
                    "name",study.ParamTable.Name(row), ...
                    "distribution",distribution);
            end
        end

        function value=sampleFloat(obj,study,trial,name,low,high,options)
            value=obj.IndependentSampler.sampleFloat( ...
                study,trial,name,low,high,options);
        end

        function value=sampleInteger(obj,study,trial,name,low,high)
            value=obj.IndependentSampler.sampleInteger( ...
                study,trial,name,low,high);
        end

        function value=sampleCategorical(obj,study,trial,name,choices)
            value=obj.IndependentSampler.sampleCategorical( ...
                study,trial,name,choices);
        end

        function values=sampleJoint(obj,study,trial,names,lows,highs,options)
            values=zeros(1,numel(names));
            for index=1:numel(names)
                values(index)=obj.IndependentSampler.sampleFloat( ...
                    study,trial,names(index),lows(index),highs(index), ...
                    struct("Log",options.Log(index),"Step",NaN));
            end
        end

        function afterTrial(obj,study,trial)
            if ismethod(obj.IndependentSampler,"afterTrial")
                obj.IndependentSampler.afterTrial(study,trial);
            end
            obj.recordState(study,trial.Number);
        end
    end

    methods (Access=private)
        function attach(obj,study)
            changed=isempty(obj.AttachedStudy) || ...
                ~isequal(obj.AttachedStudy,study);
            if changed
                obj.AttachedStudy=study;
                obj.NextSampleId=0;
                obj.Restored=false;
            end
            if obj.Restored
                return
            end
            state=study.samplerState(obj.SamplerName,obj.StateSchema);
            if ~isempty(state)
                valid=isstruct(state) && isscalar(state) && ...
                    isfield(state,"schema") && isfield(state,"qmc_type") && ...
                    isfield(state,"scramble") && isfield(state,"seed") && ...
                    isfield(state,"next_sample_id") && ...
                    string(state.schema)==obj.StateSchema && ...
                    string(state.qmc_type)==obj.QMCType && ...
                    logical(state.scramble)==obj.Scramble && ...
                    double(state.seed)==obj.Seed;
                if ~valid
                    error("radia:optuna:QMCState", ...
                        "Stored QMC sampler state is invalid or incompatible.");
                end
                obj.NextSampleId=double(state.next_sample_id);
            end
            obj.Restored=true;
        end

        function recordState(obj,study,trialNumber)
            obj.attach(study);
            state=struct("schema",obj.StateSchema,"qmc_type",obj.QMCType, ...
                "scramble",obj.Scramble,"seed",obj.Seed, ...
                "next_sample_id",obj.NextSampleId);
            study.recordSamplerState(obj.SamplerName,obj.StateSchema, ...
                trialNumber,0,state);
        end

        function point=generatePoint(obj,dimension,sampleId)
            if obj.QMCType=="sobol"
                point=obj.sobolPoint(dimension,sampleId);
            else
                point=obj.haltonPoint(dimension,sampleId);
            end
            if obj.Scramble
                stream=RandStream("mt19937ar","Seed",obj.Seed);
                point=mod(point+rand(stream,1,dimension),1);
            end
        end

        function point=sobolPoint(~,dimension,sampleId)
            if dimension>32
                error("radia:optuna:QMCDimension", ...
                    "The MATLAB-native Sobol implementation supports at most 32 dimensions.");
            end
            [polynomials,initialValues]= ...
                radia.optuna.QMCSampler.sobolParameters();
            gray=bitxor(uint32(sampleId),bitshift(uint32(sampleId),-1));
            point=zeros(1,dimension);
            for dim=1:dimension
                directions=zeros(1,32,"uint32");
                if dim==1
                    for bit=1:32
                        directions(bit)=bitshift(uint32(1),32-bit);
                    end
                else
                    polynomial=uint32(polynomials(dim));
                    degree=floor(log2(double(polynomial)));
                    initial=uint32(initialValues{dim});
                    for bit=1:degree
                        directions(bit)=bitshift(initial(bit),32-bit);
                    end
                    coefficients=bitand(bitshift(polynomial,-1), ...
                        uint32(2^(degree-1)-1));
                    for bit=(degree+1):32
                        value=bitxor(directions(bit-degree), ...
                            bitshift(directions(bit-degree),-degree));
                        for offset=1:(degree-1)
                            mask=bitshift(uint32(1),degree-1-offset);
                            if bitand(coefficients,mask)~=0
                                value=bitxor(value,directions(bit-offset));
                            end
                        end
                        directions(bit)=value;
                    end
                end
                integer=uint32(0);
                for bit=1:32
                    if bitget(gray,bit)
                        integer=bitxor(integer,directions(bit));
                    end
                end
                point(dim)=double(integer)/2^32;
            end
        end

        function point=haltonPoint(~,dimension,sampleId)
            bases=radia.optuna.QMCSampler.firstPrimes(dimension);
            point=zeros(1,dimension);
            for dim=1:dimension
                index=sampleId;
                factor=1/bases(dim);
                while index>0
                    point(dim)=point(dim)+factor*mod(index,bases(dim));
                    index=floor(index/bases(dim));
                    factor=factor/bases(dim);
                end
            end
        end

        function value=untransform(obj,fraction,distribution)
            low=distribution.low;
            high=distribution.high;
            step=distribution.step;
            if distribution.log
                if distribution.kind=="integer"
                    lower=log(low-0.5*step);
                    upper=log(high+0.5*step);
                else
                    lower=log(low);
                    upper=log(high);
                end
                raw=exp(lower+fraction*(upper-lower));
            else
                halfStep=0;
                if isfinite(step), halfStep=0.5*step; end
                raw=(low-halfStep)+fraction*((high+halfStep)-(low-halfStep));
            end
            if distribution.kind=="integer" && distribution.log
                value=obj.roundTiesToEven(raw);
                value=min(max(value,low),high);
            elseif isfinite(step)
                value=low+obj.roundTiesToEven((raw-low)/step)*step;
                value=min(max(value,low),high);
            else
                value=min(raw,obj.nextDown(high));
            end
        end

        function value=roundTiesToEven(~,value)
            lower=floor(value);
            fraction=value-lower;
            value=round(value);
            if fraction==0.5
                if mod(lower,2)==0
                    value=lower;
                else
                    value=lower+1;
                end
            end
        end

        function value=nextDown(~,value)
            bits=typecast(double(value),"uint64");
            if value>0
                bits=bits-uint64(1);
            elseif value<0
                bits=bits+uint64(1);
            else
                bits=bitor(bitshift(uint64(1),63),uint64(1));
            end
            value=typecast(bits,"double");
        end
    end

    methods (Static, Access=private)
        function [polynomials,initialValues]=sobolParameters()
            % Joe-Kuo direction-number seeds used by SciPy/Optuna.
            polynomials=[1 3 7 11 13 19 25 37 41 47 55 59 61 67 91 97 ...
                103 109 115 131 137 143 145 157 167 171 185 191 193 203 211 213];
            initialValues={[],1,[1 3],[1 3 1],[1 1 1],[1 1 3 3], ...
                [1 3 5 13],[1 1 5 5 17],[1 1 5 5 5],[1 1 7 11 19], ...
                [1 1 5 1 1],[1 1 1 3 11],[1 3 5 5 31], ...
                [1 3 3 9 7 49],[1 1 1 15 21 21],[1 3 1 13 27 49], ...
                [1 1 1 15 7 5],[1 3 1 15 13 25],[1 1 5 5 19 61], ...
                [1 3 7 11 23 15 103],[1 3 7 13 13 15 69], ...
                [1 1 3 13 7 35 63],[1 3 5 9 1 25 53], ...
                [1 3 1 13 9 35 107],[1 3 1 5 27 61 31], ...
                [1 1 5 11 19 41 61],[1 3 5 3 3 13 69], ...
                [1 1 7 13 1 19 1],[1 3 7 5 13 19 59], ...
                [1 1 3 9 25 29 41],[1 3 5 13 23 1 55], ...
                [1 3 7 3 13 59 17]};
        end

        function primes=firstPrimes(count)
            primes=zeros(1,count);
            candidate=2;
            found=0;
            while found<count
                prime=true;
                for divisor=2:floor(sqrt(candidate))
                    if mod(candidate,divisor)==0
                        prime=false;
                        break
                    end
                end
                if prime
                    found=found+1;
                    primes(found)=candidate;
                end
                candidate=candidate+1;
            end
        end
    end
end
