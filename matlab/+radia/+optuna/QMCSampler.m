classdef QMCSampler < radia.optuna.BaseSampler
    %QMCSAMPLER Sobol or Halton sampling for a fixed numeric search space.

    properties (SetAccess=private)
        QMCType (1,1) string = "sobol"
        Scramble (1,1) logical = false
        Seed (1,1) double = 0
        IndependentSampler
        WarnAsynchronousSeeding (1,1) logical = true
        WarnIndependentSampling (1,1) logical = true
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
                options.Seed double = double.empty(1,0)
                options.IndependentSampler = []
                options.WarnAsynchronousSeeding (1,1) logical = true
                options.WarnIndependentSampling (1,1) logical = true
            end
            qmcType=lower(options.QMCType);
            if ~ismember(qmcType,["sobol","halton"])
                error("radia:optuna:QMCType", ...
                    "QMCType must be 'sobol' or 'halton'.");
            end
            obj.QMCType=qmcType;
            obj.Scramble=options.Scramble;
            unseeded=isempty(options.Seed);
            obj.Seed=radia.optuna.internal.resolveSeed(options.Seed);
            obj.WarnAsynchronousSeeding=options.WarnAsynchronousSeeding;
            obj.WarnIndependentSampling=options.WarnIndependentSampling;
            if unseeded && obj.Scramble && obj.WarnAsynchronousSeeding
                warning("radia:optuna:QMCAsynchronousSeeding", ...
                    "Scrambled QMC with seed=None may use a different sequence in each parallel worker.");
            end
            if isempty(options.IndependentSampler)
                obj.IndependentSampler=radia.optuna.RandomSampler(options.Seed);
            else
                if ~isa(options.IndependentSampler,"radia.optuna.BaseSampler")
                    error("radia:optuna:QMCIndependentSampler", ...
                        "IndependentSampler must derive from radia.optuna.BaseSampler.");
                end
                obj.IndependentSampler=options.IndependentSampler;
            end
        end

        function reseed_rng(obj)
            obj.IndependentSampler.reseed_rng();
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
            obj.warnIndependent(study,trial,name);
            value=obj.IndependentSampler.sampleFloat( ...
                study,trial,name,low,high,options);
        end

        function value=sampleInteger(obj,study,trial,name,low,high)
            obj.warnIndependent(study,trial,name);
            value=obj.IndependentSampler.sampleInteger( ...
                study,trial,name,low,high);
        end

        function value=sampleCategorical(obj,study,trial,name,choices)
            obj.warnIndependent(study,trial,name);
            value=obj.IndependentSampler.sampleCategorical( ...
                study,trial,name,choices);
        end

        function values=sampleJoint(obj,study,trial,names,lows,highs,options)
            values=zeros(1,numel(names));
            for index=1:numel(names)
                obj.warnIndependent(study,trial,names(index));
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

        function points=unitPoints(obj,dimension,count,options)
            %UNITPOINTS Generate a contiguous QMC block without a Study.
            arguments
                obj
                dimension (1,1) double {mustBeInteger,mustBePositive}
                count (1,1) double {mustBeInteger,mustBeNonnegative}
                options.StartIndex (1,1) double ...
                    {mustBeInteger,mustBeNonnegative} = 0
            end
            sampleIds=options.StartIndex+(0:count-1);
            points=obj.generatePoints(dimension,sampleIds);
        end
    end

    methods (Access=private)
        function warnIndependent(obj,study,trial,name)
            if ~obj.WarnIndependentSampling
                return
            end
            prior=study.TrialTable.TrialNumber<trial.Number & ...
                ismember(study.TrialTable.State, ...
                ["RUNNING","WAITING","COMPLETE","PRUNED"]);
            if any(prior)
                warning("radia:optuna:QMCIndependentSampling", ...
                    "Parameter '%s' in trial %d is sampled independently by %s.", ...
                    name,trial.Number,class(obj.IndependentSampler));
            end
        end

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
            point=obj.generatePoints(dimension,sampleId);
        end

        function points=generatePoints(obj,dimension,sampleIds)
            sampleIds=reshape(double(sampleIds),[],1);
            points=zeros(numel(sampleIds),dimension);
            if obj.Scramble
                % Optuna 4.9 delegates scrambled Sobol/Halton generation to
                % scipy.stats.qmc.  Its seeded PCG64 scrambling is part of
                % the observable proposal sequence, so use that same oracle
                % instead of maintaining a second scramble implementation.
                points=obj.scipyScrambledPoints(dimension,sampleIds);
                return
            end
            if obj.QMCType~="sobol"
                for row=1:numel(sampleIds)
                    points(row,:)=obj.haltonPoint( ...
                        dimension,sampleIds(row),obj.Scramble,obj.Seed);
                end
                return
            end
            if dimension>32
                error("radia:optuna:QMCDimension", ...
                    "The MATLAB-native Sobol implementation supports at most 32 dimensions.");
            end
            [polynomials,initialValues]= ...
                radia.optuna.QMCSampler.sobolParameters();
            directions=obj.sobolDirections( ...
                dimension,polynomials,initialValues);
            shift=zeros(1,dimension,"uint32");
            if obj.Scramble
                [directions,shift]=obj.scrambleSobolDirections( ...
                    directions,obj.Seed);
            end
            for row=1:numel(sampleIds)
                gray=bitxor(uint32(sampleIds(row)), ...
                    bitshift(uint32(sampleIds(row)),-1));
                for dim=1:dimension
                    integer=shift(dim);
                    for bit=1:32
                        if bitget(gray,bit)
                            integer=bitxor(integer,directions(dim,bit));
                        end
                    end
                    points(row,dim)=double(integer)/2^32;
                end
            end
        end

        function points=scipyScrambledPoints(obj,dimension,sampleIds)
            if isempty(sampleIds)
                points=zeros(0,dimension);
                return
            end
            if any(sampleIds<0 | sampleIds~=floor(sampleIds))
                error("radia:optuna:QMCIndex", ...
                    "QMC sample IDs must be nonnegative integers.");
            end
            try
                arguments=pyargs("d",int32(dimension), ...
                    "scramble",true,"seed",int64(obj.Seed));
                maximum=max(sampleIds);
                if obj.QMCType=="sobol"
                    engine=py.scipy.stats.qmc.Sobol(arguments);
                    exponent=ceil(log2(maximum+1));
                    generated=double(engine.random_base2(int32(exponent)));
                else
                    engine=py.scipy.stats.qmc.Halton(arguments);
                    generated=double(engine.random(int64(maximum+1)));
                end
                points=generated(sampleIds+1,:);
            catch exception
                cause=MException("radia:optuna:QMCSciPy", ...
                    "Scrambled QMCSampler requires the configured Python " + ...
                    "environment with scipy.stats.qmc (Optuna 4.9 oracle).");
                cause=addCause(cause,exception);
                throw(cause);
            end
        end

        function point=sobolPoint(obj,dimension,sampleId,scrambled,seed)
            if dimension>32
                error("radia:optuna:QMCDimension", ...
                    "The MATLAB-native Sobol implementation supports at most 32 dimensions.");
            end
            [polynomials,initialValues]= ...
                radia.optuna.QMCSampler.sobolParameters();
            directions=obj.sobolDirections( ...
                dimension,polynomials,initialValues);
            shift=zeros(1,dimension,"uint32");
            if scrambled
                [directions,shift]=obj.scrambleSobolDirections( ...
                    directions,seed);
            end
            gray=bitxor(uint32(sampleId),bitshift(uint32(sampleId),-1));
            point=zeros(1,dimension);
            for dim=1:dimension
                integer=shift(dim);
                for bit=1:32
                    if bitget(gray,bit)
                        integer=bitxor(integer,directions(dim,bit));
                    end
                end
                point(dim)=double(integer)/2^32;
            end
        end

        function directions=sobolDirections(~,dimension,polynomials,initialValues)
            directions=zeros(dimension,32,"uint32");
            for bit=1:32
                directions(1,bit)=bitshift(uint32(1),32-bit);
            end
            for dim=2:dimension
                polynomial=uint32(polynomials(dim));
                degree=floor(log2(double(polynomial)));
                initial=uint32(initialValues{dim});
                for bit=1:degree
                    directions(dim,bit)=bitshift(initial(bit),32-bit);
                end
                coefficients=bitand(bitshift(polynomial,-1), ...
                    uint32(2^(degree-1)-1));
                for bit=(degree+1):32
                    value=bitxor(directions(dim,bit-degree), ...
                        bitshift(directions(dim,bit-degree),-degree));
                    for offset=1:(degree-1)
                        mask=bitshift(uint32(1),degree-1-offset);
                        if bitand(coefficients,mask)~=0
                            value=bitxor(value,directions(dim,bit-offset));
                        end
                    end
                    directions(dim,bit)=value;
                end
            end
        end

        function [directions,shift]=scrambleSobolDirections(~,directions,seed)
            % Matousek linear-matrix scramble plus a digital XOR shift.
            dimension=size(directions,1);
            stream=RandStream("mt19937ar","Seed",seed);
            shift=zeros(1,dimension,"uint32");
            for dim=1:dimension
                shiftBits=rand(stream,1,32)>=0.5;
                shift(dim)=radia.optuna.QMCSampler.packBits(shiftBits);
                lower=false(32,32);
                for row=1:32
                    lower(row,1:row)=rand(stream,1,row)>=0.5;
                    lower(row,row)=true;
                end
                for column=1:32
                    % LMS is triangular in significance order (MSB first).
                    input=logical(bitget( ...
                        directions(dim,column),32:-1:1));
                    output=false(1,32);
                    for row=1:32
                        output(row)=mod(sum( ...
                            lower(row,1:row) & input(1:row)),2)==1;
                    end
                    directions(dim,column)= ...
                        radia.optuna.QMCSampler.packBits(fliplr(output));
                end
            end
        end

        function point=haltonPoint(~,dimension,sampleId,scrambled,seed)
            bases=radia.optuna.QMCSampler.firstPrimes(dimension);
            point=zeros(1,dimension);
            stream=RandStream("mt19937ar","Seed",seed);
            for dim=1:dimension
                permutations=cell(1,0);
                if scrambled
                    count=ceil(54/log2(bases(dim)))-1;
                    permutations=cell(1,count);
                    for digit=1:count
                        permutations{digit}=randperm(stream,bases(dim))-1;
                    end
                end
                index=sampleId;
                factor=1/bases(dim);
                digit=1;
                while index>0 || (scrambled && digit<=numel(permutations))
                    remainder=mod(index,bases(dim));
                    if scrambled
                        remainder=permutations{digit}(remainder+1);
                    end
                    point(dim)=point(dim)+factor*remainder;
                    index=floor(index/bases(dim));
                    factor=factor/bases(dim);
                    digit=digit+1;
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
            value=radia.optuna.internal.UpstreamNumerics.roundTiesToEven( ...
                value);
        end

        function value=nextDown(~,value)
            value=radia.optuna.internal.UpstreamNumerics.nextDown(value);
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

        function value=packBits(bits)
            value=uint32(0);
            for bit=1:min(32,numel(bits))
                if bits(bit)
                    value=bitor(value,bitshift(uint32(1),bit-1));
                end
            end
        end
    end
end
