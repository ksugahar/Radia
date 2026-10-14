function seed = resolveSeed(value)
%RESOLVESEED Validate an explicit uint32 seed or draw fresh process entropy.
arguments
    value double = double.empty(1,0)
end
persistent recent
if isempty(value)
    if usejava("jvm")
        secureRandom=javaObject("java.security.SecureRandom");
        candidate=double(typecast(int32(secureRandom.nextInt()),"uint32"));
    else
        % This private stream does not mutate MATLAB's global RNG. The
        % clock-seeded fallback is for no-JVM MATLAB processes only.
        stream=RandStream("mt19937ar",Seed="shuffle");
        candidate=double(stream.Seed);
    end

    % SecureRandom collisions are already negligible. This bounded guard
    % also prevents clock-seed duplicates in a tight no-JVM constructor loop.
    while any(recent==candidate)
        candidate=mod(candidate+1,double(intmax("uint32"))+1);
    end
    recent(end+1)=candidate;
    if numel(recent)>1024
        recent=recent(end-1023:end);
    end
    seed=candidate;
    return
end
if ~(isscalar(value) && isfinite(value) && value==floor(value) && ...
        value>=0 && value<=double(intmax("uint32")))
    error("radia:optuna:Seed", ...
        "Seed must be empty (Optuna seed=None) or a uint32 integer.");
end
seed=double(value);
end
