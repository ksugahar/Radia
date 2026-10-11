function result=solveHCurlJouleLP(design,objectiveGradient,moveLimit,options)
%SOLVEHCURLJOULELP Trust-region LP in arbitrary analytic HCurl design modes.
arguments
    design (:,1) double {mustBeFinite}
    objectiveGradient (:,1) double {mustBeFinite}
    moveLimit double {mustBeNonnegative,mustBeFinite}
    options.LowerBounds double=-inf(size(design))
    options.UpperBounds double=inf(size(design))
    options.VolumeGradient double=double.empty
    options.VolumeLimit double=double.empty
end
n=numel(design);
if numel(objectiveGradient)~=n
    error("radia:topopt:Shape","objectiveGradient must match design.");
end
move=localExpand(moveLimit,n,"moveLimit");
lower=localExpand(options.LowerBounds,n,"LowerBounds");
upper=localExpand(options.UpperBounds,n,"UpperBounds");
lb=max(-move,lower-design); ub=min(move,upper-design);
Aineq=[]; bineq=[];
if ~isempty(options.VolumeGradient)||~isempty(options.VolumeLimit)
    if isempty(options.VolumeGradient)||~isscalar(options.VolumeLimit)
        error("radia:topopt:Volume", ...
            "VolumeGradient and scalar VolumeLimit must be supplied together.");
    end
    volumeGradient=options.VolumeGradient(:);
    if numel(volumeGradient)~=n
        error("radia:topopt:Shape","VolumeGradient must match design.");
    end
    Aineq=volumeGradient'; bineq=options.VolumeLimit;
end
settings=optimoptions("linprog","Display","none","Algorithm","dual-simplex-highs");
[delta,~,exitflag,output]=linprog(objectiveGradient,Aineq,bineq,[],[],lb,ub,settings);
if exitflag<=0
    error("radia:topopt:LPFailed","HCurl Joule LP failed: %s",output.message);
end
result=struct("schema","radia.hcurl.topopt.joule-lp/v1", ...
    "design",design+delta,"delta",delta, ...
    "objective_gradient",objectiveGradient,"status",string(output.message));
end

function value=localExpand(value,n,name)
value=value(:);
if isscalar(value), value=repmat(value,n,1); end
if numel(value)~=n
    error("radia:topopt:Shape","%s must be scalar or match design.",name);
end
end
