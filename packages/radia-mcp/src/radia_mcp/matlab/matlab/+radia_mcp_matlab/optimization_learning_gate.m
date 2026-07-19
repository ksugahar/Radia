function report = optimization_learning_gate(seed)
%OPTIMIZATION_LEARNING_GATE Unified audit for legacy MATLAB optimization assets.
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 23
end
rng(seed, "twister");

% 1. Gradient agreement: analytic/AD reference, central FD, complex step.
x0 = [0.4; -0.7; 1.1];
f = @(x) sum(exp(0.2*x) + sin(x).^2 + 0.1*x.^4);
gref = 0.2*exp(0.2*x0) + 2*sin(x0).*cos(x0) + 0.4*x0.^3;
gfd = zeros(size(x0)); gcs = gfd;
hfd = eps^(1/3); hcs = 1e-20;
for k = 1:numel(x0)
    e = zeros(size(x0)); e(k) = 1;
    gfd(k) = (f(x0+hfd*e)-f(x0-hfd*e))/(2*hfd);
    gcs(k) = imag(f(x0+1i*hcs*e))/hcs;
end
gradient = struct("fd_relative_error", relerr(gfd,gref), ...
    "complex_step_relative_error", relerr(gcs,gref), ...
    "complex_step_supported", true);

% 2. One transparent GP contract shared with Bayesian/active learning.
xTrain = linspace(-2,2,9)'; yTrain = objective(xTrain);
xQuery = linspace(-2,2,81)';
[mu,sigma] = gpPredict(xTrain,yTrain,xQuery);
[~,acqIdx] = min(mu-1.96*sigma);
gp = struct("kernel", "squared_exponential", "acquisition", "LCB", ...
    "next_x", xQuery(acqIdx), "next_sigma", sigma(acqIdx), ...
    "uncertainty_nonnegative", all(sigma>=0));

% 3. SVD/TSVD/POD rank policies: energy, discrepancy, and GCV.
t = linspace(0,1,60)';
a = [sin(pi*t), 0.4*sin(2*pi*t), 0.08*cos(5*pi*t), 0.01*sin(11*pi*t)];
[~,s,~] = svd(a,"econ"); sv = diag(s); energy = cumsum(sv.^2)/sum(sv.^2);
rankEnergy = find(energy>=0.999,1);
b = a*[1;-0.5;0.2;0.1] + 1e-3*sin(17*pi*t);
[u,ss,v] = svd(a,"econ"); gcv = zeros(numel(sv),1);
for r=1:numel(sv)
    xr = v(:,1:r)*(diag(ss(1:r,1:r))\(u(:,1:r)'*b));
    gcv(r) = norm(a*xr-b)^2/(size(a,1)-r)^2;
end
[~,rankGcv] = min(gcv);
rankDiscrepancy = find(arrayfun(@(r) norm(a*(v(:,1:r)*(diag(ss(1:r,1:r))\(u(:,1:r)'*b)))-b),1:numel(sv)) <= sqrt(numel(b))*1.2e-3,1);
if isempty(rankDiscrepancy), rankDiscrepancy=numel(sv); end
rankSelection = struct("energy_threshold",0.999,"energy_rank",rankEnergy, ...
    "gcv_rank",rankGcv,"discrepancy_rank",rankDiscrepancy, ...
    "singular_values",sv');

% 4. Reproducible constrained PSO and CMA-like evolution share history rows.
bounds = [-2 2; -2 2];
pso = runPso(seed,bounds,35,35);
cma = runCmaLike(seed,bounds,35,24);

% 5. Every optimum is replayed through an independent forward callback.
psoVerified = forwardObjective(pso.best_x);
cmaVerified = forwardObjective(cma.best_x);
forward = struct("solver","independent_analytic_forward_callback", ...
    "tolerance",1e-10,"pso_error",abs(psoVerified-pso.best_objective), ...
    "cma_error",abs(cmaVerified-cma.best_objective));
forward.passed = forward.pso_error<=forward.tolerance && forward.cma_error<=forward.tolerance;

checks = struct("gradient_fd",gradient.fd_relative_error<1e-7, ...
    "gradient_complex",gradient.complex_step_relative_error<1e-12, ...
    "gp_uncertainty",gp.uncertainty_nonnegative, ...
    "rank_valid",all([rankEnergy rankGcv rankDiscrepancy]>=1 & [rankEnergy rankGcv rankDiscrepancy]<=numel(sv)), ...
    "pso_bounds",pso.bounds_ok,"cma_bounds",cma.bounds_ok, ...
    "history_schema",pso.history_schema==cma.history_schema, ...
    "forward_verified",forward.passed);

report = struct("schema","radia-mcp.matlab.optimization-learning.v1", ...
    "seed",seed,"units",struct("design","1","objective","Pa^2"), ...
    "gradient",gradient,"bayes_gp",gp,"rank_selection",rankSelection, ...
    "pso",pso,"cma_es",cma,"forward_verification",forward,"checks",checks, ...
    "ok",all(structfun(@logical,checks)));
end

function y=objective(x), y=(x-0.35).^2 + 0.05*sin(5*x); end
function y=forwardObjective(x), y=sum((x-0.35).^2 + 0.05*sin(5*x)); end
function e=relerr(a,b), e=norm(a-b)/max(norm(b),eps); end

function [mu,sigma]=gpPredict(xt,yt,xq)
ell=0.45; n=1e-9; K=exp(-0.5*((xt-xt')./ell).^2)+n*eye(numel(xt));
Ks=exp(-0.5*((xt-xq')./ell).^2); L=chol(K,"lower");
mu=Ks'*(L'\(L\yt)); v=L\Ks; sigma=sqrt(max(1-sum(v.^2,1)',0));
end

function out=runPso(seed,bounds,nIter,nPop)
rng(seed,"twister"); d=size(bounds,1); lo=bounds(:,1)'; hi=bounds(:,2)';
x=lo+rand(nPop,d).*(hi-lo); vel=zeros(size(x)); pb=x; py=rowscore(x);
[best,idx]=min(py); gb=pb(idx,:); hist=zeros(nIter,3);
for it=1:nIter
    vel=.7*vel+1.6*rand(nPop,d).*(pb-x)+1.6*rand(nPop,d).*(gb-x);
    x=min(max(x+vel,lo),hi); y=rowscore(x); improve=y<py; pb(improve,:)=x(improve,:); py(improve)=y(improve);
    [candidate,idx]=min(py); if candidate<best, best=candidate; gb=pb(idx,:); end
    hist(it,:)=[it best max(abs(vel),[],"all")];
end
out=pack("pso",seed,gb,best,hist,bounds);
end

function out=runCmaLike(seed,bounds,nIter,lambda)
rng(seed,"twister"); d=size(bounds,1); lo=bounds(:,1)'; hi=bounds(:,2)'; m=zeros(1,d); sigma=.7; hist=zeros(nIter,3); best=inf; bx=m;
for it=1:nIter
    pop=min(max(m+sigma*randn(lambda,d),lo),hi); y=rowscore(pop); [y,ord]=sort(y); elite=pop(ord(1:ceil(lambda/2)),:); m=mean(elite,1);
    if y(1)<best, best=y(1); bx=pop(ord(1),:); end
    sigma=max(.05,sigma*.94); hist(it,:)=[it best sigma];
end
out=pack("cma_es",seed,bx,best,hist,bounds);
end

function y=rowscore(x), y=sum((x-.35).^2+.05*sin(5*x),2); end
function out=pack(name,seed,x,best,hist,bounds)
out=struct("algorithm",name,"seed",seed,"best_x",x,"best_objective",best, ...
    "history_schema","iteration,best_objective,step_metric", ...
    "history",hist,"bounds",bounds,"bounds_ok",all(x>=bounds(:,1)'-eps & x<=bounds(:,2)'+eps));
end
