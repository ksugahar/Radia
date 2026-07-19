function report = gan_design_gate(seed)
%GAN_DESIGN_GATE Detect mode collapse before CAE solver promotion.
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 43
end
rng(seed,"twister"); n=600;
centers=[-1 0;1 0;0 1.25];
real=sampleMixture(centers,repmat(n/3,1,3),.12);
good=sampleMixture(centers,repmat(n/3,1,3),.12);
collapsed=sampleMixture(centers(1,:),n,.12);

goodCoverage=modeCoverage(good,centers,.45);
collapsedCoverage=modeCoverage(collapsed,centers,.45);
goodAccuracy=discriminatorAccuracy(real,good,seed+1);
collapsedAccuracy=discriminatorAccuracy(real,collapsed,seed+2);

qoi=forwardSolver(good);
% A solver-derived acceptance band deliberately separates feasible and
% infeasible generated modes; discriminator quality alone cannot promote all.
safe=qoi<1.0;
checks=struct( ...
    "good_covers_all_modes",goodCoverage==size(centers,1), ...
    "collapsed_control_rejected",collapsedCoverage==1, ...
    "good_fools_discriminator",goodAccuracy<.65, ...
    "collapse_detected_by_discriminator",collapsedAccuracy>.65, ...
    "forward_solver_finite",all(isfinite(qoi)), ...
    "solver_accepts_some_not_all",any(safe)&&any(~safe));

lessons=[struct("video_id","PnwkEwccPEE","topic","gan","public_url","https://www.youtube.com/watch?v=PnwkEwccPEE"); ...
    struct("video_id","Ul5gVsx6dRI","topic","gan_game_intro","public_url","https://www.youtube.com/watch?v=Ul5gVsx6dRI")];
report=struct("schema","radia-mcp.matlab.gan-design-gate.v1", ...
    "lessons",lessons,"principle","generator and discriminator co-evolution needs collapse diagnostics", ...
    "seed",seed,"units",struct("design","1","qoi","Pa^2"), ...
    "expected_modes",size(centers,1),"good_mode_coverage",goodCoverage, ...
    "collapsed_mode_coverage",collapsedCoverage, ...
    "good_discriminator_accuracy",goodAccuracy, ...
    "collapsed_discriminator_accuracy",collapsedAccuracy, ...
    "solver_accepted_count",sum(safe),"candidate_count",n, ...
    "generated_candidate_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks, ...
    "ok",all(structfun(@logical,checks)));
end

function x=sampleMixture(centers,counts,sigma)
x=[];
for k=1:size(centers,1)
    x=[x;centers(k,:)+sigma*randn(counts(k),2)]; %#ok<AGROW>
end
x=x(randperm(size(x,1)),:);
end
function n=modeCoverage(x,c,r)
d=zeros(size(x,1),size(c,1));
for k=1:size(c,1), d(:,k)=vecnorm(x-c(k,:),2,2); end
n=sum(sum(d<r,1)>=20);
end
function acc=discriminatorAccuracy(real,fake,seed)
rng(seed,"twister"); n=min(size(real,1),size(fake,1));
x=[real(1:n,:);fake(1:n,:)]; y=[ones(n,1);zeros(n,1)];
ord=randperm(2*n); tr=ord(1:round(1.5*n)); te=ord(round(1.5*n)+1:end);
xt=[ones(numel(tr),1) x(tr,:)]; w=zeros(3,1);
for it=1:300
    p=1./(1+exp(-xt*w)); w=w-.25*(xt'*(p-y(tr)))/numel(tr);
end
pt=1./(1+exp(-[ones(numel(te),1) x(te,:)]*w));
acc=mean((pt>=.5)==y(te));
end
function q=forwardSolver(x), q=sum((x-[.15 .2]).^2,2)+.05*sin(3*x(:,1)); end
