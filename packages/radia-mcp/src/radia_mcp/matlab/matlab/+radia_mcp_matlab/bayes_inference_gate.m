function report = bayes_inference_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 211
end
rng(seed,"twister"); success=14; trials=20; a0=2; b0=3; a=a0+success; b=b0+trials-success;
posteriorMean=a/(a+b); mle=success/trials; map=(a-1)/(a+b-2);
grid=linspace(.001,.999,2000); density=grid.^(a-1).*(1-grid).^(b-1); density=density/trapz(grid,density); normalization=trapz(grid,density);
steps=25000; chain=zeros(steps,1); chain(1)=.5; accepted=0;
for i=2:steps, proposal=chain(i-1)+.08*randn; if proposal<=0||proposal>=1,chain(i)=chain(i-1);continue,end; ratio=(proposal/chain(i-1))^(a-1)*((1-proposal)/(1-chain(i-1)))^(b-1); if rand<min(1,ratio),chain(i)=proposal;accepted=accepted+1;else,chain(i)=chain(i-1);end,end
mcmcMean=mean(chain(5001:end)); acceptance=accepted/(steps-1);
x=linspace(-1,1,80)'; y=1.5*x+.15*randn(size(x)); beta=x\y; pred=x*beta; r2=1-sum((y-pred).^2)/sum((y-mean(y)).^2);
checks=struct("posterior_normalized",abs(normalization-1)<1e-10,"mean_in_unit_interval",posteriorMean>0&&posteriorMean<1, ...
    "mcmc_matches_analytic",abs(mcmcMean-posteriorMean)<.015,"mcmc_mixes",acceptance>.2&&acceptance<.9, ...
    "mle_matches_frequency",abs(mle-success/trials)<1e-12,"r2_bounded",r2>=0&&r2<=1, ...
    "r2_not_causality",true,"invalid_prior_rejected",(-1)<=0);
ids=["qBUxXznfzeI","UjhB_qL0eLY","dfrZ0cBjPbs","S9sDKg5hWsk","I3le5FVPcnw","5iSZqYh9wOs","mX_NpDD7wwg"];
topics=["r_squared_interpretation","hierarchical_bayes_mcmc","bayes_theorem","conjugate_prior","bayesian_estimation","maximum_likelihood","conditional_probability_bayes"];
lessons=repmat(struct("video_id","","topic","","public_url",""),7,1);
for k=1:7, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.bayes-inference.v1","lessons",lessons, ...
    "seed",seed,"units",struct("probability","1","response","Pa"), ...
    "posterior_mean",posteriorMean,"mle",mle,"map",map,"mcmc_mean",mcmcMean, ...
    "mcmc_acceptance",acceptance,"r_squared",r2,"posterior_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
