function report = gradient_boost_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 191
end
rng(seed,"twister"); x=linspace(-1,1,300)'; y=x.^2+.2*sin(5*x); pred=mean(y)*ones(size(y)); eta=.2; loss=zeros(20,1);
for m=1:20
    residual=y-pred; [stump,threshold]=fitStump(x,residual); pred=pred+eta*stump; loss(m)=mean((y-pred).^2);
end
[badStep,~]=fitStump(x,y-mean(y)); badPred=mean(y)+3*badStep; badLoss=mean((y-badPred).^2);
checks=struct("loss_decreases",loss(end)<loss(1),"mostly_monotone",sum(diff(loss)>1e-12)==0, ...
    "overstep_negative_worse",badLoss>loss(1),"threshold_in_domain",threshold>=min(x)&&threshold<=max(x), ...
    "prediction_finite",all(isfinite(pred)));
report=struct("schema","radia-mcp.matlab.gradient-boost.v1", ...
    "lesson",struct("video_id","u0IIqeNZOXY","topic","gradient_boosted_trees","public_url","https://www.youtube.com/watch?v=u0IIqeNZOXY"), ...
    "seed",seed,"units",struct("feature","1","response","Pa"),"learning_rate",eta,"loss_history",loss', ...
    "overstep_loss",badLoss,"boosted_prediction_is_ground_truth",false,"promotion_requires_forward_solver",true, ...
    "checks",checks,"ok",all(structfun(@logical,checks)));
end
function [best,thr]=fitStump(x,r)
candidates=x(2:end-1); bestLoss=inf; best=zeros(size(r)); thr=candidates(1);
for k=1:10:numel(candidates), t=candidates(k); left=x<=t; p=zeros(size(r)); p(left)=mean(r(left)); p(~left)=mean(r(~left)); l=mean((r-p).^2); if l<bestLoss,bestLoss=l;best=p;thr=t;end,end
end
