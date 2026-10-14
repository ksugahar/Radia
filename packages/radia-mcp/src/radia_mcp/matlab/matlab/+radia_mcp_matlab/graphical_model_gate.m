function report = graphical_model_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 173
end
rng(seed,"twister"); px=[.55 .45]; pyGivenX=[.8 .2;.25 .75]; pzGivenY=[.9 .1;.2 .8]; joint=zeros(2,2,2);
for x=1:2, for y=1:2, for z=1:2, joint(x,y,z)=px(x)*pyGivenX(x,y)*pzGivenY(y,z); end, end, end
condError=0;
for y=1:2
    slice=squeeze(joint(:,y,:)); slice=slice/sum(slice,"all"); pxgy=sum(slice,2); pzgy=sum(slice,1); condError=max(condError,norm(slice-pxgy*pzgy,"fro"));
end
precision=[1 -.4 0;-.4 1 -.3;0 -.3 1]; covariance=inv(precision); samples=randn(30000,3)*chol(covariance);
r=corrcoef(samples); partial13=(r(1,3)-r(1,2)*r(2,3))/sqrt((1-r(1,2)^2)*(1-r(2,3)^2));
cycle=[0 1 0;0 0 1;1 0 0]; cycleRejected=trace(cycle^3)>0;
checks=struct("joint_normalized",abs(sum(joint,"all")-1)<1e-12,"conditional_independence",condError<1e-12, ...
    "precision_encodes_missing_edge",precision(1,3)==0,"sample_partial_correlation_small",abs(partial13)<.03, ...
    "cyclic_dag_rejected",cycleRejected);
lessons=[struct("video_id","knCbMFQJXxY","topic","bayesian_network","public_url","https://www.youtube.com/watch?v=knCbMFQJXxY"); ...
    struct("video_id","hh_KPDZ1D2Y","topic","graphical_model_spurious_correlation","public_url","https://www.youtube.com/watch?v=hh_KPDZ1D2Y")];
report=struct("schema","radia-mcp.matlab.graphical-model.v1","lessons",lessons, ...
    "seed",seed,"units",struct("probability","1","partial_correlation","1"), ...
    "conditional_independence_error",condError,"sample_partial_correlation",partial13, ...
    "causal_graph_is_solver_truth",false,"promotion_requires_forward_solver",true, ...
    "checks",checks,"ok",all(structfun(@logical,checks)));
end
