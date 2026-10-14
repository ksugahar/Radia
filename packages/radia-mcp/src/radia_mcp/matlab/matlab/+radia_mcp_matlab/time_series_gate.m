function report = time_series_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 181
end
rng(seed,"twister"); a=[.55 .25;0 .45]; n=4000; y=zeros(n,2); noise=.2*randn(n,2);
for t=2:n, y(t,:)=(a*y(t-1,:)')'+noise(t,:); end
x=y(1:end-1,:); target=y(2:end,:); ahat=(x\target)'; fitError=norm(ahat-a,"fro");
irf=zeros(8,2); shock=[1;0]; for h=0:7, irf(h+1,:)=(a^h*shock)'; end
full=target(:,1)-x*(x\target(:,1)); restricted=target(:,1)-x(:,1)*(x(:,1)\target(:,1)); grangerGain=1-sum(full.^2)/sum(restricted.^2);
innovation=target-x*ahat'; ac1=corrManual(innovation(1:end-1,1),innovation(2:end,1));
unstable=[1.05 0;0 .5]; checks=struct("var_recovered",fitError<.03,"stable_var",max(abs(eig(a)))<1, ...
    "unstable_negative_rejected",max(abs(eig(unstable)))>=1,"impulse_decays",norm(irf(end,:))<norm(irf(1,:)), ...
    "granger_gain_positive",grangerGain>.01,"innovation_nearly_white",abs(ac1)<.04);
ids=["e3bcZQzSuIc","fTql4nO8Dnw","Ja7OZzu1s7w","P_204LIcsFk","d0EGcXZlpJ4","HbAGHRbuvQ0"];
topics=["impulse_response","var_granger","arma_error","recurrence_eigendecomposition","recurrence_characteristic_equation","subscriber_forecast"];
lessons=repmat(struct("video_id","","topic","","public_url",""),6,1);
for k=1:6, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.time-series.v1","lessons",lessons, ...
    "seed",seed,"units",struct("signal","Pa","lag","sample"),"var_fit_error",fitError, ...
    "granger_sse_gain",grangerGain,"innovation_lag1_correlation",ac1,"impulse_response",irf, ...
    "causal_score_is_ground_truth",false,"promotion_requires_forward_solver",true, ...
    "checks",checks,"ok",all(structfun(@logical,checks)));
end
function r=corrManual(x,y), x=x-mean(x); y=y-mean(y); r=(x'*y)/(norm(x)*norm(y)); end
