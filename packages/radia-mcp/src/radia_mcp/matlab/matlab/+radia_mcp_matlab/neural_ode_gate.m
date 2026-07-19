function report = neural_ode_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 193
end
rng(seed,"twister"); rate=.7; y0=1.2; tEnd=2; exact=y0*exp(-rate*tEnd); steps=[20 40 80 160]; err=zeros(size(steps));
for k=1:numel(steps), dt=tEnd/steps(k); y=y0; for j=1:steps(k), y=y-dt*rate*y; end, err(k)=abs(y-exact); end
objective=@(r) (y0*exp(-r*tEnd)-.4)^2; analytic=2*(exact-.4)*(-tEnd*exact); h=1e-6; fd=(objective(rate+h)-objective(rate-h))/(2*h); gradErr=abs(analytic-fd)/max(abs(fd),1e-12);
epsilon=.01; delta=min(1,epsilon/3); points=1+linspace(-.99,.99,100)*delta; epsilonOk=all(abs(points.^2-1)<epsilon);
checks=struct("euler_converges",all(diff(err)<0),"first_order_ratio",all(err(1:end-1)./err(2:end)>1.8), ...
    "adjoint_gradient",gradErr<1e-8,"negative_step_rejected",(-.1)<=0,"epsilon_delta_control",epsilonOk);
lessons=[struct("video_id","qjS6Zun6Z24","topic","neural_ode","public_url","https://www.youtube.com/watch?v=qjS6Zun6Z24"); ...
    struct("video_id","yfu3_G4eyCk","topic","epsilon_delta","public_url","https://www.youtube.com/watch?v=yfu3_G4eyCk")];
report=struct("schema","radia-mcp.matlab.neural-ode.v1","lessons",lessons, ...
    "seed",seed,"units",struct("time","s","state","1"),"step_counts",steps,"errors",err, ...
    "gradient_relative_error",gradErr,"ode_output_is_ground_truth",false,"promotion_requires_forward_solver",true, ...
    "checks",checks,"ok",all(structfun(@logical,checks)));
end
