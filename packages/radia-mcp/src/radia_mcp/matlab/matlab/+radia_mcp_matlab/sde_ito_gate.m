function report = sde_ito_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 179
end
rng(seed,"twister"); paths=8000; steps=500; tEnd=1; dt=tEnd/steps; w=zeros(paths,1); integral=zeros(paths,1);
for k=1:steps, dw=sqrt(dt)*randn(paths,1); integral=integral+w.*dw; w=w+dw; end
itoError=sqrt(mean((w.^2-(2*integral+tEnd)).^2)); isometryEmp=mean(integral.^2); isometryExact=tEnd^2/2;
isometrySe=std(integral.^2)/sqrt(paths);
s=100; strike=95; rate=.03; sigma=.2; tau=.75; d1=(log(s/strike)+(rate+.5*sigma^2)*tau)/(sigma*sqrt(tau)); d2=d1-sigma*sqrt(tau);
n1=normalCdf(d1); n2=normalCdf(d2); phi=exp(-.5*d1^2)/sqrt(2*pi); price=s*n1-strike*exp(-rate*tau)*n2;
delta=n1; gamma=phi/(s*sigma*sqrt(tau)); theta=-(s*phi*sigma/(2*sqrt(tau))+rate*strike*exp(-rate*tau)*n2);
pdeResidual=theta+.5*sigma^2*s^2*gamma+rate*s*delta-rate*price;
checks=struct("ito_identity_converged",itoError<.09,"ito_isometry",abs(isometryEmp-isometryExact)<3*isometrySe, ...
    "option_price_positive",price>0,"black_scholes_pde",abs(pdeResidual)<1e-10, ...
    "invalid_volatility_rejected",(-.1)<=0);
ids=["_XwVB11x034","PrYPWD8nEmI","H1UmUtppWg8","CFrochh5UOg","HHOyFzn6NZ0","mrExmReKrcM","X2f1VsGCdDE","3rtFlKhGj3c","NE1W0wJH8q8"];
topics=["ito_integral_final","ito_integral","ito_formula","black_scholes_solution","black_scholes_derivation","ito_formula_intro","ito_integral_intro","sde_vs_ode","brownian_motion"];
lessons=repmat(struct("video_id","","topic","","public_url",""),9,1);
for k=1:9, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.sde-ito.v1","lessons",lessons, ...
    "seed",seed,"units",struct("time","s","state","1","price","currency"), ...
    "ito_rms_error",itoError,"ito_isometry_empirical",isometryEmp,"ito_isometry_exact",isometryExact, ...
    "ito_isometry_standard_error",isometrySe, ...
    "black_scholes_price",price,"pde_residual",pdeResidual,"stochastic_result_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
function p=normalCdf(x), p=.5*erfc(-x/sqrt(2)); end
