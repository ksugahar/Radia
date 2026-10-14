function report = quantification_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 127
end
rng(seed,"twister"); a=[1 1 0 0;1 0 1 0;0 1 0 1;0 0 1 1;1 1 1 0;0 1 1 1];
p=a/sum(a,"all"); pr=sum(p,2); pc=sum(p,1); residual=(p-pr*pc)./sqrt(pr*pc);
[u,s,v]=svd(residual,"econ"); reconstruction=norm(residual-u*s*v',"fro");
top=v(:,1); lambda=s(1,1)^2; gram=residual'*residual; lagrangeResidual=norm(gram*top-lambda*top);
rq=zeros(500,1); for i=1:500, x=randn(4,1); x=x/norm(x); rq(i)=x'*gram*x; end
rowScore=(diag(1./sqrt(pr))*u(:,1))*s(1,1); colScore=(diag(1./sqrt(pc))*v(:,1))*s(1,1);
checks=struct("probability_normalized",abs(sum(p,"all")-1)<1e-12, ...
    "svd_reconstruction",reconstruction<1e-12,"lagrange_eigen_equation",lagrangeResidual<1e-12, ...
    "top_rayleigh_maximal",lambda>=max(rq)-1e-12,"scores_finite",all(isfinite([rowScore;colScore])));
ids=["pAOV6RWt-xM","dWSdlX8-YMc","Adm9vv61mjg","49-dUOnzyTo","SfUtRvsjcqs"];
topics=["quantification_iii_roundtrip","quantification_iii_svd","lagrange_optimization","preference_correlation","quantification_iii_application"];
lessons=repmat(struct("video_id","","topic","","public_url",""),5,1);
for ii=1:5, lessons(ii)=struct("video_id",ids(ii),"topic",topics(ii),"public_url","https://www.youtube.com/watch?v="+ids(ii)); end
report=struct("schema","radia-mcp.matlab.quantification.v1","lessons",lessons, ...
    "seed",seed,"units",struct("preference","1","score","1"),"singular_values",diag(s)', ...
    "reconstruction_error",reconstruction,"lagrange_residual",lagrangeResidual, ...
    "row_scores",rowScore',"column_scores",colScore',"score_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
