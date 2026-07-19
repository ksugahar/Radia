function report = quantification_ii_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 137
end
rng(seed,"twister"); x=[randn(80,3)+[-1 0 .3];randn(80,3)+[1 .2 -.3]]; label=[ones(80,1);2*ones(80,1)];
mu=mean(x); sw=zeros(3); sb=zeros(3);
for g=1:2, xg=x(label==g,:); mg=mean(xg); sw=sw+(xg-mg)'*(xg-mg); sb=sb+size(xg,1)*(mg-mu)'*(mg-mu); end
sw=sw+1e-9*eye(3); [w,d]=eig(sb,sw,"vector"); [lambda,idx]=max(real(d)); direction=real(w(:,idx)); direction=direction/sqrt(direction'*sw*direction);
lagrangeResidual=norm(sb*direction-lambda*sw*direction); ratio=lambda; randomRatio=zeros(500,1);
for k=1:500, a=randn(3,1); randomRatio(k)=(a'*sb*a)/(a'*sw*a); end
checks=struct("covariances_symmetric",norm(sw-sw',"fro")<1e-12&&norm(sb-sb',"fro")<1e-12, ...
    "generalized_eigen_residual",lagrangeResidual<1e-8,"ratio_maximal",ratio>=max(randomRatio)-1e-10, ...
    "normalized_constraint",abs(direction'*sw*direction-1)<1e-10);
ids=["GCPEMJe6epY","aXx9LScCHzY","XjirFbccIkQ","7gy5jplHExI","2-E4XiHQEcM","gqrrsAxM3r0"];
topics=["ratio_maximization","orthogonal_diagonalization","lagrange_eigenproblem","data_covariance","lagrange_meaning","quantification_ii_application"];
lessons=repmat(struct("video_id","","topic","","public_url",""),6,1);
for k=1:6, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.quantification-ii.v1","lessons",lessons, ...
    "seed",seed,"units",struct("feature","1","ratio","1"),"max_ratio",ratio, ...
    "max_random_ratio",max(randomRatio),"lagrange_residual",lagrangeResidual, ...
    "direction_is_ground_truth",false,"promotion_requires_forward_solver",true, ...
    "checks",checks,"ok",all(structfun(@logical,checks)));
end
