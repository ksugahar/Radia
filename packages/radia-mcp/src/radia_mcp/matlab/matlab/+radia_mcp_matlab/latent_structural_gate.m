function report = latent_structural_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 167
end
rng(seed,"twister"); groups=12; per=5; trueGroup=.5*randn(groups,1); obs=trueGroup+randn(groups,per); raw=mean(obs,2);
noiseVar=1/per; priorVar=.25; posterior=(raw/noiseVar)/(1/noiseVar+1/priorVar);
rawMse=mean((raw-trueGroup).^2); posteriorMse=mean((posterior-trueGroup).^2);
b=[0 .25;0 0]; gamma=[.8;.4]; sigmaE=diag([.2 .1]); transform=inv(eye(2)-b); implied=transform*(gamma*gamma'+sigmaE)*transform';
samples=randn(50000,1)*gamma'+randn(50000,2)*chol(sigmaE); y=samples/(eye(2)-b)'; sampleCov=cov(y); semError=norm(sampleCov-implied,"fro")/norm(implied,"fro");
category=[1;1;2;2;3;3]; response=[1;1.2;2;1.8;.4;.6]; onehot=full(sparse(1:6,category,1,6,3)); beta=onehot\response;
unstable=[0 1.2;1.2 0]; checks=struct("hierarchical_shrinkage_improves",posteriorMse<rawMse, ...
    "sem_stable",max(abs(eig(b)))<1,"sem_covariance_matches",semError<.03, ...
    "unstable_sem_rejected",max(abs(eig(unstable)))>=1,"category_means_recovered",norm(beta-[1.1;1.9;.5])<1e-12);
ids=["DOnSapmaev4","GCuT6NY4ins","ry63Qw5S69k"];
topics=["hierarchical_bayes","structural_equation_model","quantification_i"];
lessons=repmat(struct("video_id","","topic","","public_url",""),3,1);
for k=1:3, lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.latent-structural.v1","lessons",lessons, ...
    "seed",seed,"units",struct("latent","1","response","Pa"),"raw_group_mse",rawMse, ...
    "posterior_group_mse",posteriorMse,"sem_covariance_relative_error",semError, ...
    "category_coefficients",beta',"latent_result_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
