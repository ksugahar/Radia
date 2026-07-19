function report = energy_model_gate(seed)
%ENERGY_MODEL_GATE Exact small-RBM checks for Boltzmann/DBN teaching.
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 61
end
rng(seed,"twister"); nv=3; nh=2;
v=dec2bin(0:2^nv-1)-'0'; h=dec2bin(0:2^nh-1)-'0';
w=.12*randn(nv,nh); bv=zeros(1,nv); bh=zeros(1,nh);
data=[1 1 0;1 1 0;1 0 1;1 1 0;0 0 1;1 0 1];
[ll,grad,probV,posterior]=exactRbm(w,bv,bh,v,h,data);
step=.08; [llAfter,~,~,~]=exactRbm(w+step*grad,bv,bh,v,h,data);

epsv=1e-6; wp=w; wm=w; wp(1,1)=wp(1,1)+epsv; wm(1,1)=wm(1,1)-epsv;
lp=exactRbm(wp,bv,bh,v,h,data); lm=exactRbm(wm,bv,bh,v,h,data);
fd=(lp-lm)/(2*epsv); gradientError=abs(fd-grad(1,1))/max(abs(fd),1e-12);

hiddenFeatures=posterior(dataToIndex(data),:);
wLayer2=randn(nh,1)*.1; layer2Score=hiddenFeatures*wLayer2;
[~,modeIdx]=max(probV); candidate=v(modeIdx,:); solverQoi=sum((candidate-[1 1 0]).^2);
checks=struct("probability_normalized",abs(sum(probV)-1)<1e-12, ...
    "probability_nonnegative",all(probV>=0),"gradient_matches_fd",gradientError<1e-6, ...
    "learning_step_improves_likelihood",llAfter>ll, ...
    "dbn_handoff_shape",isequal(size(hiddenFeatures),[size(data,1) nh]), ...
    "second_layer_finite",all(isfinite(layer2Score)), ...
    "forward_candidate_finite",isfinite(solverQoi));

report=struct("schema","radia-mcp.matlab.energy-model.v1", ...
    "lessons",[struct("video_id","FW2tEPyMc7A","topic","boltzmann_overview", ...
    "public_url","https://www.youtube.com/watch?v=FW2tEPyMc7A"); ...
    struct("video_id","Yv5NQMLXZe4","topic","boltzmann_learning", ...
    "public_url","https://www.youtube.com/watch?v=Yv5NQMLXZe4"); ...
    struct("video_id","ZWhOR26jge0","topic","deep_belief_network", ...
    "public_url","https://www.youtube.com/watch?v=ZWhOR26jge0")], ...
    "seed",seed,"units",struct("energy","1","solver_qoi","Pa^2"), ...
    "log_likelihood_before",ll,"log_likelihood_after",llAfter, ...
    "gradient_fd_relative_error",gradientError,"visible_probability",probV', ...
    "dbn_hidden_feature_shape",size(hiddenFeatures),"solver_candidate",candidate, ...
    "solver_qoi",solverQoi,"generated_candidate_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end

function [ll,grad,pv,post]=exactRbm(w,bv,bh,v,h,data)
nv=size(v,1); nhs=size(h,1); joint=zeros(nv,nhs);
for i=1:nv
    for j=1:nhs, joint(i,j)=exp(bv*v(i,:)'+bh*h(j,:)'+v(i,:)*w*h(j,:)'); end
end
joint=joint/sum(joint,"all"); pv=sum(joint,2);
postState=joint./max(pv,eps);
% Convert probabilities over the 2^nh hidden states into expectations of
% the nh hidden units before using them as DBN features/positive phase.
post=postState*h;
idx=dataToIndex(data); ll=mean(log(max(pv(idx),realmin)));
dataVh=zeros(size(w));
for k=1:size(data,1), dataVh=dataVh+data(k,:)'*post(idx(k),:); end
dataVh=dataVh/size(data,1); modelVh=zeros(size(w));
for i=1:nv, for j=1:nhs, modelVh=modelVh+joint(i,j)*(v(i,:)'*h(j,:)); end, end
grad=dataVh-modelVh;
end
function idx=dataToIndex(data), idx=1+data*(2.^(size(data,2)-1:-1:0))'; end
