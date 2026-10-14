function report = actor_critic_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 83
end
rng(seed,"twister"); theta=[.2;-.1;.05]; reward=[-.4;.8;.2]; p=softmax(theta); expected=p'*reward;
grad=p.*(reward-expected); h=1e-6; fd=zeros(3,1);
for k=1:3, e=zeros(3,1); e(k)=1; fd(k)=(softmax(theta+h*e)'*reward-softmax(theta-h*e)'*reward)/(2*h); end
gradErr=norm(grad-fd)/norm(fd);
n=5000; u=rand(n,1); acts=1+(u>p(1))+(u>p(1)+p(2));
returns=reward(acts)+.3*randn(n,1); score=repmat(-p',n,1); score=subsSet(score,acts);
raw=score.*returns; baseline=expected; centered=score.*(returns-baseline);
varianceRatio=mean(var(centered,0,1))/mean(var(raw,0,1));
gamma=.95; lambda=.8; trace=zeros(3,1); traces=zeros(3,20);
for t=1:20, trace=gamma*lambda*trace+grad; traces(:,t)=trace; end
dqnTarget=.5+gamma*max([.1 .4 .2]);
checks=struct("policy_normalized",abs(sum(p)-1)<1e-12,"gradient_matches_fd",gradErr<1e-8, ...
    "baseline_reduces_variance",varianceRatio<1,"eligibility_finite",all(isfinite(traces),"all"), ...
    "dqn_target_finite",isfinite(dqnTarget));
ids=["tqvOn1JZRXM","hE2A-IxHWc4","pXJigtsFiaA","JjSlfur0GvI","H_cx9IZseK0","0KzrjmlJ42A","-IRybgttzMo", ...
    "qyWev-rj5gg","2k5y3blCMbw","qyzfDt84OWM","LdL3UrYSDLw","E1r0vkv4lkk","uvxmmIWnJS4","TfoOugDRajM"];
topics=["dpg_theorem","dqn","actor_critic_trace","actor_critic","reinforce","policy_gradient_infinite","policy_gradient_proof", ...
    "policy_gradient_overview","deep_rl_intro","eligibility_trace_rta","eligibility_trace_derivation","forward_backward_trace","backward_td_lambda","td_lambda"];
lessons=repmat(struct("video_id","","topic","","public_url",""),numel(ids),1);
for k=1:numel(ids), lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.actor-critic.v1","lessons",lessons, ...
    "seed",seed,"units",struct("reward","-Pa^2","gradient","1"),"policy",p', ...
    "gradient_relative_error",gradErr,"baseline_variance_ratio",varianceRatio, ...
    "dqn_target",dqnTarget,"policy_candidate_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
function p=softmax(x), e=exp(x-max(x)); p=e/sum(e); end
function s=subsSet(s,acts)
for i=1:numel(acts), s(i,acts(i))=s(i,acts(i))+1; end
end
