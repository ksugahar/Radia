function report = nice_flow_gate(seed)
%NICE_FLOW_GATE Reversible additive coupling for CAE design distributions.
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 37
end
rng(seed,"twister"); n=320;
xa=randn(n,2); wTrue=[.8 -.35; .25 .65]; xb=xa*wTrue+.18*randn(n,2);
x=[xa xb]; w=xa\xb; z=[xa xb-xa*w]; xRound=[z(:,1:2) z(:,3:4)+z(:,1:2)*w];
roundtrip=max(abs(xRound-x),[],"all");

% Additive coupling Jacobian is triangular with determinant one.
j=[eye(2) zeros(2); -w' eye(2)]; detAnalytic=det(j);
h=1e-6; x0=x(1,:)'; jfd=zeros(4);
for k=1:4
    e=zeros(4,1); e(k)=1;
    jfd(:,k)=(forwardMap(x0+h*e,w)-forwardMap(x0-h*e,w))/(2*h);
end

nllBefore=diagonalGaussianNll(x); nllAfter=diagonalGaussianNll(z);
zs=mean(z,1)+randn(64,4).*std(z,0,1);
% Keep an explicit OOD negative control; random sampling alone is not a
% deterministic rejection test.
zs(end,:)=mean(z,1)+6*std(z,0,1);
xs=[zs(:,1:2) zs(:,3:4)+zs(:,1:2)*w];
mahal=sum(((zs-mean(z,1))./std(z,0,1)).^2,2); accepted=mahal<=13.28; % chi2(4), ~0.99

qDirect=forwardQoi(xs); zReplay=[xs(:,1:2) xs(:,3:4)-xs(:,1:2)*w];
xReplay=[zReplay(:,1:2) zReplay(:,3:4)+zReplay(:,1:2)*w];
qReplay=forwardQoi(xReplay); replayErr=max(abs(qDirect-qReplay));
checks=struct("roundtrip",roundtrip<1e-12,"jacobian_unit",abs(detAnalytic-1)<1e-12, ...
    "jacobian_fd",norm(j-jfd,"fro")<1e-8,"likelihood_improves",nllAfter<nllBefore, ...
    "ood_filter_active",any(~accepted)&&any(accepted),"forward_replay",replayErr<1e-12);
report=struct("schema","radia-mcp.matlab.nice-flow-gate.v1", ...
    "lesson",struct("video_id","EiUx0AyyDMU","public_url", ...
    "https://www.youtube.com/watch?v=EiUx0AyyDMU","principle", ...
    "invertible coupling with controlled Jacobian and maximum likelihood"), ...
    "seed",seed,"units",struct("design","1","qoi","Pa^2"), ...
    "coupling_matrix",w,"roundtrip_max_abs",roundtrip, ...
    "jacobian_determinant",detAnalytic,"jacobian_fd_error",norm(j-jfd,"fro"), ...
    "mean_nll_before",nllBefore,"mean_nll_after",nllAfter, ...
    "candidate_count",size(xs,1),"accepted_count",sum(accepted), ...
    "forward_replay_error",replayErr,"generated_candidate_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end

function z=forwardMap(x,w), z=[x(1:2);x(3:4)-w'*x(1:2)]; end
function v=diagonalGaussianNll(x)
s=max(std(x,0,1),eps); c=x-mean(x,1); v=mean(.5*sum((c./s).^2+log(2*pi*s.^2),2));
end
function q=forwardQoi(x), q=sum((x-[.2 -.1 .3 .4]).^2,2); end
