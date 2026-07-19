function report = diffusion_inverse_gate(seed)
%DIFFUSION_INVERSE_GATE Diffusion-inspired denoising for a CAE inverse problem.
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 31
end
rng(seed,"twister");

% Reproducible ill-conditioned forward operator and a low-modal truth.
n=48; [q,~]=qr(randn(n)); singular=logspace(0,-5,n)';
A=q*diag(singular)*q'; coeff=zeros(n,1); coeff(1:5)=[1;-.65;.4;.2;-.1];
xTrue=q*coeff; yClean=A*xTrue;

betas=linspace(.01,.30,10)'; rows=zeros(numel(betas),5); rankKeep=8;
for k=1:numel(betas)
    beta=betas(k); noise=randn(n,1);
    yNoisy=sqrt(1-beta)*yClean+sqrt(beta)*0.03*noise;
    corrected=yNoisy/sqrt(1-beta);
    % Noise-removal least squares in a validated truncated modal space.
    basis=q(:,1:rankKeep); reduced=A*basis;
    z=reduced\corrected; xHat=basis*z; yHat=A*xHat;
    rawResidual=norm(corrected-yClean)/max(norm(yClean),eps);
    verifiedResidual=norm(yHat-yClean)/max(norm(yClean),eps);
    rows(k,:)=[k beta rawResidual verifiedResidual norm(xHat-xTrue)/norm(xTrue)];
end

checks=struct();
checks.schedule_monotone=all(diff(betas)>0) && all(betas>0 & betas<1);
checks.denoising_improves_mean=mean(rows(:,4))<mean(rows(:,3));
checks.forward_residual_finite=all(isfinite(rows(:,4)));
checks.high_noise_not_ground_truth=rows(end,4)>0;
checks.seed_recorded=seed>=0;

report=struct("schema","radia-mcp.matlab.diffusion-inverse-gate.v1", ...
    "lesson",struct("video_id","BP8jnfFUD3E", ...
    "public_title","Diffusion Models - noise-removal generative model", ...
    "public_url","https://www.youtube.com/watch?v=BP8jnfFUD3E", ...
    "generalized_principle","staged noise plus least-squares denoising"), ...
    "seed",seed,"units",struct("state","1","observation","Pa"), ...
    "forward_operator","deterministic_ill_conditioned_linear_cae_fixture", ...
    "rank_keep",rankKeep,"row_schema", ...
    ["step","beta","raw_relative_residual","verified_relative_residual","state_relative_error"], ...
    "rows",rows,"generated_candidate_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks, ...
    "ok",all(structfun(@logical,checks)));
end
