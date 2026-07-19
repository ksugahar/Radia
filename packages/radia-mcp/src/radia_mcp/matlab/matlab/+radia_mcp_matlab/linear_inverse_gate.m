function report = linear_inverse_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 79
end
rng(seed,"twister"); a=randn(8,4); xTrue=[1;-.5;.2;.8]; b=a*xTrue+.01*randn(8,1); ap=pinv(a); x=ap*b;
mp=[norm(a*ap*a-a,"fro"),norm(ap*a*ap-ap,"fro"),norm((a*ap)'-(a*ap),"fro"),norm((ap*a)'-(ap*a),"fro")];
normalResidual=norm(a'*(a*x-b));
[u,s,v]=svd(a,"econ"); r=3; xt=v(:,1:r)*(s(1:r,1:r)\(u(:,1:r)'*b)); tsvdResidual=norm(a*xt-b);

xdata=randn(300,3); ydata=xdata*[.8 -.2;.1 .7;-.4 .3]+.1*randn(300,2);
cx=cov(xdata)+1e-9*eye(3); cy=cov(ydata)+1e-9*eye(2); cxy=cov([xdata ydata]); cxy=cxy(1:3,4:5);
[wx,dx]=eig(cx); invsx=wx*diag(1./sqrt(diag(dx)))*wx'; [wy,dy]=eig(cy); invsy=wy*diag(1./sqrt(diag(dy)))*wy';
sv=svd(invsx*cxy*invsy); ccaMax=sv(1);
checks=struct("moore_penrose_identities",max(mp)<1e-10,"normal_equation",normalResidual<1e-10, ...
    "least_squares_matches_backslash",norm(x-a\b)<1e-10,"tsvd_finite",isfinite(tsvdResidual), ...
    "cca_bounded",ccaMax>=0&&ccaMax<=1+1e-8,"transpose_shape",isequal(size(a'),[4 8]));
ids=["cQrK6maa-fc","EReuJVxYc80","UKSpRjSDB9M","Lm8egVsK2tQ","nQgepN94V6E","Bs6M5iky7a0", ...
    "c0OfA7AjLjA","8mqH8aZdceM","uUpCeCjkhEk","c179m_cO6H4"];
topics=["multiple_regression","pseudoinverse_formula","pseudoinverse_direction","transpose","cca","svd_formula", ...
    "relation_components_2","relation_calculation","relation_decomposition","singular_value_scaling"];
lessons=repmat(struct("video_id","","topic","","public_url",""),numel(ids),1);
for k=1:numel(ids), lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.linear-inverse.v1","lessons",lessons, ...
    "seed",seed,"units",struct("matrix","1","observation","Pa"),"moore_penrose_errors",mp, ...
    "normal_residual",normalResidual,"tsvd_rank",r,"tsvd_residual",tsvdResidual, ...
    "cca_max_correlation",ccaMax,"inverse_candidate_is_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
