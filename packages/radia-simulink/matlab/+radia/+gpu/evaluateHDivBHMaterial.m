function out=evaluateHDivBHMaterial(E,weights,law,m,wantTangent,deviceOutput)
% GPU material energy/gradient/tensor Hessian from NGSolve-owned sparse maps.
% E is a three-cell array of sparse gpuArray mapped evaluations.
if nargin<6,deviceOutput=false;end
validateattributes(deviceOutput,{'logical'},{'scalar'});
assert(isa(E{1},'gpuArray') && isa(weights,'gpuArray'));
nq=numel(weights);if isa(m,'gpuArray'),mDevice=m(:);else,mDevice=gpuArray(m(:));end;v=[E{1}*mDevice,E{2}*mDevice,E{3}*mDevice];
a=sqrt(sum(v.^2,2));clipped=min(a,law.cap);
if isfield(law,'device')
    knots=law.device.knots;c=law.device.coefficients;base=law.device.integral_base;
else
    knots=gpuArray(law.knots);c=gpuArray(law.coefficients);base=gpuArray(law.integral_base);
end
lo=ones(nq,1,'gpuArray');hi=(numel(law.knots)-1)*lo;
for k=1:ceil(log2(numel(law.knots)))+1
    mid=floor((lo+hi+1)/2);take=clipped>=knots(mid);lo=take.*mid+(~take).*lo;hi=take.*hi+(~take).*(mid-1);
end
x=clipped-knots(lo);cc=c(lo,:);
h=((cc(:,1).*x+cc(:,2)).*x+cc(:,3)).*x+cc(:,4);
diffnu=(3*cc(:,1).*x+2*cc(:,2)).*x+cc(:,3);
wco=base(lo)+cc(:,1).*x.^4/4+cc(:,2).*x.^3/3+cc(:,3).*x.^2/2+cc(:,4).*x;
over=a>law.cap;extra=max(a-law.cap,0);h=h+law.barrier_slope*extra;
diffnu=over*law.barrier_slope+(~over).*diffnu;
wco=wco+law.last_field*extra+.5*law.barrier_slope*extra.^2;
sec=h./max(a,realmin);sec=(a>0).*sec+(a==0)*law.initial_reluctivity;
assert(gather(all(isfinite(sec) & isfinite(diffnu) & sec>0 & diffnu>0)),'Invalid inverse BH tangent');
field=v.*sec;gradient=E{1}'*(weights.*field(:,1))+E{2}'*(weights.*field(:,2))+E{3}'*(weights.*field(:,3));
out=struct('gradient',gradient,'energy',sum(weights.*wco),'maximum_magnetization',max(a));
if ~deviceOutput,out.gradient=gather(out.gradient);out.energy=gather(out.energy);out.maximum_magnetization=gather(out.maximum_magnetization);end
if wantTangent
    nd=size(E{1},2);W=gpuArray(sparse(nd,nd));direction=v./max(a,realmin);
    for i=1:3
        for j=i:3
            d=(diffnu-sec).*direction(:,i).*direction(:,j);if i==j,d=d+sec;end
            term=E{i}'*spdiags(weights.*d,0,nq,nq)*E{j};
            if i==j,W=W+term;else,W=W+term+term';end
        end
    end
    % Assemble the declared symmetric bilinear form from its upper triangle.
    % This is not a repair of an externally supplied non-symmetric matrix.
    W=triu(W)+triu(W,1)';if deviceOutput,out.tangent=W;else,out.tangent=gather(W);end
end
end
