function q=affineHDivQuadrature(order)
% Existing CPU near rule: symmetric degree 5 for p1/p2, Gauss-Duffy5 for p3.
arguments
    order (1,1) double {mustBeMember(order,[1 2 3])}
end
q=radia.gpu.curvedHDivQuadrature(1);
if order==3
    n=5;k=(1:n-1)';beta=k./sqrt(4*k.^2-1);
    [V,D]=eig(diag(beta,1)+diag(beta,-1));[x,ix]=sort(diag(D));x=(x+1)/2;w=V(1,ix)'.^2;
    tet=zeros(n^3,4);tri=zeros(n^2,4);t=0;s=0;
    for a=1:n
        for b=1:n
            u=x(a);v=x(b);s=s+1;tri(s,:)=[u,v*(1-u),0,w(a)*w(b)*(1-u)];
            for c=1:n
                t=t+1;tet(t,:)=[u,v*(1-u),x(c)*(1-u)*(1-v),w(a)*w(b)*w(c)*(1-u)^2*(1-v)];
            end
        end
    end
    q.tet_rule=tet;q.tri_rule=tri;
end
% CPU affine Gram monomials use mesh-vertex coordinates. Native P2 samples
% use NGSolve coordinates; map the integration nodes, preserving their order.
g=q.tet_rule(:,1:3);q.tet_rule(:,1:3)=[1-sum(g,2),g(:,1:2)];
g=q.tri_rule(:,1:2);q.tri_rule(:,1:2)=[1-sum(g,2),g(:,1)];
end
