function test_affine_hdiv_quadrature()
% CPU P2 near-rule degree and simplex moments, including transformed ordering.
q1=radia.gpu.affineHDivQuadrature(1);
assert(size(q1.tet_rule,1)==15);
q3=radia.gpu.affineHDivQuadrature(3);
assert(size(q3.tet_rule,1)==125);
q=radia.gpu.affineHDivQuadrature(2);
assert(size(q.tet_rule,1)==81,'P2 must use the symmetric degree-10 CPU rule');
for i=0:10
    for j=0:10-i
        expected=factorial(i)*factorial(j)/factorial(i+j+2);
        actual=sum(q.tri_rule(:,4).*q.tri_rule(:,1).^i.*q.tri_rule(:,2).^j);
        assert(abs(actual-expected)<=1e-15+2e-13*expected);
        for k=0:10-i-j
            expected=factorial(i)*factorial(j)*factorial(k)/factorial(i+j+k+3);
            actual=sum(q.tet_rule(:,4).*q.tet_rule(:,1).^i.*q.tet_rule(:,2).^j.*q.tet_rule(:,3).^k);
            assert(abs(actual-expected)<=1e-15+2e-13*expected);
        end
    end
end
end
