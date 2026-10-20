function H = planarChargeField(Xq, Q, P)
%PLANARCHARGEFIELD Evaluate a 2D planar charge-cloud H field.
H = radia.internal.callMex('radia.PlanarChargeField', double(Xq), double(Q), double(P));
end
