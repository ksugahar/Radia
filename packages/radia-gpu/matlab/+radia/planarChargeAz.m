function Az = planarChargeAz(Xq, Q, P)
%PLANARCHARGEAZ Evaluate the 2D planar charge-cloud A_z potential.
Az = radia.internal.callMex('radia.PlanarChargeAz', double(Xq), double(Q), double(P));
end
