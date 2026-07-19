function z = skinImpedance(s, sigma, mu)
%SKINIMPEDANCE Return the MQS half-space surface impedance.

z = radia.internal.callMex('hybrid_vim.skin_impedance', double(s), sigma, mu);
end
