function z = sibcTerminationImpedance(s, kSibc, d)
%SIBCTERMINATIONIMPEDANCE Return the Schur termination impedance.

z = radia.internal.callMex( ...
    'hybrid_vim.sibc_termination_impedance', double(s), kSibc, d);
end
