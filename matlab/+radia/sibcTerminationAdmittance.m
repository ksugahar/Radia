function y = sibcTerminationAdmittance(s, kSibc, d)
%SIBCTERMINATIONADMITTANCE Return the Schur termination admittance.

y = radia.internal.callMex( ...
    'hybrid_vim.sibc_termination_admittance', double(s), kSibc, d);
end
