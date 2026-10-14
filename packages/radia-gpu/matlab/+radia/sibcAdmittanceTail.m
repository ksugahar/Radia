function y = sibcAdmittanceTail(s, surfaceMeasure, sigma, mu)
%SIBCADMITTANCETAIL Return the unresolved SIBC admittance contribution.

y = radia.internal.callMex( ...
    'hybrid_vim.sibc_admittance_tail', double(s), surfaceMeasure, sigma, mu);
end
