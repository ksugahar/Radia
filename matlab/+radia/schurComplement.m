function S = schurComplement(Kkk, Kke, Kek, Kee)
%SCHURCOMPLEMENT Eliminate the second block of a mixed-Galerkin operator.

S = radia.internal.callMex( ...
    'hybrid_vim.schur', double(Kkk), double(Kke), double(Kek), double(Kee));
end
