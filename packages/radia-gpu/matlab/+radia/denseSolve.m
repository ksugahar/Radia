function x = denseSolve(A, b)
%DENSESOLVE Solve a small dense real or complex mixed-Galerkin block.

x = radia.internal.callMex('hybrid_vim.solve', double(A), double(b));
end
