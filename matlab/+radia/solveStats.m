function stats = solveStats()
%SOLVESTATS Return statistics from the most recent native solve.

stats = radia.internal.callMex('radia.GetSolveStats');
end
