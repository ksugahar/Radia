function statesOut = MatHysCommitBatch(material, B, states)
%MATHYSCOMMITBATCH Commit Play hysteresis states for many rows.
statesOut = radia.internal.callMex('radia.MatHysCommitBatch', double(material), ...
    double(B), double(states));
end
