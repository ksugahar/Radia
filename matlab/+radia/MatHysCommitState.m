function MatHysCommitState(material)
%MATHYSCOMMITSTATE Commit a hysteresis material state.

radia.internal.callMex('radia.MatHysCommitState', double(material));
end
