function MatHysRestoreState(material, state)
%MATHYSRESTORESTATE Restore a hysteresis material state.

radia.internal.callMex('radia.MatHysRestoreState', double(material), double(state));
end
