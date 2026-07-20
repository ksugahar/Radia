function state = MatHysSaveState(material)
%MATHYSSAVESTATE Save a hysteresis material state.

state = radia.internal.callMex('radia.MatHysSaveState', double(material));
end
