function H = MatHysForwardBatch(material, B, states)
%MATHYSFORWARDBATCH Evaluate Play hysteresis response for many rows.
H = radia.internal.callMex('radia.MatHysForwardBatch', double(material), ...
    double(B), double(states));
end
