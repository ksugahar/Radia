function hIrr = MatHysIrreversible(material, B)
%MATHYSIRREVERSIBLE Evaluate the irreversible field contribution.

hIrr = radia.internal.callMex('radia.MatHysIrreversible', double(material), double(B));
end
