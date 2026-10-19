function shape = FldFrcShpRtg(center, dimensions)
%FLDFRCSHPRTG Create a rectangular Maxwell-tensor force surface.

shape = radia.internal.callMex('radia.FldFrcShpRtg', double(center), double(dimensions));
end
