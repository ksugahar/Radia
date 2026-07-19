function transform = TrfTrsl(vector)
%TRFTRSL Create a translation transform.

transform = radia.internal.callMex('radia.TrfTrsl', double(vector));
end
