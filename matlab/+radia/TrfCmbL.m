function transform = TrfCmbL(original, added)
%TRFCMBL Compose a transform on the left.

transform = radia.internal.callMex( ...
    'radia.TrfCmbL', double(original), double(added));
end
