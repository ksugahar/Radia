function transform = TrfCmbR(original, added)
%TRFCMBR Compose a transform on the right.

transform = radia.internal.callMex( ...
    'radia.TrfCmbR', double(original), double(added));
end
