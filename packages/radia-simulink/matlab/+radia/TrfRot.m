function transform = TrfRot(point, axis, angle)
%TRFROT Create a rotation transform.

transform = radia.internal.callMex( ...
    'radia.TrfRot', double(point), double(axis), double(angle));
end
