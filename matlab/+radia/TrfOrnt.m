function object = TrfOrnt(source, transform)
%TRFORNT Apply a transform to a Radia object.

object = radia.internal.callMex( ...
    'radia.TrfOrnt', double(source), double(transform));
end
