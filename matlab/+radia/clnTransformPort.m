function reduced = clnTransformPort(Q, port)
%CLNTRANSFORMPORT Transform a loop-space port vector into CLN coordinates.
reduced = radia.internal.callMex('cln.transform_port', double(Q), double(port));
end
