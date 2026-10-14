function names = commands()
%COMMANDS Return the native commands implemented by radia_mex.

names = radia.internal.callMex('api.commands');
end
