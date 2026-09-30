function handle = nextIHNativeHandle()
%NEXTIHNATIVEHANDLE Process-local, non-serializable IH runtime identity.
% Locked state survives ordinary clear functions and radia_mex reloads.
% Never save/transfer handles between MATLAB processes. Explicitly unlocking
% or modifying this internal allocator while handles exist is unsupported.
mlock
persistent counter
if isempty(counter)
    counter = uint64(0);
end
limit = bitshift(uint64(1),63) - uint64(1);
if counter >= limit
    error('radia:ih:HandleExhausted','IH handle identities exhausted in this MATLAB process.');
end
counter = counter + uint64(1);
handle = bitor(bitshift(uint64(1),63),counter);
end
