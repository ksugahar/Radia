function varargout = callMex(varargin)
%CALLMEX Reject the retired CPU Radia gateway, including stale installed copies.
error('radia:mex:Retired', ...
    'CPU radia_mex has been removed. Use explicit Python CPU workflows or supported radia-gpu mexcuda operations.');
end
