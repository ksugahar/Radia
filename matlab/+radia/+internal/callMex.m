function varargout = callMex(varargin)
%CALLMEX Enter the checked native MEX boundary.

radia.setup();
[varargout{1:nargout}] = radia_mex(varargin{:});
end
