function varargout = callMex(varargin)
%CALLMEX Enter the checked native MEX boundary.

radia.setup();
if ~isempty(varargin) && isstring(varargin{1}) && isscalar(varargin{1})
    varargin{1} = char(varargin{1});
end
[varargout{1:nargout}] = radia_mex(varargin{:});
end
