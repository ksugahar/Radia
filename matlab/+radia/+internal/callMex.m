function varargout = callMex(varargin)
%CALLMEX Enter the checked native MEX boundary.

% File-generation settings belong to explicit setup/model initialization.
% Keep runtime checks on every call; do not cache a configured flag here.
radia.setup(ConfigureSimulinkFileGeneration=false);
if ~isempty(varargin) && isstring(varargin{1}) && isscalar(varargin{1})
    varargin{1} = char(varargin{1});
end
[varargout{1:nargout}] = radia_mex(varargin{:});
end
