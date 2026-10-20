function result = Solve(object, precision, maxIter, method, image)
%SOLVE Solve a Radia magnetostatic interaction problem.

if nargin < 4 || isempty(method)
    method = 0;
end
if nargin < 5
    image = "";
end
result = radia.internal.callMex( ...
    'radia.Solve', double(object), double(precision), double(maxIter), ...
    double(method), char(string(image)));
end
