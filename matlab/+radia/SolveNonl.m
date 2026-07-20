function result = SolveNonl(object, precision, maxIter, method, nonlMethod, image)
%SOLVENONL Solve a nonlinear magnetostatic problem with a selected method.

if nargin < 6
    image = "";
end
result = radia.internal.callMex('radia.SolveNonl', double(object), ...
    double(precision), double(maxIter), double(method), double(nonlMethod), ...
    char(string(image)));
end
