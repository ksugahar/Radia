function handle = BuildMatrix(object, image)
%BUILDMATRIX Build and cache a Radia interaction matrix.

if nargin < 2
    image = "";
end
handle = radia.internal.callMex('radia.BuildMatrix', double(object), char(string(image)));
end
