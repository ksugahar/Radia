function force = planarMaxwellForceCircle(Xq, Q, Rc, Hext, center, n)
%PLANARMAXWELLFORCECIRCLE Maxwell-stress force of a planar charge cloud.
if nargin < 4 || isempty(Hext), Hext = [0, 0]; end
if nargin < 5 || isempty(center), center = [0, 0]; end
if nargin < 6 || isempty(n), n = 1440; end
force = radia.internal.callMex('radia.PlanarMaxwellForceCircle', ...
    double(Xq), double(Q), double(Rc), double(center(1)), double(center(2)), ...
    double(n), double(Hext(1)), double(Hext(2)));
end
