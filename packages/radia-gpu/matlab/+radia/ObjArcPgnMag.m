function object = ObjArcPgnMag(center, axis, vertices, phiRange, nseg, symNoSym, magnetization)
%OBJARCPGNMAG Create a finite-length arc magnet with polygonal cross-section.

if nargin < 6
    symNoSym = "nosym";
end
if nargin < 7
    magnetization = zeros(1, 3);
end
object = radia.internal.callMex('radia.ObjArcPgnMag', double(center), ...
    char(string(axis)), double(vertices), double(phiRange), double(nseg), ...
    char(string(symNoSym)), double(magnetization));
end
