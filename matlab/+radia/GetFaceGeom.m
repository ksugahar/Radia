function [geometry, dof] = GetFaceGeom(handle)
%GETFACEGEOM Return per-DoF face geometry for a cached interaction matrix.

[geometry, dof] = radia.internal.callMex('radia.GetFaceGeom', double(handle));
end
