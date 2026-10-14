function [matrix, dof] = GetInteractMatrix(handle)
%GETINTERACTMATRIX Return a cached Radia interaction matrix.

[matrix, dof] = radia.internal.callMex('radia.GetInteractMatrix', double(handle));
end
