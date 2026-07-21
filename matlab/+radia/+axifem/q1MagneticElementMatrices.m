function result = q1MagneticElementMatrices(ra, rb, za, zb, permeability, conductivity)
%Q1MAGNETICELEMENTMATRICES Native Q1 Henrotte magnetic element matrices.
% Node order is (ra,za), (rb,za), (rb,zb), (ra,zb). The returned
% stiffness and sigma_mass fields use nodal A_phi (V-DOF) values.
arguments
    ra (1,1) double {mustBeFinite,mustBeNonnegative}
    rb (1,1) double {mustBeFinite,mustBePositive}
    za (1,1) double {mustBeFinite}
    zb (1,1) double {mustBeFinite}
    permeability (1,1) double {mustBeFinite,mustBePositive}
    conductivity (1,1) double {mustBeFinite,mustBeNonnegative}
end
if rb <= ra
    error("radia:axifem:RadiusOrder", "rb must be greater than ra.");
end
if zb <= za
    error("radia:axifem:AxialOrder", "zb must be greater than za.");
end
result = radia.internal.callMex("axifem.q1_magnetic_element_matrices", ...
    ra, rb, za, zb, permeability, conductivity);
end
