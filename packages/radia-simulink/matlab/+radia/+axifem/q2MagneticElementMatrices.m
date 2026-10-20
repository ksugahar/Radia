function result = q2MagneticElementMatrices(ra, rb, za, zb, mu, sigma)
%Q2MAGNETICELEMENTMATRICES Native Q2 axisymmetric magnetic element matrices.
%   RESULT = radia.axifem.q2MagneticElementMatrices(RA,RB,ZA,ZB,MU,SIGMA)
%   returns the 9-by-9 Henrotte V-DOF stiffness and conductivity mass
%   matrices. Axis-touching elements retain nine local entries; the three
%   r=0 rows and columns are exactly zero and must be constrained globally.

arguments
    ra (1,1) double {mustBeFinite,mustBeNonnegative}
    rb (1,1) double {mustBeFinite,mustBePositive}
    za (1,1) double {mustBeFinite}
    zb (1,1) double {mustBeFinite}
    mu (1,1) double {mustBeFinite,mustBePositive}
    sigma (1,1) double {mustBeFinite,mustBeNonnegative}
end
result = radia.internal.callMex( ...
    'axifem.q2_magnetic_element_matrices', ra, rb, za, zb, mu, sigma);
end
