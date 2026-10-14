function model = makeAxiEddyElementModel(ra, rb, za, zb, mu, sigma, options)
%MAKEAXIEDDYELEMENTMODEL Build a native Q2 axisymmetric eddy state model.
%   The semi-discrete element equation is
%       M_sigma * dA_phi/dt + K * A_phi = source * current(t).
%   Element matrices come from the same C++ kernel used by radia.axifem and
%   radia_mex. The returned exact-ZOH model can run in Simulink without a
%   Python call or finite-element assembly at each sample.

arguments
    ra (1,1) double {mustBeFinite, mustBeNonnegative}
    rb (1,1) double {mustBeFinite, mustBePositive}
    za (1,1) double {mustBeFinite}
    zb (1,1) double {mustBeFinite}
    mu (1,1) double {mustBeFinite, mustBePositive}
    sigma (1,1) double {mustBeFinite, mustBePositive}
    options.SampleTime_s (1,1) double {mustBeFinite, mustBePositive} = 1.0e-5
    options.DirichletDofs double = 1:8
    options.Source double = [0; 0; 0; 0; 0; 0; 0; 0; 1]
    options.InitialState double = []
end

if rb <= ra || zb <= za
    error("radia:simulink:AxiEddyGeometry", ...
        "The element requires 0 <= ra < rb and za < zb.");
end

element = radia.axifem.q2MagneticElementMatrices( ...
    ra, rb, za, zb, mu, sigma);
K = 0.5 * (double(element.stiffness) + double(element.stiffness).');
M = 0.5 * (double(element.sigma_mass) + double(element.sigma_mass).');

dirichlet = unique(double(options.DirichletDofs(:).'), "sorted");
if any(dirichlet ~= fix(dirichlet)) || any(dirichlet < 1) || any(dirichlet > 9)
    error("radia:simulink:AxiEddyBoundary", ...
        "DirichletDofs must contain unique integer Q2 node indices from 1 to 9.");
end
if element.axis_touching && ~all(ismember([1, 4, 8], dirichlet))
    error("radia:simulink:AxiEddyAxis", ...
        "Axis-touching Q2 elements require Dirichlet DOFs 1, 4, and 8.");
end
free = setdiff(1:9, dirichlet, "stable");
if isempty(free)
    error("radia:simulink:AxiEddyBoundary", ...
        "At least one Q2 degree of freedom must remain free.");
end

source = double(options.Source(:));
if numel(source) ~= 9 || any(~isfinite(source))
    error("radia:simulink:AxiEddySource", ...
        "Source must contain nine finite Q2 nodal load values.");
end
Mff = M(free, free);
Kff = K(free, free);
if min(eig(Mff)) <= 0 || rcond(Mff) <= eps
    error("radia:simulink:AxiEddyMass", ...
        "The free conductive mass matrix must be positive definite.");
end

A = -(Mff \ Kff);
B = Mff \ source(free);
C = zeros(9, numel(free));
C(free, :) = eye(numel(free));
D = zeros(9, 1);
if isempty(options.InitialState)
    x0 = zeros(numel(free), 1);
else
    x0 = double(options.InitialState(:));
    if numel(x0) ~= numel(free) || any(~isfinite(x0))
        error("radia:simulink:AxiEddyState", ...
            "InitialState must contain one finite value per free DOF.");
    end
end

augmented = [A, B; zeros(1, size(A, 1) + 1)];
discrete = expm(augmented * options.SampleTime_s);
n = size(A, 1);
Ad = discrete(1:n, 1:n);
Bd = discrete(1:n, n + 1);

model = struct( ...
    "schema", "radia.axifem.q2_eddy.state_space.v1", ...
    "backend", "shared-native-q2-mex", ...
    "element", struct("ra_m", ra, "rb_m", rb, "za_m", za, "zb_m", zb, ...
        "mu_H_per_m", mu, "sigma_S_per_m", sigma, ...
        "axis_touching", logical(element.axis_touching)), ...
    "stiffness", K, "sigma_mass", M, ...
    "dirichlet_dofs", dirichlet, "free_dofs", free, "source", source, ...
    "A", A, "B", B, "C", C, "D", D, ...
    "Ad", Ad, "Bd", Bd, "Cd", C, "Dd", D, ...
    "x0", x0, "sample_time_s", options.SampleTime_s, ...
    "state_order", n, "input_count", 1, "output_count", 9, ...
    "input_convention", "impressed-current amplitude multiplying the nodal source", ...
    "output_convention", "nine Q2 nodal A_phi values; constrained nodes are zero", ...
    "python_per_step", false);
end
