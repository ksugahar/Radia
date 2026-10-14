function result = affineCellSelfEnergyShapeDerivative(cellType, nodes, velocities)
%AFFINECELLSELFENERGYSHAPEDERIVATIVE Analytic TET/HEX/WEDGE self derivative.
arguments
    cellType (1,1) string
    nodes (:,3) double
    velocities (:,:,3) double
end
cellType = lower(cellType);
validTypes = ["tet", "hex", "wedge"];
if ~any(cellType == validTypes)
    error("radia:topopt:CellType", ...
        "cellType must be tet, hex, or wedge.");
end
expectedNodes = [4, 8, 6];
expected = expectedNodes(cellType == validTypes);
if size(nodes,1) ~= expected || size(velocities,2) ~= expected
    error("radia:topopt:CellShape", ...
        "nodes and velocities do not match the selected element family.");
end
result = radia.internal.callMex( ...
    'hdiv.affine_cell_self_energy_shape_derivative', ...
    char(cellType), nodes, velocities);
end
