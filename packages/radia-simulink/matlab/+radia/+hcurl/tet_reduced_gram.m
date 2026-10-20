function gram = tet_reduced_gram(cell_verts, exponents, coefficients, ...
        n_modes, ref_points, ref_weights)
%TET_REDUCED_GRAM Assemble the reduced HCurl Gram on tetrahedra.
%   Canonical snake_case MATLAB name matching the MEX/Python contract.

gram = radia.tetHCurlReducedGram(cell_verts, exponents, coefficients, ...
    n_modes, ref_points, ref_weights);
end
