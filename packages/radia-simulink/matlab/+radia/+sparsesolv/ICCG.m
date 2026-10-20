function [solution, info] = ICCG(matrix, rhs, options)
%ICCG Native IC-preconditioned CG solve of matrix * x = rhs.
%   [x, info] = radia.sparsesolv.ICCG(A, b) solves on the free DOFs of the
%   matrix's FE space with the shared SparseSolv ICCG contract:
%     - diagonal scaling (S A S) y = S b, x = S y, when DiagonalScaling;
%     - stop when the recursive relative residual of that system is below
%       Tolerance (0 MaxIterations means 2*n);
%     - an exactly zero rhs returns x = 0; an initial guess that already
%       meets Tolerance returns without building the IC factor;
%     - IC shift alpha >= 1 on rows with Re(a_ii) > 0; AutoShift adds 0.01
%       while a pivot has Re(d) < 1e-6 |a_ii| and alpha < 5, and fails at
%       the limit; a zero or non-finite pivot is an error;
%     - stagnation stop when more than DivergenceCount iterations neither
%       set a new best nor stay below best*DivergenceThreshold;
%     - SaveBestResult returns the best iterate, initial guess included.
%   x is a new radia.ngsolve.Vector. info has converged, iterations,
%   best_iteration, final_residual (scaled-system recursive residual),
%   true_residual (||b - A x|| / ||b|| on the original free-DOF system),
%   actual_shift (0 when no IC factor was built) and residual_history
%   (empty unless SaveResidualHistory). Defaults equal the Python
%   SparseSolvSolver defaults.
%   Non-convergence is reported in info, not raised; invalid options,
%   non-finite input and IC breakdown raise radia:mex:Exception.
%   Conjugate=true selects Hermitian positive-definite ICCG (L D L^H).
%   The default false retains complex-symmetric products (L D L^T).
%   Hermitian structure and positive diagonal are checked; the caller must
%   ensure positive definiteness. Encountered non-positive curvature raises.
arguments
    matrix (1,1) radia.ngsolve.Matrix
    rhs (1,1) radia.ngsolve.Vector
    options.InitialGuess = []
    options.Tolerance (1,1) double {mustBeFinite,mustBePositive} = 1e-8
    options.MaxIterations (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.Shift (1,1) double {mustBeFinite,mustBeGreaterThanOrEqual(options.Shift,1)} = 1
    options.Conjugate (1,1) logical = false
    options.AutoShift (1,1) logical = true
    options.DiagonalScaling (1,1) logical = true
    options.SaveBestResult (1,1) logical = true
    options.SaveResidualHistory (1,1) logical = false
    options.DivergenceCheck (1,1) logical = true
    options.DivergenceThreshold (1,1) double {mustBeFinite,mustBePositive} = 10
    options.DivergenceCount (1,1) double {mustBeInteger,mustBeNonnegative} = 10
    options.UseABMC (1,1) logical = false
    options.ABMCBlockSize (1,1) double {mustBeInteger,mustBePositive} = 4
    options.ABMCNumColors (1,1) double {mustBeInteger,mustBePositive} = 4
    options.ABMCReorderSpMV (1,1) logical = false
    options.ABMCUseRCM (1,1) logical = false
end
if isempty(options.InitialGuess)
    guess = uint64(0);
else
    if ~isa(options.InitialGuess, 'radia.ngsolve.Vector') || ~isscalar(options.InitialGuess)
        error('radia:sparsesolv:InvalidInitialGuess', ...
            'InitialGuess must be empty or a scalar radia.ngsolve.Vector');
    end
    guess = options.InitialGuess.nativeHandle();
end
native = struct( ...
    'tolerance', options.Tolerance, ...
    'max_iterations', options.MaxIterations, ...
    'shift', options.Shift, ...
    'conjugate', options.Conjugate, ...
    'auto_shift', options.AutoShift, ...
    'diagonal_scaling', options.DiagonalScaling, ...
    'save_best_result', options.SaveBestResult, ...
    'save_residual_history', options.SaveResidualHistory, ...
    'divergence_check', options.DivergenceCheck, ...
    'divergence_threshold', options.DivergenceThreshold, ...
    'divergence_count', options.DivergenceCount, ...
    'use_abmc', options.UseABMC, ...
    'abmc_block_size', options.ABMCBlockSize, ...
    'abmc_num_colors', options.ABMCNumColors, ...
    'abmc_reorder_spmv', options.ABMCReorderSpMV, ...
    'abmc_use_rcm', options.ABMCUseRCM);
[h, info] = radia.internal.callMex('sparsesolv.iccg', matrix.nativeHandle(), ...
    rhs.nativeHandle(), guess, native);
solution = radia.ngsolve.Vector.fromNativeHandle(h);
end
