classdef test_ngsolve_sparsecholesky < matlab.unittest.TestCase
    % Explicit native factorization, numerical parity, and repeated RHS reuse.
    properties (TestParameter)
        scalar = struct('real', 2.5, 'complex', 2.5 + 0.5i)
        order = {2, 4}
    end
    methods (TestClassSetup)
        function setupPath(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'matlab')));
            radia.setup(Force=true);
        end
    end
    methods (Test)
        function inverseResidualAndParity(testCase, scalar, order)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            mesh = radia.ngsolve.Mesh.create(fullfile(root,'tests','fixtures','beam','affine_field_tetra.vol'));
            testCase.addTeardown(@() delete(mesh));
            space = radia.ngsolve.FESpace.create(mesh, 'h1', order, Complex=~isreal(scalar));
            testCase.addTeardown(@() delete(space));
            coefficient = radia.ngsolve.CoefficientFunction.constant(scalar);
            testCase.addTeardown(@() delete(coefficient));
            form = radia.ngsolve.BilinearForm.createFromCoefficient(space, 'mass', coefficient);
            testCase.addTeardown(@() delete(form));
            matrix = form.matrix();
            testCase.addTeardown(@() delete(matrix));
            A = matrix.sparse();
            inverse = matrix.inverse();
            testCase.addTeardown(@() delete(inverse));
            testCase.verifyTrue(startsWith(string(inverse.info().kind),'sparsecholesky('));
            rhs = matrix.vector();
            testCase.addTeardown(@() delete(rhs));
            b = scalar * (1:matrix.Rows).';
            rhs.setValues(b);
            x = inverse.matvec(rhs);
            testCase.addTeardown(@() delete(x));
            testCase.verifyLessThan(norm(A*x.values()-b)/norm(b),1e-10);
            testCase.verifyLessThan(norm(x.values()-A\b)/norm(A\b),1e-9);
            rhs.setValues(flipud(b));
            y = inverse.matvec(rhs);
            testCase.addTeardown(@() delete(y));
            testCase.verifyLessThan(norm(A*y.values()-flipud(b))/norm(b),1e-10);
        end
    end
end
