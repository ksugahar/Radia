classdef test_hacapk_rank_cap < matlab.unittest.TestCase
    % A compression rank budget must not relax the requested kernel accuracy.
    properties (TestParameter)
        rankCap = {1, 2}
    end

    methods (TestClassSetup)
        function setupNative(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'matlab')));
            radia.setup(Force=true);
        end
    end

    methods (Test)
        function cappedMatrixMatchesLaplaceKernel(testCase, rankCap)
            index = (1:64)';
            points = [sin(sqrt(2)*index), cos(sqrt(3)*index), ...
                sin(sqrt(5)*index)];
            weights = linspace(0.5, 1.5, 64)';
            epsilon = 1e-3;
            gram = radia.HACApKChargeGram.from_sampled_laplace( ...
                points, weights, epsilon, Build=false);
            testCase.addTeardown(@() delete(gram));

            built = gram.build(AcaEps=1e-14, LeafSize=4, Eta=2, ...
                MaxRank=rankCap);
            basis = eye(64);
            columns = arrayfun(@(j) gram.matvecSym(basis(:, j)), ...
                1:64, UniformOutput=false);
            actual = cat(2, columns{:});
            delta = reshape(points, 64, 1, 3) - reshape(points, 1, 64, 3);
            expected = (weights*weights') ./ ...
                (4*pi*sqrt(sum(delta.^2, 3) + epsilon^2));
            stats = gram.info();

            testCase.verifyTrue(built);
            testCase.verifyLessThanOrEqual(stats.max_rank, rankCap);
            testCase.verifyEqual(actual, expected, ...
                RelTol=2e-12, AbsTol=1e-14);
        end
    end
end
