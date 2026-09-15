classdef TestSpecleash < matlab.unittest.TestCase
    properties
        Fixtures
    end

    methods (TestClassSetup)
        function setPaths(testCase)
            here = fileparts(mfilename('fullpath'));
            matlabDir = fileparts(here);
            addpath(matlabDir);
            testCase.Fixtures = fullfile(fileparts(matlabDir), 'tests', 'fixtures');
        end
    end

    methods (Test)
        function identicalCurvesAreZero(testCase)
            P = [0 0; 0.5 1; 1 0];
            testCase.verifyEqual(specleash_frechet(P, P), 0);
        end

        function firstPointSeparationCounts(testCase)
            x = linspace(0, 1, 10).';
            a = zeros(10, 1);
            a(1) = 1;
            testCase.verifyEqual(specleash_frechet([x a], [x zeros(10, 1)]), 1, 'AbsTol', 1e-12);
        end

        function eitherCurveMayWait(testCase)
            x = linspace(0, 1, 11).';
            a = double((0:10).' >= 6);
            b = double((0:10).' >= 5);
            testCase.verifyEqual(specleash_frechet([x a], [x b]), 0.1, 'AbsTol', 1e-12);
            testCase.verifyEqual(specleash_frechet([x b], [x a]), 0.1, 'AbsTol', 1e-12);
        end

        function pathIsAValidCoupling(testCase)
            rng(3);
            P = rand(7, 2);
            Q = rand(9, 2);
            [d, path, links] = specleash_frechet(P, Q);
            testCase.verifyEqual(path(1, :), [1 1]);
            testCase.verifyEqual(path(end, :), [7 9]);
            testCase.verifyTrue(all(ismember(diff(path), [1 0; 0 1; 1 1], 'rows')));
            testCase.verifyEqual(max(links), d, 'AbsTol', 1e-12);
        end

        function matchesPythonPairs(testCase)
            e = jsondecode(fileread(fullfile(testCase.Fixtures, 'expected.json')));
            for k = 1:numel(e.pairs)
                p = e.pairs(k);
                P = readmatrix(fullfile(testCase.Fixtures, p.p), 'FileType', 'text');
                Q = readmatrix(fullfile(testCase.Fixtures, p.q), 'FileType', 'text');
                testCase.verifyEqual(specleash_frechet(P, Q), p.distance, 'AbsTol', 1e-9, p.p);
            end
        end

        function matchesPythonComparisons(testCase)
            e = jsondecode(fileread(fullfile(testCase.Fixtures, 'expected.json')));
            for k = 1:numel(e.comparisons)
                c = e.comparisons(k);
                P = readmatrix(fullfile(testCase.Fixtures, c.a), 'FileType', 'text');
                Q = readmatrix(fullfile(testCase.Fixtures, c.b), 'FileType', 'text');
                d = specleash_compare(P, Q, c.wl_range.', 'NormaliseAt', c.normalise_at, 'XScale', c.x_scale);
                testCase.verifyEqual(d, c.distance, 'AbsTol', 1e-9, c.a);
            end
        end

        function partialCoverageIsRefused(testCase)
            P = [0.5 1; 1.0 1.1; 2.0 1.3];
            Q = [0.8 1; 1.5 1.2; 2.0 1.4];
            testCase.verifyError(@() specleash_compare(P, Q, [0.6 2.0]), 'specleash:coverage');
        end
    end
end
