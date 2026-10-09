classdef TestGetTformQCMetrics < matlab.unittest.TestCase
    %TESTGETTFORMQCMETRICS Unit tests for the shared transform QC helper.
    %
    %   The helper exists because MATLAB exposes two incompatible matrix
    %   conventions for 2-D transforms. Reading the translation from the wrong
    %   element index silently reports [0 0] for every modern transform class,
    %   which in turn disables the runaway-registration QC guard.

    properties (Constant)
        Scale = 1.2
        RotationDeg = 10
        TranslationXY = [3 -4]
    end

    methods (Access = private)
        function verifyExpectedGeometry(testCase, qc, context)
            testCase.verifyEqual(double(qc.translationXY), ...
                TestGetTformQCMetrics.TranslationXY, 'AbsTol', 1e-9, ...
                sprintf('Wrong translation for %s.', context));
            testCase.verifyEqual(double(qc.rotationDeg), ...
                TestGetTformQCMetrics.RotationDeg, 'AbsTol', 1e-9, ...
                sprintf('Wrong rotation for %s.', context));
            testCase.verifyEqual(double(qc.scaleXY), ...
                TestGetTformQCMetrics.Scale * [1 1], 'AbsTol', 1e-9, ...
                sprintf('Wrong scale for %s.', context));
            testCase.verifyEqual(double(qc.determinant), ...
                TestGetTformQCMetrics.Scale^2, 'AbsTol', 1e-9, ...
                sprintf('Wrong determinant for %s.', context));
        end
    end

    methods (Test)
        function testModernPremultiplyTransform(testCase)
            % simtform2d is what imregtform(..., 'similarity', ...) returns on
            % R2024b. Its translation lives in A(1:2,3), not A(3,1:2).
            tform = simtform2d(TestGetTformQCMetrics.Scale, ...
                TestGetTformQCMetrics.RotationDeg, ...
                TestGetTformQCMetrics.TranslationXY);

            qc = getTformQCMetrics(tform);

            testCase.verifyExpectedGeometry(qc, 'simtform2d');
            testCase.verifyEqual(qc.convention, 'premultiply');
            testCase.verifyEqual(qc.transformType, 'similarity');
        end

        function testLegacyPostmultiplyTransform(testCase)
            % The same geometry expressed as a legacy affine2d must produce
            % identical metrics, including the sign of the rotation.
            modern = simtform2d(TestGetTformQCMetrics.Scale, ...
                TestGetTformQCMetrics.RotationDeg, ...
                TestGetTformQCMetrics.TranslationXY);
            legacy = affine2d(modern.A');

            qc = getTformQCMetrics(legacy);

            testCase.verifyExpectedGeometry(qc, 'affine2d');
            testCase.verifyEqual(qc.convention, 'postmultiply');
            testCase.verifyEqual(qc.transformType, 'affine');
        end

        function testConventionsAgree(testCase)
            modern = simtform2d(TestGetTformQCMetrics.Scale, ...
                TestGetTformQCMetrics.RotationDeg, ...
                TestGetTformQCMetrics.TranslationXY);

            qcModern = getTformQCMetrics(modern);
            qcLegacy = getTformQCMetrics(affine2d(modern.A'));

            testCase.verifyEqual(qcLegacy.translationXY, qcModern.translationXY, 'AbsTol', 1e-9);
            testCase.verifyEqual(qcLegacy.rotationDeg, qcModern.rotationDeg, 'AbsTol', 1e-9);
            testCase.verifyEqual(qcLegacy.scaleXY, qcModern.scaleXY, 'AbsTol', 1e-9);
        end

        function testNumericPostmultiplyMatrix(testCase)
            modern = simtform2d(TestGetTformQCMetrics.Scale, ...
                TestGetTformQCMetrics.RotationDeg, ...
                TestGetTformQCMetrics.TranslationXY);

            qc = getTformQCMetrics(modern.A');

            testCase.verifyExpectedGeometry(qc, 'numeric postmultiply matrix');
            testCase.verifyEqual(qc.convention, 'postmultiply');
        end

        function testNumericPremultiplyMatrix(testCase)
            modern = simtform2d(TestGetTformQCMetrics.Scale, ...
                TestGetTformQCMetrics.RotationDeg, ...
                TestGetTformQCMetrics.TranslationXY);

            qc = getTformQCMetrics(modern.A);

            testCase.verifyExpectedGeometry(qc, 'numeric premultiply matrix');
        end

        function testTranslationOnlyTransform(testCase)
            tform = transltform2d([7 11]);

            qc = getTformQCMetrics(tform);

            testCase.verifyEqual(double(qc.translationXY), [7 11], 'AbsTol', 1e-9);
            testCase.verifyEqual(double(qc.rotationDeg), 0, 'AbsTol', 1e-9);
            testCase.verifyEqual(qc.transformType, 'translation');
        end

        function testRigidTransformType(testCase)
            qc = getTformQCMetrics(rigidtform2d(15, [1 2]));

            testCase.verifyEqual(qc.transformType, 'rigid');
            testCase.verifyEqual(double(qc.translationXY), [1 2], 'AbsTol', 1e-9);
            testCase.verifyEqual(double(qc.rotationDeg), 15, 'AbsTol', 1e-9);
        end

        function testUnsupportedInputReturnsNaN(testCase)
            qc = getTformQCMetrics(struct('notATransform', true));

            testCase.verifyTrue(all(isnan(qc.translationXY)));
            testCase.verifyTrue(isnan(qc.rotationDeg));
            testCase.verifyTrue(all(isnan(qc.scaleXY)));
            testCase.verifyTrue(isnan(qc.determinant));
            testCase.verifyEqual(qc.transformType, 'unknown');
        end

        function testNonFiniteMatrixReturnsNaN(testCase)
            qc = getTformQCMetrics([1 0 NaN; 0 1 0; 0 0 1]);

            testCase.verifyTrue(all(isnan(qc.translationXY)));
            testCase.verifyTrue(isnan(qc.determinant));
        end
    end
end
