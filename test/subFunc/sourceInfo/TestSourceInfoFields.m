classdef TestSourceInfoFields < matlab.unittest.TestCase
    %TESTSOURCEINFOFIELDS The closed whitelist of injectable per-data metadata (Phase 6b-2).

    methods (Test)
        function whitelistIsClosedAndNamed(testCase)
            f = sourceInfoFields();
            testCase.verifyEqual({f.field}, {'frameRateHz', 'dimNames', 'exposureMsec'});
            testCase.verifyEqual({f.nvName}, {'FrameRateHz', 'DimNames', 'ExposureMsec'});
        end

        function frameRateRule(testCase)
            v = iRule('frameRateHz');
            testCase.verifyTrue(v(25));
            testCase.verifyFalse(v(NaN));
            testCase.verifyFalse(v([]));
            testCase.verifyFalse(v(0));
            testCase.verifyFalse(v(-1));
            testCase.verifyFalse(v(Inf));
            testCase.verifyFalse(v([1 2]));
            testCase.verifyFalse(v('25'));
        end

        function exposureRule(testCase)
            v = iRule('exposureMsec');
            testCase.verifyTrue(v(0));
            testCase.verifyTrue(v(3.5));
            testCase.verifyFalse(v(NaN));
            testCase.verifyFalse(v([]));
            testCase.verifyFalse(v(Inf));
        end

        function layoutRule(testCase)
            v = iRule('dimNames');
            testCase.verifyTrue(v({'Y', 'X'}));
            testCase.verifyTrue(v({'Y', 'X', 'T', 'E'}));
            testCase.verifyTrue(v(["Y", "X", "F"]));
            testCase.verifyFalse(v({'X', 'Y', 'T'}));
            testCase.verifyFalse(v({'Y', 'X', 'E', 'T'}));
            testCase.verifyFalse(v({'Y', 'X', 'T', 'T'}));
            testCase.verifyFalse(v({'Y', 'X', 'Z'}));
            testCase.verifyFalse(v({}));
            testCase.verifyFalse(v('YXT'));
            testCase.verifyEqual(PipelineManager.isValidDatLayout({'Y', 'X', 'T'}), true);
        end
    end
end

function rule = iRule(field)
f = sourceInfoFields();
rule = f(strcmp({f.field}, field)).isValid;
end
