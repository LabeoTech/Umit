classdef TestAssertDatLayout < matlab.unittest.TestCase
    %TESTASSERTDATLAYOUT Consumers refuse unsupported .dat layouts (Phase 6b-1).

    methods (Test)
        function acceptedLayoutsPass(testCase)
            Info = struct('dimNames', {{'Y', 'X', 'T'}}, 'filePath', 'a.dat');
            testCase.verifyWarningFree(@() assertDatLayout(Info, {{'Y', 'X', 'T'}}, 'fx'));
            yx = struct('dimNames', {{'Y', 'X'}}, 'filePath', 'b.dat');
            testCase.verifyWarningFree(@() assertDatLayout(yx, {{'Y', 'X', 'T'}, {'Y', 'X'}}, 'fx'));
            testCase.verifyWarningFree(@() assertDatLayout(yx, {["Y", "X"]}, "fx"));
        end

        function refusedLayoutRaisesNamedError(testCase)
            Info = struct('dimNames', {{'Y', 'X', 'T', 'E'}}, 'filePath', 'C:\data\ev.dat');
            try
                assertDatLayout(Info, {{'Y', 'X', 'T'}, {'Y', 'X'}}, 'myFunction');
                testCase.verifyFail('expected an error');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:myFunction:unsupportedLayout');
                testCase.verifySubstring(ME.message, 'C:\data\ev.dat');
                testCase.verifySubstring(ME.message, '{Y,X,T,E}');
                testCase.verifySubstring(ME.message, '{Y,X,T} or {Y,X}');
                testCase.verifySubstring(ME.message, 'Event-split .dat input is not supported');
            end
        end

        function nonEventLayoutHasNoEventWording(testCase)
            Info = struct('dimNames', {{'Y', 'X', 'F'}}, 'filePath', 'f.dat');
            try
                assertDatLayout(Info, {{'Y', 'X', 'T'}}, 'fx');
                testCase.verifyFail('expected an error');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:fx:unsupportedLayout');
                testCase.verifyFalse(contains(ME.message, 'Event-split'));
            end
        end

        function invalidInputsAreRejected(testCase)
            Info = struct('dimNames', {{'Y', 'X', 'T'}});
            testCase.verifyError(@() assertDatLayout(struct('a', 1), {{'Y', 'X', 'T'}}, 'fx'), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
            testCase.verifyError(@() assertDatLayout(Info, {}, 'fx'), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
            testCase.verifyError(@() assertDatLayout(Info, {{'Y', 'X', 'T'}}, ''), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
        end
    end
end
