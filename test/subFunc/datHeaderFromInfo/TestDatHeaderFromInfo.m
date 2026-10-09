classdef TestDatHeaderFromInfo < matlab.unittest.TestCase
    %TESTDATHEADERFROMINFO Output header descriptions from a source Info.
    %
    %   .dat header Phase 4c-2a: datHeaderFromInfo copies class, axes,
    %   sizes, frame rate, and exposure from the source Info, applies
    %   overrides, and sanitizes the channel name as saveData does.

    methods (Test)
        function fieldsComeFromInfo(testCase)
            info = iInfo();
            hdr = datHeaderFromInfo(info, 'GSR');
            testCase.verifyEqual(hdr.dataClass, 'uint16');
            testCase.verifyEqual(hdr.dimNames, {'Y', 'X', 'T'});
            testCase.verifyEqual(hdr.dimSizes, [3 4 5]);
            testCase.verifyEqual(hdr.frameRateHz, 20);
            testCase.verifyEqual(hdr.exposureMsec, 2.5);
            testCase.verifyEqual(hdr.channelName, 'GSR');
        end

        function missingRateAndExposureAreNaN(testCase)
            info = rmfield(iInfo(), {'frameRateHz', 'exposureMsec'});
            hdr = datHeaderFromInfo(info, 'x', 'dimNames', {'Y', 'X', 'E'});
            testCase.verifyTrue(isnan(hdr.frameRateHz));
            testCase.verifyTrue(isnan(hdr.exposureMsec));
        end

        function overridesReplaceInfoValues(testCase)
            hdr = datHeaderFromInfo(iInfo(), 'x', 'dataClass', 'single', ...
                'dimNames', {'Y', 'X', 'T', 'E'}, 'dimSizes', [3 4 5 2], ...
                'frameRateHz', 7, 'exposureMsec', NaN);
            testCase.verifyEqual(hdr.dataClass, 'single');
            testCase.verifyEqual(hdr.dimNames, {'Y', 'X', 'T', 'E'});
            testCase.verifyEqual(hdr.dimSizes, [3 4 5 2]);
            testCase.verifyEqual(hdr.frameRateHz, 7);
            testCase.verifyTrue(isnan(hdr.exposureMsec));
        end

        function channelNameIsSanitized(testCase)
            [~, codes] = datHeaderSchema(1);
            maxChars = codes.constants.channelNameMaxChars;
            name = ['caf' char(233) '_' repmat('z', 1, 60)];
            hdr = datHeaderFromInfo(iInfo(), name);
            expected = strrep(name, char(233), '_');
            testCase.verifyEqual(hdr.channelName, expected(1:maxChars));
            testCase.verifyEqual(datHeaderFromInfo(iInfo(), "str").channelName, 'str');
        end

        function missingLayoutFieldErrors(testCase)
            testCase.verifyError(@() datHeaderFromInfo(rmfield(iInfo(), 'dimSizes'), 'x'), ...
                'Umitoolbox:datHeaderFromInfo:invalidInfo');
        end

        function roundTripThroughCreate(testCase)
            folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = fullfile(folder, 'out.dat');
            hdr = datHeaderFromInfo(iInfo(), 'out', 'dataClass', 'single');
            h = spatialSlabIO('create', f, hdr);
            spatialSlabIO('write', h, 1:4, ones(3, 4, 5, 'single'));
            spatialSlabIO('finalize', h);

            got = readDatHeader(f);
            testCase.verifyEqual(got.dataClass, 'single');
            testCase.verifyEqual(got.dimSizes, [3 4 5]);
            testCase.verifyEqual(got.frameRateHz, 20);
            testCase.verifyEqual(got.exposureMsec, 2.5);
            testCase.verifyEqual(got.channelName, 'out');
            testCase.verifyTrue(got.writeComplete);
        end
    end
end

function info = iInfo()
info = struct('filePath', 'C:\none\in.dat', 'format', 'header', 'dataOffset', 512, ...
    'dataClass', 'uint16', 'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [3 4 5], ...
    'frameRateHz', 20, 'exposureMsec', 2.5, 'channelName', 'in', 'writeComplete', true);
end
