classdef TestSaveDataHeader < matlab.unittest.TestCase
    %TESTSAVEDATAHEADER saveData writes headered .dat files.
    %
    %   .dat header Phase 4c-1, updated in Phase 6a: covers the header of
    %   new files, the required explicit axes (DimNames) and every
    %   supported layout, the frame-rate rule (FrameRateHz > Info, no
    %   AcqInfos.mat fallback, NaN without a T axis), the exposure and
    %   channelName rules, overwriting, and last-axis append with its
    %   checks.

    properties (Constant)
        YXT = {'Y', 'X', 'T'}
    end

    properties
        Folder
    end

    properties (TestParameter)
        layoutCase = struct( ...
            'yx', struct('names', {{'Y', 'X'}}, 'size', [4 3], 'rate', NaN), ...
            'yxt_t1', struct('names', {{'Y', 'X', 'T'}}, 'size', [4 3], 'rate', 10), ...
            'yxt', struct('names', {{'Y', 'X', 'T'}}, 'size', [4 3 5], 'rate', 10), ...
            'yxte', struct('names', {{'Y', 'X', 'T', 'E'}}, 'size', [4 3 5 2], 'rate', 10), ...
            'yxte_e1', struct('names', {{'Y', 'X', 'T', 'E'}}, 'size', [4 3 5], 'rate', 10), ...
            'yxe', struct('names', {{'Y', 'X', 'E'}}, 'size', [4 3 2], 'rate', NaN), ...
            'yxf', struct('names', {{'Y', 'X', 'F'}}, 'size', [4 3 6], 'rate', NaN))
        badNames = struct( ...
            'wrongOrder', {{'Y', 'X', 'E', 'T'}}, ...
            'duplicate', {{'Y', 'X', 'T', 'T'}}, ...
            'unknownAxis', {{'Y', 'X', 'Z'}}, ...
            'missingY', {{'X', 'T'}}, ...
            'xFirst', {{'X', 'Y', 'T'}}, ...
            'tooFewForArray', {{'Y', 'X'}})
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            testCase.Folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        function newFileHasHeader(testCase)
            data = iData(6, 5, 7);
            outFile = saveData(fullfile(testCase.Folder, 'green.whatever'), data, ...
                'DimNames', testCase.YXT, 'FrameRateHz', 20);

            testCase.verifyEqual(outFile, 'green.dat');
            f = fullfile(testCase.Folder, 'green.dat');
            testCase.verifyTrue(isDatWithHeader(f));
            hdr = readDatHeader(f);
            testCase.verifyEqual(hdr.dataClass, 'single');
            testCase.verifyEqual(hdr.dimNames, {'Y', 'X', 'T'});
            testCase.verifyEqual(hdr.dimSizes, [6 5 7]);
            testCase.verifyEqual(hdr.frameRateHz, 20);
            testCase.verifyTrue(isnan(hdr.exposureMsec));
            testCase.verifyEqual(hdr.channelName, 'green');
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(loadData(f), data);
            testCase.verifyFalse(isfile(fullfile(testCase.Folder, 'AcqInfos.mat')));
        end

        % ------------------------------------------------ layouts (Phase 6a)
        function everyLayoutRoundTrips(testCase, layoutCase)
            data = reshape(single(1:prod(layoutCase.size)), [layoutCase.size 1]);
            f = fullfile(testCase.Folder, 'layout.dat');
            saveData(f, data, 'DimNames', layoutCase.names, 'FrameRateHz', 10, ...
                'Info', struct('exposureMsec', 2));

            hdr = readDatHeader(f);
            nAxes = numel(layoutCase.names);
            testCase.verifyEqual(hdr.dimNames, layoutCase.names);
            testCase.verifyEqual(hdr.dimSizes, size(data, 1:nAxes));
            if isnan(layoutCase.rate)
                testCase.verifyTrue(isnan(hdr.frameRateHz), 'no T axis: NaN frame rate');
            else
                testCase.verifyEqual(hdr.frameRateHz, layoutCase.rate);
            end
            testCase.verifyEqual(hdr.exposureMsec, 2);
            testCase.verifyTrue(hdr.writeComplete);
            Info = loadMetaData(f);
            testCase.verifyEqual(Info.dimNames, layoutCase.names);
            testCase.verifyEqual(loadData(f), data);
        end

        function singleFrameAsYXTHasOneFrame(testCase)
            % DFR-20260929-003: a 2-D array is accepted as Y-X-T with T = 1.
            data = iData(4, 3, 1);
            f = fullfile(testCase.Folder, 'one.dat');
            saveData(f, data, 'DimNames', testCase.YXT, 'FrameRateHz', 10);
            testCase.verifyEqual(readDatHeader(f).dimSizes, [4 3 1]);
            testCase.verifyEqual(datAxisSize(loadMetaData(f), 'T'), 1);
            testCase.verifyEqual(loadData(f), data);
        end

        function missingDimNamesErrors(testCase)
            f = fullfile(testCase.Folder, 'nodims');
            testCase.verifyError(@() saveData(f, iData(4, 3, 5), 'FrameRateHz', 10), ...
                'Umitoolbox:saveData:missingDimNames');
            testCase.verifyFalse(isfile([f '.dat']));
        end

        function invalidDimNamesError(testCase, badNames)
            f = fullfile(testCase.Folder, 'bad.dat');
            saveData(f, iData(4, 3, 5), 'DimNames', testCase.YXT, 'FrameRateHz', 10);
            before = iBytes(f);
            testCase.verifyError(@() saveData(f, iData(4, 3, 5), 'DimNames', badNames, ...
                'FrameRateHz', 10), 'Umitoolbox:saveData:invalidDimNames');
            testCase.verifyEqual(iBytes(f), before, 'the existing file must be unchanged');
        end

        function umtStructNeedsNoDimNames(testCase)
            umt = genUMTStruct(iData(4, 3, 5), 'kind', 'image', 'dimNames', {'Y', 'X', 'T'});
            outFile = saveData(fullfile(testCase.Folder, 'derived'), umt);
            testCase.verifyEqual(outFile, 'derived.umt');
            testCase.verifyTrue(isfile(fullfile(testCase.Folder, 'derived.umt')));
        end

        % ------------------------------------------------ frame rate
        function frameRatePrecedence(testCase)
            data = iData(4, 3, 5);
            info = struct('frameRateHz', 25, 'exposureMsec', 4);

            saveData(fullfile(testCase.Folder, 'a'), data, 'DimNames', testCase.YXT, ...
                'FrameRateHz', 40, 'Info', info);
            saveData(fullfile(testCase.Folder, 'b'), data, 'DimNames', testCase.YXT, 'Info', info);

            testCase.verifyEqual(iRate(testCase.Folder, 'a'), 40);
            testCase.verifyEqual(iRate(testCase.Folder, 'b'), 25);
        end

        function noAcqInfosFrameRateFallback(testCase)
            % AcqInfos.mat describes the raw acquisition (temporal binning
            % changes the imported rate): it is never used (Phase 6a).
            iWriteAcqInfos(testCase.Folder, 10);
            data = iData(4, 3, 5);
            for k = 1:2
                if k == 1
                    args = {};
                else
                    args = {'Info', struct('frameRateHz', NaN, 'exposureMsec', 1)};
                end
                f = fullfile(testCase.Folder, sprintf('norate%d', k));
                testCase.verifyError(@() saveData(f, data, 'DimNames', testCase.YXT, args{:}), ...
                    'Umitoolbox:saveData:missingFrameRate');
                testCase.verifyFalse(isfile([f '.dat']));
            end
        end

        function frameRateIsIgnoredWithoutTAxis(testCase)
            f = fullfile(testCase.Folder, 'map.dat');
            saveData(f, iData(4, 3, 2), 'DimNames', {'Y', 'X', 'E'}, 'FrameRateHz', 40, ...
                'Info', struct('frameRateHz', 25));
            testCase.verifyTrue(isnan(readDatHeader(f).frameRateHz));
            f2 = fullfile(testCase.Folder, 'frame.dat');
            saveData(f2, iData(4, 3, 1), 'DimNames', {'Y', 'X'});
            testCase.verifyTrue(isnan(readDatHeader(f2).frameRateHz));
        end

        function acqInfoStreamArgumentIsRejected(testCase)
            % The positional AcqInfoStream argument was removed (Phase 6a);
            % saveData never writes AcqInfos.mat.
            acq = struct('Height', 4, 'Width', 3, 'FrameRateHz', 12.5);
            f = fullfile(testCase.Folder, 'x');
            testCase.verifyError(@() saveData(f, iData(4, 3, 5), acq), ?MException);
            testCase.verifyFalse(isfile(fullfile(testCase.Folder, 'AcqInfos.mat')));
            testCase.verifyFalse(isfile([f '.dat']));
        end

        % ------------------------------------------------ exposure and name
        function exposureIsInheritedFromInfo(testCase)
            saveData(fullfile(testCase.Folder, 'e'), iData(4, 3, 5), 'DimNames', testCase.YXT, ...
                'Info', struct('frameRateHz', 10, 'exposureMsec', 7.5));
            hdr = readDatHeader(fullfile(testCase.Folder, 'e.dat'));
            testCase.verifyEqual(hdr.exposureMsec, 7.5);
        end

        function longOrNonAsciiNameIsSanitized(testCase)
            base = ['PMTMP_20260929_' repmat('x', 1, 40) char(233)];
            saveData(fullfile(testCase.Folder, base), iData(4, 3, 5), 'DimNames', testCase.YXT, ...
                'FrameRateHz', 10);
            hdr = testCase.verifyWarningFree(@() readDatHeader(fullfile(testCase.Folder, [base '.dat'])));
            [~, codes] = datHeaderSchema(1);
            maxChars = codes.constants.channelNameMaxChars;
            testCase.verifyEqual(hdr.channelName, base(1:maxChars));
        end

        function channelNameOptionOverridesBaseName(testCase)
            saveData(fullfile(testCase.Folder, 'PMTMP_x'), iData(4, 3, 5), 'DimNames', testCase.YXT, ...
                'FrameRateHz', 10, 'ChannelName', 'GSR');
            testCase.verifyEqual(readDatHeader(fullfile(testCase.Folder, 'PMTMP_x.dat')).channelName, 'GSR');
        end

        function noAcqInfosSizeOrTimelineChecks(testCase)
            iWriteAcqInfos(testCase.Folder, 10);   % Height 4, Width 3, Length 5
            data = iData(9, 8, 13);
            saveData(fullfile(testCase.Folder, 'other'), data, 'DimNames', testCase.YXT, ...
                'FrameRateHz', 10);
            testCase.verifyEqual(loadData(fullfile(testCase.Folder, 'other.dat')), data);
        end

        function overwritesHeaderlessAndHeaderedFiles(testCase)
            f = fullfile(testCase.Folder, 'ow.dat');
            fid = fopen(f, 'w');
            fwrite(fid, zeros(1, 100, 'single'), 'single');
            fclose(fid);
            data = iData(4, 3, 5);
            saveData(f, data, 'DimNames', testCase.YXT, 'FrameRateHz', 10);
            testCase.verifyEqual(loadData(f), data);

            data2 = iData(2, 2, 3) + 100;
            saveData(f, data2, 'DimNames', testCase.YXT, 'FrameRateHz', 30);
            testCase.verifyEqual(loadData(f), data2);
            testCase.verifyEqual(readDatHeader(f).frameRateHz, 30);
        end

        % ------------------------------------------------ append
        function appendCreatesMissingFile(testCase)
            f = fullfile(testCase.Folder, 'app.dat');
            data = iData(4, 3, 5);
            saveData(f, data, 'DimNames', testCase.YXT, 'FrameRateHz', 10, 'Append', true);
            testCase.verifyEqual(loadData(f), data);
            testCase.verifyTrue(readDatHeader(f).writeComplete);
        end

        function appendGrowsLastAxis(testCase)
            f = fullfile(testCase.Folder, 'grow.dat');
            a = iData(4, 3, 5);
            b = iData(4, 3, 2) + 50;
            c = iData(4, 3, 1) + 90;
            info = struct('frameRateHz', 10, 'exposureMsec', 3);
            saveData(f, a, 'DimNames', testCase.YXT, 'Info', info);
            saveData(f, b, 'DimNames', testCase.YXT, 'Info', info, 'Append', true);
            saveData(f, cat(3, c, c), 'DimNames', testCase.YXT, 'Info', info, 'Append', true);
            saveData(f, c, 'DimNames', testCase.YXT, 'Info', info, 'Append', true);   % 2-D block, T = 1

            hdr = readDatHeader(f);
            testCase.verifyEqual(hdr.dimSizes, [4 3 10]);
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(hdr.exposureMsec, 3);
            testCase.verifyEqual(loadData(f), cat(3, a, b, c, c, c));
        end

        function appendGrowsEventAxis(testCase)
            names = {'Y', 'X', 'T', 'E'};
            f = fullfile(testCase.Folder, 'events.dat');
            a = reshape(single(1:4*3*5*2), 4, 3, 5, 2);
            b = reshape(single(1:4*3*5), 4, 3, 5) + 500;   % one more event (E = 1)
            saveData(f, a, 'DimNames', names, 'FrameRateHz', 10);
            saveData(f, b, 'DimNames', names, 'FrameRateHz', 10, 'Append', true);

            hdr = readDatHeader(f);
            testCase.verifyEqual(hdr.dimNames, names);
            testCase.verifyEqual(hdr.dimSizes, [4 3 5 3]);
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(loadData(f), cat(4, a, b));
        end

        function appendMismatchesAreRefused(testCase)
            f = fullfile(testCase.Folder, 'mis.dat');
            saveData(f, iData(4, 3, 5), 'DimNames', testCase.YXT, 'FrameRateHz', 10);
            before = iBytes(f);

            id = 'Umitoolbox:saveData:appendMismatch';
            testCase.verifyError(@() saveData(f, iData(5, 3, 2), 'DimNames', testCase.YXT, ...
                'FrameRateHz', 10, 'Append', true), id);
            testCase.verifyError(@() saveData(f, iData(4, 2, 2), 'DimNames', testCase.YXT, ...
                'FrameRateHz', 10, 'Append', true), id);
            testCase.verifyError(@() saveData(f, iData(4, 3, 2), 'DimNames', testCase.YXT, ...
                'FrameRateHz', 20, 'Append', true), id);
            testCase.verifyError(@() saveData(f, iData(4, 3, 2), 'DimNames', {'Y', 'X', 'E'}, ...
                'Append', true), id);
            testCase.verifyEqual(iBytes(f), before);
        end

        function appendOfDifferentFrameSizeIsRefused(testCase)
            % A 5 x 5 frame cannot be appended to a 100 x 100 x 10 file.
            f = fullfile(testCase.Folder, 'big.dat');
            saveData(f, zeros(100, 100, 10, 'single'), 'DimNames', testCase.YXT, 'FrameRateHz', 10);
            before = iBytes(f);
            try
                saveData(f, zeros(5, 5, 'single'), 'DimNames', testCase.YXT, ...
                    'FrameRateHz', 10, 'Append', true);
                testCase.verifyFail('the append should have been refused');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:saveData:appendMismatch');
                testCase.verifySubstring(ME.message, 'axis Y is 100 in the file but 5 in the data');
            end
            testCase.verifyEqual(iBytes(f), before);
        end

        function appendOfDifferentTInEventFileIsRefused(testCase)
            names = {'Y', 'X', 'T', 'E'};
            f = fullfile(testCase.Folder, 'ev.dat');
            saveData(f, zeros(4, 3, 5, 2, 'single'), 'DimNames', names, 'FrameRateHz', 10);
            before = iBytes(f);
            try
                saveData(f, zeros(4, 3, 6, 'single'), 'DimNames', names, ...
                    'FrameRateHz', 10, 'Append', true);
                testCase.verifyFail('the append should have been refused');
            catch ME
                testCase.verifyEqual(ME.identifier, 'Umitoolbox:saveData:appendMismatch');
                testCase.verifySubstring(ME.message, 'axis T is 5 in the file but 6 in the data');
            end
            testCase.verifyEqual(iBytes(f), before);
        end

        function appendToSingleFrameFileIsRefused(testCase)
            f = fullfile(testCase.Folder, 'frame.dat');
            saveData(f, iData(4, 3, 1), 'DimNames', {'Y', 'X'});
            before = iBytes(f);
            testCase.verifyError(@() saveData(f, iData(4, 3, 1), 'DimNames', {'Y', 'X'}, ...
                'Append', true), 'Umitoolbox:saveData:appendMismatch');
            testCase.verifyEqual(iBytes(f), before);
        end

        function appendToUnfinishedFileIsRefused(testCase)
            f = fullfile(testCase.Folder, 'unfinished.dat');
            hdr = struct('dataClass', 'single', 'frameRateHz', 10, 'exposureMsec', NaN, ...
                'channelName', 'unfinished', 'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [4 3 5]);
            h = spatialSlabIO('create', f, hdr);
            spatialSlabIO('close', h);
            before = iBytes(f);

            testCase.verifyError(@() saveData(f, iData(4, 3, 2), 'DimNames', testCase.YXT, ...
                'FrameRateHz', 10, 'Append', true), 'Umitoolbox:saveData:appendMismatch');
            testCase.verifyEqual(iBytes(f), before);
        end

        function appendToHeaderlessIsRefused(testCase)
            f = fullfile(testCase.Folder, 'legacy.dat');
            fid = fopen(f, 'w');
            fwrite(fid, iData(4, 3, 5), 'single');
            fclose(fid);
            before = iBytes(f);

            testCase.verifyError(@() saveData(f, iData(4, 3, 2), 'DimNames', testCase.YXT, ...
                'FrameRateHz', 10, 'Append', true), 'Umitoolbox:saveData:appendToHeaderless');
            testCase.verifyEqual(iBytes(f), before);
        end
    end
end

function data = iData(ny, nx, nt)
data = reshape(single(1:ny * nx * nt), ny, nx, nt);
end

function iWriteAcqInfos(folder, rate)
AcqInfoStream = struct('Height', 4, 'Width', 3, 'Length', 5, 'FrameRateHz', rate);
save(fullfile(folder, 'AcqInfos.mat'), 'AcqInfoStream');
end

function r = iRate(folder, base)
r = readDatHeader(fullfile(folder, [base '.dat'])).frameRateHz;
end

function b = iBytes(f)
fid = fopen(f, 'r');
b = fread(fid, inf, '*uint8');
fclose(fid);
end
