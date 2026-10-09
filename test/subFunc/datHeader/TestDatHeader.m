classdef TestDatHeader < matlab.unittest.TestCase
    %TESTDATHEADER Version-1 .dat header schema, encoder, decoder, reader, validator.
    %
    %   Protects the on-disk header format defined in
    %   docs/dev/dat-header-spec.md. Byte offsets, codes, and fixed values
    %   are hard-coded here from the spec on purpose, so datHeaderSchema is
    %   checked against the spec rather than against itself.
    %
    %   Reference files live in fixtures/v1 and are produced once by
    %   fixtures/makeDatFixturesV1.m; tests only read them.

    properties (Constant)
        Magic = uint8([137 85 77 68 13 10])
    end

    properties
        FixtureFolder
    end

    properties (TestParameter)
        fixtureCase = struct( ...
            'yxt_single', struct('file', 'yxt_single.dat', 'dataClass', 'single', ...
                'frameRateHz', 10, 'exposureMsec', 5, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [6 5 4], ...
                'writeComplete', true, 'values', single(1:120)), ...
            'yx_uint16', struct('file', 'yx_uint16.dat', 'dataClass', 'uint16', ...
                'frameRateHz', NaN, 'exposureMsec', NaN, 'channelName', 'amplitude', ...
                'dimNames', {{'Y', 'X'}}, 'dimSizes', [6 5], ...
                'writeComplete', true, 'values', uint16(1:30)), ...
            'yxte_single', struct('file', 'yxte_single.dat', 'dataClass', 'single', ...
                'frameRateHz', 10, 'exposureMsec', 5, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'T', 'E'}}, 'dimSizes', [4 3 5 2], ...
                'writeComplete', true, 'values', single(1:120)), ...
            'yxe_single', struct('file', 'yxe_single.dat', 'dataClass', 'single', ...
                'frameRateHz', NaN, 'exposureMsec', NaN, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'E'}}, 'dimSizes', [4 3 2], ...
                'writeComplete', true, 'values', single(1:24)), ...
            'yxt_incomplete', struct('file', 'yxt_incomplete.dat', 'dataClass', 'single', ...
                'frameRateHz', 10, 'exposureMsec', 5, 'channelName', 'green', ...
                'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [6 5 4], ...
                'writeComplete', false, 'values', single(1:120)))

        invalidCase = struct( ...
            'unknownDataClass',   struct('field', 'dataClass',   'value', 'logical'), ...
            'lengthMismatch',     struct('field', 'dimSizes',    'value', [6 5]), ...
            'unassignedAxis',     struct('field', 'dimNames',    'value', {{'Y', 'X', 'Z'}}), ...
            'missingY',           struct('field', 'dimNames',    'value', {{'X', 'T', 'E'}}), ...
            'missingX',           struct('field', 'dimNames',    'value', {{'Y', 'T', 'E'}}), ...
            'duplicateAxis',      struct('field', 'dimNames',    'value', {{'Y', 'X', 'X'}}), ...
            'outOfSlotOrder',     struct('field', 'dimNames',    'value', {{'X', 'Y', 'T'}}), ...
            'zeroSize',           struct('field', 'dimSizes',    'value', [6 0 4]), ...
            'nonIntegerSize',     struct('field', 'dimSizes',    'value', [6 5 4.5]), ...
            'sizeAboveUint32',    struct('field', 'dimSizes',    'value', [6 5 2^32]), ...
            'tAxisNaNRate',       struct('field', 'frameRateHz', 'value', NaN), ...
            'tAxisZeroRate',      struct('field', 'frameRateHz', 'value', 0), ...
            'tAxisInfiniteRate',  struct('field', 'frameRateHz', 'value', Inf))
    end

    methods (TestClassSetup)
        function locateFixtures(testCase)
            testCase.FixtureFolder = fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'v1');
            testCase.assertTrue(isfolder(testCase.FixtureFolder), ...
                'Reference files missing: run fixtures/makeDatFixturesV1 once.');
        end
    end

    methods (Test)
        % ------------------------------------------------------------ schema
        function schemaRejectsUnknownVersion(testCase)
            testCase.verifyError(@() datHeaderSchema(2), ...
                'Umitoolbox:datHeaderSchema:unknownHeaderVersion');
        end

        function layoutCoversHeaderContiguously(testCase)
            layout = datHeaderSchema(1);
            offsets = [layout.offset];
            nBytes = [layout.nBytes];
            testCase.verifyEqual(offsets(1), 0);
            testCase.verifyEqual(offsets(2:end), offsets(1:end-1) + nBytes(1:end-1));
            testCase.verifyEqual(sum(nBytes), 512);
        end

        function codeTablesMatchSpec(testCase)
            [~, codes] = datHeaderSchema(1);
            testCase.verifyEqual([codes.dataClass.code], 1:7);
            testCase.verifyEqual({codes.dataClass.className}, ...
                {'single', 'double', 'uint8', 'uint16', 'int16', 'uint32', 'int32'});
            testCase.verifyEqual([codes.dataClass.bytesPerValue], [4 8 1 2 2 4 4]);
            testCase.verifyEqual([codes.axis.slot], 0:7);
            testCase.verifyEqual({codes.axis(1:5).name}, {'Y', 'X', 'T', 'E', 'F'});
            testCase.verifyEqual([codes.axis.reserved], [false(1, 5) true(1, 3)]);
        end

        % ---------------------------------------------------------- encoding
        function encodeReturns512ByteRow(testCase)
            bytes = encodeDatHeader(testCase.yxtHeader());
            testCase.verifyClass(bytes, 'uint8');
            testCase.verifySize(bytes, [1 512]);
        end

        function byteOffsetsMatchSpec(testCase)
            b = encodeDatHeader(testCase.yxtHeader(true));
            % Offsets below are zero-based spec offsets + 1.
            testCase.verifyEqual(b(1:6), testCase.Magic, 'magic');
            testCase.verifyEqual(b(7:8), uint8([2 1]), 'byteOrderMark, little-endian');
            testCase.verifyEqual(b(9:10), uint8([1 0]), 'headerVersion');
            testCase.verifyEqual(b(11:12), uint8([1 0]), 'flags, write complete');
            testCase.verifyEqual(b(13:16), uint8([0 2 0 0]), 'headerLength 512');
            testCase.verifyEqual(b(17), uint8(1), 'dtype single');
            testCase.verifyEqual(b(18:21), uint8([0 0 32 65]), 'frameRateHz 10 (float32 LE)');
            testCase.verifyEqual(b(22:25), uint8([0 0 160 64]), 'exposureMsec 5 (float32 LE)');
            testCase.verifyEqual(b(26:30), uint8('green'), 'channelName');
            testCase.verifyEqual(b(31:57), zeros(1, 27, 'uint8'), 'channelName padding');
            testCase.verifyEqual(b(58:480), zeros(1, 423, 'uint8'), 'channel reserved');
            testCase.verifyEqual(b(481:484), uint8([6 0 0 0]), 'slot 0 Y');
            testCase.verifyEqual(b(485:488), uint8([5 0 0 0]), 'slot 1 X');
            testCase.verifyEqual(b(489:492), uint8([4 0 0 0]), 'slot 2 T');
            testCase.verifyEqual(b(493:512), zeros(1, 20, 'uint8'), 'slots 3-7');
        end

        function roundTripPreservesDescription(testCase, fixtureCase)
            decoded = decodeDatHeader(encodeDatHeader(testCase.headerFromCase(fixtureCase)));
            testCase.verifyDescriptionMatches(decoded, fixtureCase);
        end

        function absentMiddleAxisEncodesZeroSlot(testCase)
            h = testCase.yxtHeader();
            h.dimNames = {'Y', 'X', 'E'};
            h.dimSizes = [4 3 2];
            h.frameRateHz = NaN;
            b = encodeDatHeader(h);
            testCase.verifyEqual(b(489:492), zeros(1, 4, 'uint8'), 'T slot is 0');
            testCase.verifyEqual(b(493:496), uint8([2 0 0 0]), 'E stays in slot 3');
            decoded = decodeDatHeader(b);
            testCase.verifyEqual(decoded.dimNames, {'Y', 'X', 'E'});
            testCase.verifyEqual(decoded.dimSizes, [4 3 2]);
        end

        function writeCompleteOnlyChangesFlagBit(testCase)
            bytesOff = encodeDatHeader(testCase.yxtHeader(false));
            bytesOn = encodeDatHeader(testCase.yxtHeader(true));
            differing = find(bytesOff ~= bytesOn);
            testCase.verifyEqual(differing, 11, 'only the low byte of flags differs');
            testCase.verifyEqual(bitxor(bytesOff(11), bytesOn(11)), uint8(1));
            testCase.verifyFalse(decodeDatHeader(bytesOff).writeComplete);
            testCase.verifyTrue(decodeDatHeader(bytesOn).writeComplete);

            h = rmfield(testCase.yxtHeader(), 'writeComplete');
            testCase.verifyEqual(encodeDatHeader(h), bytesOff, 'omitted field encodes false');
        end

        function bigEndianHeaderDecodes(testCase)
            b = zeros(1, 512, 'uint8');
            b(1:6) = testCase.Magic;
            b(7:8) = [1 2];                    % byte-order mark, big-endian
            b(9:10) = [0 1];                   % headerVersion 1
            b(11:12) = [0 1];                  % flags: write complete
            b(13:16) = [0 0 2 0];              % headerLength 512
            b(17) = 4;                         % dtype uint16
            b(18:21) = [65 32 0 0];            % frameRateHz 10
            b(22:25) = [64 160 0 0];           % exposureMsec 5
            b(26:30) = uint8('green');
            b(481:484) = [0 0 0 6];            % Y
            b(485:488) = [0 0 0 5];            % X
            b(489:492) = [0 0 1 44];           % T = 300
            hdr = decodeDatHeader(b);
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(hdr.dataClass, 'uint16');
            testCase.verifyEqual(hdr.frameRateHz, 10);
            testCase.verifyEqual(hdr.exposureMsec, 5);
            testCase.verifyEqual(hdr.channelName, 'green');
            testCase.verifyEqual(hdr.dimNames, {'Y', 'X', 'T'});
            testCase.verifyEqual(hdr.dimSizes, [6 5 300]);
            testCase.verifyEqual(hdr.expectedDataBytes, 6 * 5 * 300 * 2);
        end

        % ------------------------------------------------------ channel name
        function longChannelNameIsTruncatedWithWarning(testCase)
            h = testCase.yxtHeader();
            h.channelName = repmat('a', 1, 32);
            b = testCase.verifyWarning(@() encodeDatHeader(h), ...
                'Umitoolbox:encodeDatHeader:channelNameTruncated');
            testCase.verifyEqual(decodeDatHeader(b).channelName, repmat('a', 1, 31));
        end

        function channelNameOf31CharsIsStoredUnchanged(testCase)
            h = testCase.yxtHeader();
            h.channelName = repmat('b', 1, 31);
            b = testCase.verifyWarningFree(@() encodeDatHeader(h));
            testCase.verifyEqual(decodeDatHeader(b).channelName, h.channelName);
        end

        function nonAsciiChannelNameIsRejected(testCase)
            h = testCase.yxtHeader();
            h.channelName = ['gr' char(233) 'en'];
            testCase.verifyError(@() encodeDatHeader(h), ...
                'Umitoolbox:encodeDatHeader:invalidChannelName');
        end

        function emptyChannelNameIsValid(testCase)
            h = testCase.yxtHeader();
            h.channelName = '';
            b = encodeDatHeader(h);
            testCase.verifyEqual(b(26:57), zeros(1, 32, 'uint8'));
            testCase.verifyEqual(decodeDatHeader(b).channelName, '');
        end

        % ---------------------------------------------------------- validator
        function validatorRejectsInvalidDescription(testCase, invalidCase)
            h = testCase.yxtHeader();
            h.(invalidCase.field) = invalidCase.value;
            testCase.verifyError(@() validateDatHeader(h), ...
                'Umitoolbox:validateDatHeader:invalidHeader');
        end

        function validatorRejectsFrameRateWithoutTAxis(testCase)
            h = testCase.yxtHeader();
            h.dimNames = {'Y', 'X'};
            h.dimSizes = [6 5];
            testCase.verifyError(@() validateDatHeader(h), ...
                'Umitoolbox:validateDatHeader:invalidHeader');
        end

        function validatorAcceptsValidDescription(testCase)
            testCase.verifyWarningFree(@() validateDatHeader(testCase.yxtHeader()));
        end

        function fileBytesMismatchErrorsWhenComplete(testCase)
            h = testCase.yxtHeader(true);
            testCase.verifyError(@() validateDatHeader(h, 'FileBytes', 512 + 479), ...
                'Umitoolbox:validateDatHeader:fileSizeMismatch');
        end

        function fileBytesMismatchWarnsWhenIncomplete(testCase)
            h = testCase.yxtHeader(false);
            testCase.verifyWarning(@() validateDatHeader(h, 'FileBytes', 512 + 100), ...
                'Umitoolbox:validateDatHeader:incompleteFile');
        end

        function fileBytesMatchIsSilent(testCase)
            h = testCase.yxtHeader(true);
            testCase.verifyWarningFree(@() validateDatHeader(h, 'FileBytes', 512 + 480));
        end

        function nonzeroReservedBytesWarn(testCase)
            b = encodeDatHeader(testCase.yxtHeader());
            b(100) = 7;
            testCase.verifyWarning(@() decodeDatHeader(b), ...
                'Umitoolbox:validateDatHeader:reservedBytesNotZero');
        end

        % ------------------------------------------------------- corruption
        function corruptHeadersAreRejected(testCase)
            valid = encodeDatHeader(testCase.yxtHeader());
            corrupt = {
                'headerVersion 2',        @(b) setBytes(b, 9:10, [2 0])
                'reserved axis slot',     @(b) setBytes(b, 501, 1)
                'fewer than 512 bytes',   @(b) b(1:300)
                'headerLength 1024',      @(b) setBytes(b, 13:16, [0 4 0 0])
                'reserved flag bit set',  @(b) setBytes(b, 11, 2)
                'invalid byte-order mark', @(b) setBytes(b, 7:8, [0 0])
                'unterminated name',      @(b) setBytes(b, 26:57, uint8('a'))};
            for k = 1:size(corrupt, 1)
                bad = corrupt{k, 2}(valid);
                testCase.verifyError(@() decodeDatHeader(bad), ...
                    'Umitoolbox:decodeDatHeader:corruptHeader', corrupt{k, 1});
            end

            function b = setBytes(b, idx, v)
                b(idx) = v;
            end
        end

        function bytesWithoutMagicAreNotAHeader(testCase)
            testCase.verifyError(@() decodeDatHeader(typecast(single(1:128), 'uint8')), ...
                'Umitoolbox:decodeDatHeader:noHeader');
        end

        % ---------------------------------------------- detection and reading
        function detectionOnReferenceFiles(testCase, fixtureCase)
            testCase.verifyTrue(isDatWithHeader(fullfile(testCase.FixtureFolder, fixtureCase.file)));
        end

        function detectionRejectsLegacyAndShortFiles(testCase)
            testCase.verifyFalse(isDatWithHeader(testCase.legacyFile()));
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            emptyFile = fullfile(folder, 'empty.dat');
            fclose(fopen(emptyFile, 'w'));
            testCase.verifyFalse(isDatWithHeader(emptyFile));
            shortFile = fullfile(folder, 'short.dat');
            fid = fopen(shortFile, 'w');
            fwrite(fid, testCase.Magic(1:5), 'uint8');
            fclose(fid);
            testCase.verifyFalse(isDatWithHeader(shortFile));
        end

        function readDatHeaderOnReferenceFiles(testCase, fixtureCase)
            filePath = fullfile(testCase.FixtureFolder, fixtureCase.file);
            hdr = readDatHeader(filePath);
            testCase.verifyDescriptionMatches(hdr, fixtureCase);
            testCase.verifyEqual(hdr.headerVersion, 1);
            testCase.verifyEqual(hdr.dataOffset, 512);
            testCase.verifyEqual(hdr.expectedDataBytes, ...
                numel(fixtureCase.values) * hdr.bytesPerValue);

            info = dir(filePath);
            testCase.verifyWarningFree(@() validateDatHeader(hdr, 'FileBytes', info.bytes));

            fid = fopen(filePath, 'r', 'ieee-le');
            cleanupObj = onCleanup(@() fclose(fid));
            fseek(fid, hdr.dataOffset, 'bof');
            values = fread(fid, inf, ['*' hdr.dataClass]);
            clear cleanupObj
            testCase.verifyEqual(values(:).', fixtureCase.values, 'data values');
        end

        function readDatHeaderRejectsLegacyFile(testCase)
            testCase.verifyError(@() readDatHeader(testCase.legacyFile()), ...
                'Umitoolbox:readDatHeader:noHeader');
        end

        function readDatHeaderRejectsMissingFile(testCase)
            testCase.verifyError(@() readDatHeader(fullfile(tempdir, 'no_such_file.dat')), ...
                'Umitoolbox:readDatHeader:fileNotFound');
        end
    end

    methods (Access = private)
        function h = yxtHeader(~, writeComplete)
            if nargin < 2
                writeComplete = false;
            end
            h = struct('dataClass', 'single', 'frameRateHz', 10, 'exposureMsec', 5, ...
                'channelName', 'green', 'dimNames', {{'Y', 'X', 'T'}}, ...
                'dimSizes', [6 5 4], 'writeComplete', writeComplete);
        end

        function h = headerFromCase(~, c)
            h = struct('dataClass', c.dataClass, 'frameRateHz', c.frameRateHz, ...
                'exposureMsec', c.exposureMsec, 'channelName', c.channelName, ...
                'dimNames', {c.dimNames}, 'dimSizes', c.dimSizes, ...
                'writeComplete', c.writeComplete);
        end

        function verifyDescriptionMatches(testCase, hdr, c)
            testCase.verifyEqual(hdr.dataClass, c.dataClass);
            testCase.verifyEqual(hdr.frameRateHz, c.frameRateHz);
            testCase.verifyEqual(hdr.exposureMsec, c.exposureMsec);
            testCase.verifyEqual(hdr.channelName, c.channelName);
            testCase.verifyEqual(hdr.dimNames, c.dimNames);
            testCase.verifyEqual(hdr.dimSizes, c.dimSizes);
            testCase.verifyEqual(hdr.writeComplete, c.writeComplete);
        end

        function f = legacyFile(testCase)
            f = fullfile(testCase.FixtureFolder, 'legacy', 'green.dat');
        end
    end
end
