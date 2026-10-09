classdef TestSpatialSlabIO < matlab.unittest.TestCase
    %TESTSPATIALSLABIO Handle-based spatialSlabIO reads and the unchanged positional form.
    %
    %   Every handle-based read is compared with the matching part of the
    %   loadData array, on the header and legacy sidecar reference files in
    %   test/subFunc/datHeader/fixtures; the AcqInfos-bound reference file is
    %   rejected (.dat header Phase 5b). The positional
    %   form is characterized on a temporary headerless file.

    properties
        FixtureRoot
    end

    properties (TestParameter)
        fileCase = struct( ...
            'yxt_single', struct('rel', {{'v1', 'yxt_single.dat'}}, 'size', [6 5 4], 'class', 'single'), ...
            'yx_uint16', struct('rel', {{'v1', 'yx_uint16.dat'}}, 'size', [6 5], 'class', 'uint16'), ...
            'yxte_single', struct('rel', {{'v1', 'yxte_single.dat'}}, 'size', [4 3 5 2], 'class', 'single'), ...
            'yxe_single', struct('rel', {{'v1', 'yxe_single.dat'}}, 'size', [4 3 2], 'class', 'single'), ...
            'legacySidecar', struct('rel', {{'legacy_sidecar', 'green.dat'}}, 'size', [6 5 4], 'class', 'single'))
        writeCase = struct( ...
            'yxt_single', struct('dataClass', 'single', 'frameRateHz', 10, ...
                'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [6 5 4]), ...
            'yxte_single', struct('dataClass', 'single', 'frameRateHz', 10, ...
                'dimNames', {{'Y', 'X', 'T', 'E'}}, 'dimSizes', [4 3 5 2]), ...
            'yxe_single', struct('dataClass', 'single', 'frameRateHz', NaN, ...
                'dimNames', {{'Y', 'X', 'E'}}, 'dimSizes', [4 3 2]), ...
            'yx_uint16', struct('dataClass', 'uint16', 'frameRateHz', NaN, ...
                'dimNames', {{'Y', 'X'}}, 'dimSizes', [6 5]))
    end

    methods (TestClassSetup)
        function locateFixtures(testCase)
            testCase.FixtureRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                'datHeader', 'fixtures');
            testCase.assertTrue(isfile(fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat')), ...
                'Reference files missing: run makeDatFixturesV1 once.');
            testCase.assertTrue(isfile(fullfile(testCase.FixtureRoot, 'legacy_sidecar', 'green.dat')), ...
                'Legacy reference file missing: run makeLegacySidecarFixture once.');
        end
    end

    methods (Test)
        function openReportsLayout(testCase, fileCase)
            h = testCase.openCase(fileCase);
            sz = fileCase.size;
            testCase.verifyEqual(h.Ny, sz(1));
            testCase.verifyEqual(h.Nx, sz(2));
            if numel(sz) > 2
                testCase.verifyEqual(h.trailingSizes, sz(3:end));
            else
                testCase.verifyEqual(h.trailingSizes, 1);
            end
            testCase.verifyEqual(h.nFrames, prod(h.trailingSizes));
            testCase.verifyEqual(h.dataClass, fileCase.class);
        end

        function fullReadKeepsTrailingShape(testCase, fileCase)
            [h, ref] = testCase.openWithReference(fileCase);
            slab = spatialSlabIO('read', h, 1:h.Nx);
            testCase.verifyClass(slab, fileCase.class);
            testCase.verifyEqual(slab, ref);
        end

        function columnSelectionsMatchLoadData(testCase, fileCase)
            [h, ref] = testCase.openWithReference(fileCase);
            flat = reshape(ref, h.Ny, h.Nx, []);
            selections = {2:3, [1 3], [3 2 1], 2, 1:h.Nx};
            for k = 1:numel(selections)
                xIdx = selections{k};
                slab = spatialSlabIO('read', h, xIdx);
                expected = reshape(flat(:, xIdx, :), [h.Ny, numel(xIdx), h.trailingSizes]);
                testCase.verifyEqual(slab, expected, sprintf('columns [%s]', num2str(xIdx)));
            end
        end

        function frameSelectionsMatchLoadData(testCase, fileCase)
            [h, ref] = testCase.openWithReference(fileCase);
            flat = reshape(ref, h.Ny, h.Nx, []);
            n = h.nFrames;
            if n == 1
                selections = {1};
            else
                selections = {1, n, 1:n, n:-1:1, [n 1]};
            end
            if n >= 3
                selections = [selections, {2:n-1, [3 1 2]}];
            end
            for k = 1:numel(selections)
                frameIdx = selections{k};
                for xIdx = {1:h.Nx, [h.Nx 1]}
                    slab = spatialSlabIO('read', h, xIdx{1}, frameIdx);
                    testCase.verifySize(slab, [h.Ny, numel(xIdx{1}), numel(frameIdx)]);
                    testCase.verifyEqual(slab, flat(:, xIdx{1}, frameIdx), ...
                        sprintf('frames [%s], columns [%s]', num2str(frameIdx), num2str(xIdx{1})));
                end
            end
        end

        function trailingAxesShapes(testCase)
            h = testCase.openCase(testCase.fileCase.yxte_single);
            testCase.verifySize(spatialSlabIO('read', h, 1:2), [4, 2, 5, 2]);
            h2 = testCase.openCase(testCase.fileCase.yxe_single);
            testCase.verifySize(spatialSlabIO('read', h2, 1:2), [4, 2, 2]);
            h3 = testCase.openCase(testCase.fileCase.yx_uint16);
            testCase.verifySize(spatialSlabIO('read', h3, 1:2), [6, 2]);
        end

        function invalidInputIsRejected(testCase)
            h = testCase.openCase(testCase.fileCase.yxt_single);
            id = 'Umitoolbox:spatialSlabIO:invalidInput';
            testCase.verifyError(@() spatialSlabIO('peek', h, 1), id);
            testCase.verifyError(@() spatialSlabIO('read', h, 0), id);
            testCase.verifyError(@() spatialSlabIO('read', h, h.Nx + 1), id);
            testCase.verifyError(@() spatialSlabIO('read', h, [2 2]), id);
            testCase.verifyError(@() spatialSlabIO('read', h, 1, 0), id);
            testCase.verifyError(@() spatialSlabIO('read', h, 1, h.nFrames + 1), id);
            testCase.verifyError(@() spatialSlabIO('read', h, 1.5), id);
            testCase.verifyError(@() spatialSlabIO('read', struct('fid', 3), 1), id);
            testCase.verifyError(@() spatialSlabIO('close', 3), id);
        end

        function closeClosesTheFileAndIsIdempotent(testCase)
            f = fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat');
            h = spatialSlabIO('open', f);
            testCase.verifyNotEmpty(fopen(h.fid));
            spatialSlabIO('close', h);
            testCase.verifyEmpty(fopen(h.fid));
            testCase.verifyWarningFree(@() spatialSlabIO('close', h));
            testCase.verifyError(@() spatialSlabIO('read', h, 1), ...
                'Umitoolbox:spatialSlabIO:invalidInput');
        end

        function nonYXLayoutIsRefused(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = fullfile(folder, 'events.dat');
            fid = fopen(f, 'w');
            fwrite(fid, single(1:60), 'single');
            fclose(fid);
            dim_names = {'E', 'Y', 'X', 'T'};
            datSize = [2 3 5];
            datLength = 2;
            Freq = 10;
            Datatype = 'single';
            save(fullfile(folder, 'events.mat'), 'dim_names', 'datSize', 'datLength', 'Freq', 'Datatype');
            testCase.verifyError(@() spatialSlabIO('open', f), ...
                'Umitoolbox:spatialSlabIO:unsupportedLayout');
        end

        function acqInfosBoundFileIsRejected(testCase)
            f = fullfile(testCase.FixtureRoot, 'v1', 'legacy', 'green.dat');
            testCase.verifyError(@() spatialSlabIO('open', f), 'Umitoolbox:loadMetaData:acqInfosBoundUnsupported');
        end

        % ------------------------------------------------ Info reuse
        function openWithInfoMatchesPlainOpen(testCase, fileCase)
            f = fullfile(testCase.FixtureRoot, fileCase.rel{:});
            [~, Info] = evalc('loadMetaData(f)');
            hPlain = testCase.openCase(fileCase);
            hInfo = spatialSlabIO('open', f, 'Info', Info);
            testCase.addTeardown(@() spatialSlabIO('close', hInfo));

            testCase.verifyEqual(hInfo.dataOffset, hPlain.dataOffset);
            testCase.verifyEqual(hInfo.trailingSizes, hPlain.trailingSizes);
            testCase.verifyEqual(spatialSlabIO('read', hInfo, 1:hInfo.Nx), ...
                spatialSlabIO('read', hPlain, 1:hPlain.Nx));
            testCase.verifyEqual(spatialSlabIO('read', hInfo, [hInfo.Nx 1], hInfo.nFrames), ...
                spatialSlabIO('read', hPlain, [hPlain.Nx 1], hPlain.nFrames));
        end

        function openWithForeignOrIncompleteInfoIsRejected(testCase)
            fHeader = fullfile(testCase.FixtureRoot, 'v1', 'yxt_single.dat');
            fLegacy = fullfile(testCase.FixtureRoot, 'legacy_sidecar', 'green.dat');
            [~, headerInfo] = evalc('loadMetaData(fHeader)');
            id = 'Umitoolbox:spatialSlabIO:invalidInput';

            testCase.verifyError(@() spatialSlabIO('open', fLegacy, 'Info', headerInfo), id);
            testCase.verifyError(@() spatialSlabIO('open', fHeader, 'Info', rmfield(headerInfo, 'dataOffset')), id);
            testCase.verifyError(@() spatialSlabIO('open', fHeader, 'Info', 42), id);
            testCase.verifyError(@() spatialSlabIO('open', fHeader, 'Layout', headerInfo), id);
        end

        % ------------------------------------------------ write side
        function createWriteFinalizeRoundTrip(testCase, writeCase)
            f = fullfile(testCase.tempFolder(), 'out.dat');
            values = testCase.knownValues(writeCase);
            h = spatialSlabIO('create', f, testCase.headerFor(writeCase));
            Nx = h.Nx;
            % Whole file as several slabs: reversed, non-contiguous, contiguous.
            pieces = {Nx:-1:max(1, Nx - 1), 1:2:max(1, Nx - 2), []};
            written = unique([pieces{:}]);
            pieces{3} = setdiff(1:Nx, written);
            flat = reshape(values, h.Ny, Nx, []);
            for k = 1:numel(pieces)
                if isempty(pieces{k})
                    continue
                end
                slab = reshape(flat(:, pieces{k}, :), [h.Ny, numel(pieces{k}), h.trailingSizes]);
                spatialSlabIO('write', h, pieces{k}, slab);
            end
            spatialSlabIO('finalize', h);

            data = loadData(f);
            testCase.verifyClass(data, writeCase.dataClass);
            testCase.verifyEqual(data, values);
            Info = loadMetaData(f);
            testCase.verifyEqual(Info.format, 'header');
            testCase.verifyTrue(Info.writeComplete);
            testCase.verifyEqual(Info.dimNames, writeCase.dimNames);
            testCase.verifyEqual(Info.dimSizes, writeCase.dimSizes);
            testCase.verifyEqual(Info.frameRateHz, writeCase.frameRateHz);
            d = dir(f);
            testCase.verifyEqual(d.bytes, 512 + numel(values) * getByteSize(writeCase.dataClass));
        end

        function frameChunkWritesRoundTrip(testCase)
            c = testCase.writeCase.yxte_single;
            f = fullfile(testCase.tempFolder(), 'frames.dat');
            values = testCase.knownValues(c);
            h = spatialSlabIO('create', f, testCase.headerFor(c));
            flat = reshape(values, h.Ny, h.Nx, h.nFrames);
            % Frames crossing the event boundary, then scattered, then the rest.
            groups = {4:7, [10 1 9], []};
            groups{3} = setdiff(1:h.nFrames, [groups{1:2}]);
            for k = 1:numel(groups)
                spatialSlabIO('write', h, 1:h.Nx, flat(:, :, groups{k}), groups{k});
            end
            % Partial columns of one frame overwrite with the same values.
            spatialSlabIO('write', h, [3 1], flat(:, [3 1], 2), 2);
            spatialSlabIO('finalize', h);
            testCase.verifyEqual(loadData(f), values);
        end

        function fileIsIncompleteUntilFinalized(testCase)
            c = testCase.writeCase.yxt_single;
            f = fullfile(testCase.tempFolder(), 'incomplete.dat');
            h = spatialSlabIO('create', f, testCase.headerFor(c));
            testCase.verifyFalse(readDatHeader(f).writeComplete);
            testCase.verifyWarning(@() loadMetaData(f), 'Umitoolbox:loadMetaData:incompleteFile');
            d = dir(f);
            testCase.verifyEqual(d.bytes, 512 + prod(c.dimSizes) * 4, 'full size before any write');

            testCase.verifyEqual(spatialSlabIO('read', h, 1:h.Nx), zeros([h.Ny h.Nx h.trailingSizes], 'single'), ...
                'unwritten regions read back as zero');
            spatialSlabIO('write', h, 2, ones(h.Ny, 1, h.nFrames, 'single'));
            testCase.verifyEqual(spatialSlabIO('read', h, 2), ones(h.Ny, 1, h.nFrames, 'single'), ...
                'read on a created handle returns what was written');

            spatialSlabIO('finalize', h);
            testCase.verifyTrue(readDatHeader(f).writeComplete);
            testCase.verifyWarningFree(@() loadMetaData(f));
        end

        function closeWithoutFinalizeLeavesIncompleteFile(testCase)
            c = testCase.writeCase.yxt_single;
            f = fullfile(testCase.tempFolder(), 'aborted.dat');
            h = spatialSlabIO('create', f, testCase.headerFor(c));
            spatialSlabIO('close', h);
            testCase.verifyFalse(readDatHeader(f).writeComplete);
            testCase.verifyEmpty(fopen(h.fid));
        end

        function invalidWritesAreRejected(testCase)
            c = testCase.writeCase.yxt_single;
            folder = testCase.tempFolder();
            id = 'Umitoolbox:spatialSlabIO:invalidInput';

            bad = testCase.headerFor(c);
            bad.dimNames = {'X', 'Y', 'T'};
            testCase.verifyError(@() spatialSlabIO('create', fullfile(folder, 'bad.dat'), bad), ...
                'Umitoolbox:validateDatHeader:invalidHeader');
            testCase.verifyFalse(isfile(fullfile(folder, 'bad.dat')), 'an invalid header creates no file');

            h = spatialSlabIO('create', fullfile(folder, 'w.dat'), testCase.headerFor(c));
            testCase.addTeardown(@() spatialSlabIO('close', h));
            testCase.verifyError(@() spatialSlabIO('write', h, 1, zeros(h.Ny, 2, h.nFrames, 'single')), id);
            testCase.verifyError(@() spatialSlabIO('write', h, [1 1], zeros(h.Ny, 2, h.nFrames, 'single')), id);
            testCase.verifyError(@() spatialSlabIO('write', h, h.Nx + 1, zeros(h.Ny, 1, h.nFrames, 'single')), id);
            testCase.verifyError(@() spatialSlabIO('write', h, 1, zeros(h.Ny, 1, 1, 'single'), h.nFrames + 1), id);

            hRead = testCase.openCase(testCase.fileCase.yxt_single);
            testCase.verifyError(@() spatialSlabIO('write', hRead, 1, zeros(hRead.Ny, 1, hRead.nFrames, 'single')), id);
            testCase.verifyError(@() spatialSlabIO('finalize', hRead), id);

            h2 = spatialSlabIO('create', fullfile(folder, 'closed.dat'), testCase.headerFor(c));
            spatialSlabIO('finalize', h2);
            testCase.verifyError(@() spatialSlabIO('write', h2, 1, zeros(h2.Ny, 1, h2.nFrames, 'single')), id);
            testCase.verifyError(@() spatialSlabIO('finalize', h2), id);
        end

        % ------------------------------------------------ positional form
        function createLeavesLastwarnUnchanged(testCase)
            % DFR-20260929-001: 'create' describes the incomplete file with a
            % disabled warning, which must not replace the caller's lastwarn.
            hdr = struct('dataClass', 'single', 'frameRateHz', 10, ...
                'exposureMsec', NaN, 'channelName', 'lw', ...
                'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [3 4 5]);
            f = fullfile(testCase.tempFolder(), 'lastwarn.dat');

            lastwarn('caller message', 'Umitoolbox:TestSpatialSlabIO:caller');
            h = spatialSlabIO('create', f, hdr);
            spatialSlabIO('close', h);
            [msg, id] = lastwarn();
            testCase.verifyEqual(id, 'Umitoolbox:TestSpatialSlabIO:caller');
            testCase.verifyEqual(msg, 'caller message');

            bad = hdr;
            bad.dimSizes = [3 0 5];
            testCase.verifyError(@() spatialSlabIO('create', f, bad), ?MException);
            [msg, id] = lastwarn();
            testCase.verifyEqual(id, 'Umitoolbox:TestSpatialSlabIO:caller');
            testCase.verifyEqual(msg, 'caller message');
        end

        function growableAppendsRoundTrip(testCase)
            % .dat header Phase 4d: many appends of mixed sizes, including
            % single [Ny Nx] frames, read back like one array.
            f = fullfile(testCase.tempFolder(), 'grow.dat');
            h = spatialSlabIO('create', f, iGrowHeader('uint16'), 'Growable', true);
            blocks = {uint16(reshape(1:3*4*2, 3, 4, 2)), uint16(reshape(101:112, 3, 4)), ...
                uint16(reshape(201:3*4*5+200, 3, 4, 5)), uint16(ones(3, 4))};
            for k = 1:numel(blocks)
                h = spatialSlabIO('append', h, blocks{k});
            end
            testCase.verifyEqual(h.framesWritten, 9);
            spatialSlabIO('finalize', h);

            hdr = readDatHeader(f);
            testCase.verifyTrue(hdr.writeComplete);
            testCase.verifyEqual(hdr.dimSizes, [3 4 9]);
            testCase.verifyEqual(hdr.dataClass, 'uint16');
            testCase.verifyEqual(loadData(f), cat(3, blocks{:}));
            testCase.verifyEqual(dir(f).bytes, 512 + 3 * 4 * 9 * 2);
        end

        function growableFilesCanBeGrownInTurn(testCase)
            folder = testCase.tempFolder();
            fa = fullfile(folder, 'a.dat');
            fb = fullfile(folder, 'b.dat');
            ha = spatialSlabIO('create', fa, iGrowHeader('single'), 'Growable', true);
            hb = spatialSlabIO('create', fb, iGrowHeader('single'), 'Growable', true);
            a = cell(1, 3);
            b = cell(1, 3);
            for k = 1:3
                a{k} = rand(3, 4, k, 'single');
                b{k} = rand(3, 4, 4 - k, 'single');
                ha = spatialSlabIO('append', ha, a{k});
                hb = spatialSlabIO('append', hb, b{k});
            end
            spatialSlabIO('finalize', ha);
            spatialSlabIO('finalize', hb);
            testCase.verifyEqual(loadData(fa), cat(3, a{:}));
            testCase.verifyEqual(loadData(fb), cat(3, b{:}));
        end

        function growableFileIsIncompleteUntilFinalized(testCase)
            f = fullfile(testCase.tempFolder(), 'partial.dat');
            h = spatialSlabIO('create', f, iGrowHeader('single'), 'Growable', true);
            hdr = readDatHeader(f);
            testCase.verifyFalse(hdr.writeComplete);
            testCase.verifyEqual(hdr.dimSizes, [3 4 1]);
            testCase.verifyEqual(dir(f).bytes, 512);

            h = spatialSlabIO('append', h, ones(3, 4, 2, 'single'));
            hdr = readDatHeader(f);
            testCase.verifyFalse(hdr.writeComplete);
            testCase.verifyEqual(hdr.dimSizes, [3 4 2]);
            testCase.verifyWarning(@() loadMetaData(f), 'Umitoolbox:loadMetaData:incompleteFile');

            h = spatialSlabIO('append', h, 2 * ones(3, 4, 'single'));
            spatialSlabIO('close', h);
            hdr = readDatHeader(f);
            testCase.verifyFalse(hdr.writeComplete, 'close without finalize leaves it incomplete');
            testCase.verifyEqual(hdr.dimSizes, [3 4 3]);
        end

        function growableMisuseIsRejected(testCase)
            folder = testCase.tempFolder();
            h = spatialSlabIO('create', fullfile(folder, 'g.dat'), iGrowHeader('single'), 'Growable', true);
            cleanup = onCleanup(@() spatialSlabIO('close', h));
            testCase.verifyError(@() spatialSlabIO('finalize', h), 'Umitoolbox:spatialSlabIO:invalidInput');
            testCase.verifyError(@() spatialSlabIO('write', h, 1:4, ones(3, 4, 1, 'single')), ...
                'Umitoolbox:spatialSlabIO:invalidInput');
            testCase.verifyError(@() spatialSlabIO('read', h, 1), 'Umitoolbox:spatialSlabIO:invalidInput');
            testCase.verifyError(@() spatialSlabIO('append', h, ones(4, 4, 'single')), ...
                'Umitoolbox:spatialSlabIO:invalidInput');
            testCase.verifyError(@() spatialSlabIO('append', h, ones(3, 4, 2, 2, 'single')), ...
                'Umitoolbox:spatialSlabIO:invalidInput');

            normal = spatialSlabIO('create', fullfile(folder, 'n.dat'), ...
                setfield(iGrowHeader('single'), 'dimSizes', [3 4 2]));
            cleanupNormal = onCleanup(@() spatialSlabIO('close', normal));
            testCase.verifyError(@() spatialSlabIO('append', normal, ones(3, 4, 'single')), ...
                'Umitoolbox:spatialSlabIO:invalidInput');

            fourAxes = iGrowHeader('single');
            fourAxes.dimNames = {'Y', 'X', 'T', 'E'};
            fourAxes.dimSizes = [3 4 2 2];
            testCase.verifyError(@() spatialSlabIO('create', fullfile(folder, 'e.dat'), fourAxes, ...
                'Growable', true), 'Umitoolbox:spatialSlabIO:unsupportedLayout');
            testCase.verifyError(@() spatialSlabIO('create', fullfile(folder, 'o.dat'), ...
                iGrowHeader('single'), 'Growable', 'yes'), 'Umitoolbox:spatialSlabIO:invalidInput');
            clear cleanup cleanupNormal
        end

        function growableCreateLeavesLastwarnUnchanged(testCase)
            lastwarn('caller message', 'Umitoolbox:TestSpatialSlabIO:caller');
            h = spatialSlabIO('create', fullfile(testCase.tempFolder(), 'lw.dat'), ...
                iGrowHeader('single'), 'Growable', true);
            h = spatialSlabIO('append', h, ones(3, 4, 'single'));
            spatialSlabIO('finalize', h);
            [msg, id] = lastwarn();
            testCase.verifyEqual(id, 'Umitoolbox:TestSpatialSlabIO:caller');
            testCase.verifyEqual(msg, 'caller message');
        end

        function positionalReadAndWriteUnchanged(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = fullfile(folder, 'raw.dat');
            data = reshape(single(1:6*5*4), 6, 5, 4);
            fid = fopen(f, 'w');
            fwrite(fid, data, 'single');
            fclose(fid);

            fid = fopen(f, 'r+');
            cleanupObj = onCleanup(@() fclose(fid));
            slab = spatialSlabIO('read', fid, 6, 5, 4, 2:3, 'single');
            testCase.verifyClass(slab, 'single');
            testCase.verifyEqual(slab, data(:, 2:3, :));
            testCase.verifyEqual(spatialSlabIO('read', fid, 6, 5, 4, [1 4], 'single'), data(:, [1 4], :));

            newCols = -reshape(single(1:6*2*4), 6, 2, 4);
            spatialSlabIO('write', fid, 6, 5, 4, [2 5], 'single', newCols);
            expected = data;
            expected(:, [2 5], :) = newCols;
            testCase.verifyEqual(spatialSlabIO('read', fid, 6, 5, 4, 1:5, 'single'), expected);

            testCase.verifyError(@() spatialSlabIO('peek', fid, 6, 5, 4, 1, 'single'), ...
                'spatialSlabIO:InvalidMode');
            testCase.verifyError(@() spatialSlabIO('write', fid, 6, 5, 4, 1, 'single'), ...
                'spatialSlabIO:MissingInput');
            testCase.verifyError(@() spatialSlabIO('write', fid, 6, 5, 4, 1, 'single', zeros(6, 2, 4, 'single')), ...
                'spatialSlabIO:SizeMismatch');
            clear cleanupObj
        end
    end

    methods (Access = private)
        function h = openCase(testCase, c)
            f = fullfile(testCase.FixtureRoot, c.rel{:});
            h = spatialSlabIO('open', f);
            testCase.addTeardown(@() spatialSlabIO('close', h));
        end

        function folder = tempFolder(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end

        function values = knownValues(~, c)
            values = reshape(cast(1:prod(c.dimSizes), c.dataClass), c.dimSizes);
        end

        function hdr = headerFor(~, c)
            hdr = struct('dataClass', c.dataClass, 'frameRateHz', c.frameRateHz, ...
                'exposureMsec', 5, 'channelName', 'written', ...
                'dimNames', {c.dimNames}, 'dimSizes', c.dimSizes);
        end

        function [h, ref] = openWithReference(testCase, c)
            h = testCase.openCase(c);
            [~, ref] = evalc('loadData(h.filePath)');
        end
    end
end

function hdr = iGrowHeader(dataClass)
% Y-X-T header for growable-file tests; dimSizes(end) is ignored.
hdr = struct('dataClass', dataClass, 'frameRateHz', 10, 'exposureMsec', 2, ...
    'channelName', 'grow', 'dimNames', {{'Y', 'X', 'T'}}, 'dimSizes', [3 4 99]);
end
