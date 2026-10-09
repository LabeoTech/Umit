classdef TestGenVSM < matlab.unittest.TestCase
    properties
        ProjectRoot
        TempFolder
        RetinotopyUMT
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;

            testCase.RetinotopyUMT = testCase.buildRetinotopyFixture();
            deleteIfExists(fullfile(testCase.TempFolder, 'VisualCtxAreas.roi'));
        end
    end

    methods (Test)
        function testPipelineInfo(testCase)
            info = genVSM('pipelineInfo');
            testCase.verifyEqual(info.name, 'genVSM');
            testCase.verifyEqual(info.outputs(1).type, {'ProcessedData'});
            dataInput = info.inputs(strcmp({info.inputs.name}, 'retinotopyUMT'));
            testCase.verifyTrue(dataInput.supportsFile);
            testCase.verifyEqual(dataInput.dataMode, 'file');
            testCase.verifyEqual({info.outputs.name}, {'outData','roiFile'});
            testCase.verifyEqual(info.outputs(2).outputMode, 'file');
            testCase.verifyFalse(info.outputs(2).isData);
            testCase.verifyFalse(info.outputs(2).isRequired);
            testCase.verifyTrue(info.outputs(2).returnsValue);
            testCase.verifyEqual(info.outputs(2).defOutfilename, ...
                'VisualCtxAreas.roi');

            phaseSigma = info.parameters(strcmp( ...
                {info.parameters.name}, 'PhaseMapFilter_Sigma'));
            vsmSigma = info.parameters(strcmp( ...
                {info.parameters.name}, 'VSMFilter_Sigma'));
            testCase.verifyEqual(phaseSigma.allowed, [0 Inf]);
            testCase.verifyEqual(vsmSigma.allowed, [0 Inf]);
        end

        function testCreatesVSMEntryAndROIFile(testCase)
            [out, roiFile] = genVSM( ...
                iFile(testCase, testCase.RetinotopyUMT), testCase.TempFolder);
            testCase.verifyTrue(isfield(out.data, 'VSM'));
            testCase.verifyEqual(cellstr(string(out.data.VSM.dimNames)), {'Y','X'});

            testCase.verifyEqual(roiFile, 'VisualCtxAreas.roi');
            roiPath = fullfile(testCase.TempFolder, roiFile);
            testCase.verifyTrue(isfile(roiPath));

            ROIFile = loadROIFile(roiPath);
            testCase.verifyEqual(ROIFile.imageInfo.entryName, 'VSM');
            testCase.verifyEmpty(ROIFile.imageInfo.dataFile);
            testCase.verifyEqual(ROIFile.imageInfo.imageSizeYX, ...
                size(out.data.VSM.value));
            testCase.verifyEqual(ROIFile.statsImage.image, ...
                out.data.VSM.value);
            testCase.verifyNotEmpty(ROIFile.ROIs);

            expectedNames = arrayfun(@(idx) sprintf('CtxArea_%d', idx), ...
                1:numel(ROIFile.ROIs), 'UniformOutput', false);
            testCase.verifyEqual({ROIFile.ROIs.name}, expectedNames);
            for iROI = 1:numel(ROIFile.ROIs)
                roi = ROIFile.ROIs(iROI);
                testCase.verifyTrue(islogical(roi.mask));
                testCase.verifyTrue(any(roi.mask, 'all'));
                testCase.verifyClass(roi.geometry.polyshape, 'polyshape');
                testCase.verifyGreaterThanOrEqual( ...
                    size(roi.geometry.verticesXY_px, 1), 3);
                testCase.verifyEqual(roi.stats.NPixels, nnz(roi.mask));
                testCase.verifyEqual(roi.stats.areaPx2, nnz(roi.mask));
            end
        end

        function testCustomROIFileName(testCase)
            [~, roiFile] = genVSM(iFile(testCase, testCase.RetinotopyUMT), ...
                testCase.TempFolder, 'ROIFileName', 'CustomAreas');
            testCase.verifyEqual(roiFile, 'CustomAreas.roi');
            testCase.verifyTrue(isfile(fullfile(testCase.TempFolder, roiFile)));
        end

        function testNoPatchesWarnsAndLeavesExistingFile(testCase)
            noPatchUMT = iFile(testCase, testCase.buildNoPatchFixture());
            roiPath = fullfile(testCase.TempFolder, 'VisualCtxAreas.roi');
            iWriteSentinelFile(roiPath, 'existing ROI sentinel');

            testCase.verifyWarning( ...
                @() iRunNoPatchGenVSM(noPatchUMT, testCase.TempFolder), ...
                'Umitoolbox:genVSM:NoVisualPatches');
            testCase.verifyEqual(fileread(roiPath), 'existing ROI sentinel');
        end

        function testRemovedPatchOptionsAreRejected(testCase)
            testCase.verifyError(@() genVSM(iFile(testCase, testCase.RetinotopyUMT), ...
                testCase.TempFolder, 'b_CreatePatches', true), ...
                'MATLAB:InputParser:UnmatchedParameter');
            testCase.verifyError(@() genVSM(iFile(testCase, testCase.RetinotopyUMT), ...
                testCase.TempFolder, 'PatchFileName', 'old.mat'), ...
                'MATLAB:InputParser:UnmatchedParameter');
        end

        function testRejectsInvalidROIFileName(testCase)
            testCase.verifyError(@() genVSM(iFile(testCase, testCase.RetinotopyUMT), ...
                testCase.TempFolder, 'ROIFileName', 'areas.mat'), ...
                'Umitoolbox:genVSM:InvalidROIFileName');
        end

        function testBareFileNameIsResolvedInSaveFolder(testCase)
            saveData(fullfile(testCase.TempFolder, 'retinotopyMaps.umt'), testCase.RetinotopyUMT);

            out = genVSM('retinotopyMaps.umt', testCase.TempFolder);

            testCase.verifyTrue(isfield(out.data, 'VSM'));
        end

        function testRejectsUnsupportedInputs(testCase)
            sv = testCase.TempFolder;

            % UMT structs and arrays are not input forms.
            testCase.verifyError(@() genVSM(testCase.RetinotopyUMT, sv), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
            testCase.verifyError(@() genVSM(rand(4, 4, 2), sv), ...
                'MATLAB:InputParser:ArgumentFailedValidation');

            % Missing file and other extensions.
            testCase.verifyError(@() genVSM('missing.umt', sv), ...
                'Umitoolbox:genVSM:InputFileNotFound');
            matFile = fullfile(sv, 'x.mat');
            fclose(fopen(matFile, 'w'));
            testCase.verifyError(@() genVSM(matFile, sv), ...
                'Umitoolbox:genVSM:UnsupportedInputFile');
        end

        function testRejectsInputThatIsNotTheRetinotopyOutput(testCase)
            % Only the AzimuthMap/ElevationMap pair of genRetinotopyMaps.
            extra = genUMTStruct(testCase.RetinotopyUMT, ...
                'value', ones(20, 18, 2, 'single'), 'entryName', 'Other', ...
                'dimNames', {'Y','X','F'});
            testCase.verifyError(@() genVSM(iFile(testCase, extra), testCase.TempFolder), ...
                'Umitoolbox:genVSM:InvalidInput');
        end

        function testRejectsMissingAzimuthEntry(testCase)
            bad = testCase.RetinotopyUMT;
            bad.data = rmfield(bad.data, 'AzimuthMap');
            testCase.verifyError(@() genVSM(iFile(testCase, bad), testCase.TempFolder), ...
                'Umitoolbox:genVSM:MissingInput');
        end

        function testRejectsWrongDims(testCase)
            bad = testCase.RetinotopyUMT;
            bad.data.AzimuthMap.dimNames = {'Y','X','T'};
            testCase.verifyError(@() genVSM(iFile(testCase, bad), testCase.TempFolder), ...
                'Umitoolbox:genVSM:InvalidInput');
        end

        function testRejectsWrongFLength(testCase)
            bad = testCase.RetinotopyUMT;
            bad.data.ElevationMap.value = bad.data.ElevationMap.value(:,:,1);
            bad.data.ElevationMap.dimNames = {'Y','X'};
            testCase.verifyError(@() genVSM(iFile(testCase, bad), testCase.TempFolder), ...
                'Umitoolbox:genVSM:InvalidInput');
        end

        function testAsymmetricNaNMasksInPhaseMaps(testCase)
            % Each phase map must be masked against its OWN NaNs. Regression
            % guard for the copy/paste defect where the elevation map was
            % masked with isnan(phaseAz), which both left elevation NaNs
            % unreplaced and clobbered valid elevation values wherever the
            % azimuth map happened to be NaN.
            azNaN = [3 4];   % NaN in azimuth phase only
            elNaN = [9 11];  % NaN in elevation phase only

            withNaN = testCase.RetinotopyUMT;
            withNaN.data.AzimuthMap.value(azNaN(1), azNaN(2), 2) = NaN;
            withNaN.data.ElevationMap.value(elNaN(1), elNaN(2), 2) = NaN;

            % Same fixture with the intended substitution already applied.
            preSubstituted = testCase.RetinotopyUMT;
            preSubstituted.data.AzimuthMap.value(azNaN(1), azNaN(2), 2) = 1000;
            preSubstituted.data.ElevationMap.value(elNaN(1), elNaN(2), 2) = 1000;

            actual = genVSM(iFile(testCase, withNaN), testCase.TempFolder);
            expected = genVSM(iFile(testCase, preSubstituted), testCase.TempFolder);

            testCase.verifyEqual(actual.data.VSM.value, expected.data.VSM.value);

            % The elevation NaN must actually have been substituted rather
            % than propagated: a surviving NaN would be zeroed by the
            % vsm(isnan(vsm)) = 0 guard, so assert the result is finite.
            testCase.verifyTrue(all(isfinite(actual.data.VSM.value), 'all'));
        end
    end

    methods (Test)
        function testPhaseSmoothingIsCircularAcrossWrapBoundary(testCase)
            % Azimuth phase ramps through the 0/2*pi wrap (2*pi-1 ... 2*pi+1
            % rewrapped); elevation is a clean ramp, so the VSM is -1 wherever
            % the phase gradients are sane. Smoothing wrapped phase directly
            % averages 6.2 with 0.1 into a steep reversed ramp that flips the
            % VSM sign over about 2*sigma columns on each side of the wrap
            % (10 columns here). Smoothing through exp(1i*phase) leaves only
            % the unavoidable one-pixel discontinuity (2 columns).
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'Umitoolbox:genVSM:NoVisualPatches'));
            Ny = 30;
            Nx = 40;
            azPhase = repmat(single(mod(linspace(2*pi - 1, 2*pi + 1, Nx), 2*pi)), Ny, 1);
            elPhase = repmat(single(pi + linspace(-1, 1, Ny)'), 1, Nx);
            umt = iBuildPhaseUMT(azPhase, elPhase);

            out = genVSM(iFile(testCase, umt), testCase.TempFolder, ...
                'PhaseMapFilter_Sigma', 2);
            vsm = out.data.VSM.value;

            rows = 8:(Ny - 7);
            columnMean = mean(vsm(rows, :), 1);
            testCase.verifyLessThanOrEqual(nnz(columnMean > 0), 3, ...
                'Phase smoothing must not blur across the 0/2*pi wrap.');
            awayFromWrap = [1:12, 28:Nx];
            testCase.verifyEqual(double(mean(vsm(rows, awayFromWrap), 'all')), -1, ...
                'AbsTol', 1e-4);
        end

        function testDegreePhaseMapsAreSmoothedDirectly(testCase)
            % Maps calibrated to degrees of visual angle leave [0, 2*pi] and
            % have no known period, so they must be filtered as plain values.
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'Umitoolbox:genVSM:NoVisualPatches'));
            Ny = 24;
            Nx = 28;
            [X, Y] = meshgrid(linspace(-1, 1, Nx), linspace(-1, 1, Ny));
            azPhase = single(20 * X + 3 * Y);    % degrees, negative and positive
            elPhase = single(15 * Y - 4 * X);
            umt = iBuildPhaseUMT(azPhase, elPhase);
            sigma = 1.5;

            out = genVSM(iFile(testCase, umt), testCase.TempFolder, ...
                'PhaseMapFilter_Sigma', sigma);

            [gax, gay] = gradient(imgaussfilt(azPhase, sigma));
            [gex, gey] = gradient(imgaussfilt(elPhase, sigma));
            expected = sin(angle(exp(1i .* atan2(gay, gax)) .* ...
                exp(-1i .* atan2(gey, gex))));
            testCase.verifyEqual(double(out.data.VSM.value), double(expected), ...
                'AbsTol', 1e-5);
        end
    end

    methods (Access = private)
        function umt = buildRetinotopyFixture(~)
            Ny = 20;
            Nx = 18;
            az = zeros(Ny, Nx, 2, 'single');
            el = zeros(Ny, Nx, 2, 'single');

            [X, Y] = meshgrid(linspace(-1,1,Nx), linspace(-1,1,Ny));
            az(:,:,1) = single(abs(X));
            az(:,:,2) = single(pi + X);
            el(:,:,1) = single(abs(Y));
            el(:,:,2) = single(pi + Y);

            labels = struct();
            labels.F = {'Amplitude','Phase'};

            umt = genUMTStruct( ...
                az, ...
                'kind', 'image', ...
                'entryName', 'AzimuthMap', ...
                'dimNames', {'Y','X','F'}, ...
                'labels', labels);

            umt = genUMTStruct( ...
                umt, ...
                'value', el, ...
                'entryName', 'ElevationMap', ...
                'dimNames', {'Y','X','F'});
        end

        function umt = buildNoPatchFixture(testCase)
            umt = testCase.RetinotopyUMT;
            phase = umt.data.AzimuthMap.value(:,:,2);
            umt.data.ElevationMap.value(:,:,2) = phase;
        end
    end
end

function umt = iBuildPhaseUMT(azPhase, elPhase)
%IBUILDPHASEUMT Retinotopy UMT with unit amplitude and the given phase maps.
labels = struct();
labels.F = {'Amplitude','Phase'};
az = cat(3, ones(size(azPhase), 'single'), single(azPhase));
el = cat(3, ones(size(elPhase), 'single'), single(elPhase));
umt = genUMTStruct(az, 'kind', 'image', 'entryName', 'AzimuthMap', ...
    'dimNames', {'Y','X','F'}, 'labels', labels);
umt = genUMTStruct(umt, 'value', el, 'entryName', 'ElevationMap', ...
    'dimNames', {'Y','X','F'});
end

function umtFile = iFile(testCase, umt)
%IFILE Save a UMT struct as a .umt file in the test folder and return its path.
umtFile = [tempname(testCase.TempFolder) '.umt'];
saveData(umtFile, umt);
end

function iRunNoPatchGenVSM(umtFile, saveFolder)
[~, roiFile] = genVSM(umtFile, saveFolder);
assert(isempty(roiFile), ...
    'Umitoolbox:TestGenVSM:UnexpectedROIFile', ...
    'Zero-patch genVSM execution returned an ROI filename.');
end

function iWriteSentinelFile(filePath, contents)
fid = fopen(filePath, 'w');
assert(fid ~= -1, 'Failed to create sentinel ROI file.');
cleanupObj = onCleanup(@() fclose(fid));
fwrite(fid, contents, 'char');
end

function deleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end
