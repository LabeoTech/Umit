classdef testIO < matlab.unittest.TestCase
    %TESTIO Unit tests for UMT data packaging and raw/UMT I/O.
    %
    %   This suite covers:
    %       - getUMTSchema
    %       - genUMTStruct
    %       - validateUMTStruct
    %       - saveData
    %       - loadData
    %
    %   The tests validate the current .umt schema by checking:
    %       - all allowed dimension combinations for kind='image'
    %       - all allowed dimension combinations for kind='roi'
    %       - invalid dimension combinations
    %       - row-vector rejection
    %       - labels validation
    %       - append and overwrite behavior
    %       - representative .umt round-trip save/load
    %       - raw .dat round-trip save/load
    %       - AcqInfos.mat base/imported timeline validation
    %       - legacy sidecar metadata compatibility
    %       - defensive validation edge cases
    %
    %   To run:
    %       results = runtests('testIO')

    properties (Constant)
        % Set this to a folder containing .mat files with 3D arrays if desired.
        % Leave as "" to use automatic fallback behavior.
        SourceDataFolder = "D:\umit-dev\SaveFolderForDev"

        % Temporary folder name prefix used during test execution.
        TempFolderPrefix = 'tUMTIO_'
    end

    properties
        TempFolder
        BaseData    % single 3D array
        BaseName    % descriptive name of selected source
    end

    methods (TestMethodSetup)
        function createTempFolderAndLoadData(testCase)
            %CREATETEMPFOLDERANDLOADDATA Create temp folder and get base 3D data.

            testCase.TempFolder = TestIOHelpers.makeTempFolder(testCase.TempFolderPrefix);

            [testCase.BaseData, testCase.BaseName] = ...
                TestIOHelpers.getTest3DArray(testCase.SourceDataFolder);
        end
    end

    methods (TestMethodTeardown)
        function removeTempFolder(testCase)
            %REMOVETEMPFOLDER Clean up temporary folder after each test.

            TestIOHelpers.removeFolderIfExists(testCase.TempFolder);
        end
    end

    methods (Test)
        function testGetUMTSchemaVersion1(testCase)
            %TESTGETUMTSCHEMAVERSION1 Validate centralized schema contents.

            schema = getUMTSchema(1);

            testCase.verifyEqual(schema.version, 1);
            testCase.verifyEqual(schema.allowedKinds, {'image', 'roi'});
            testCase.verifyEqual(schema.allowedDims, {'Y', 'X', 'T', 'E', 'F', 'ROI', 'Pixel', 'Measure'});
            testCase.verifyTrue(isfield(schema.allowedPatterns, 'image'));
            testCase.verifyTrue(isfield(schema.allowedPatterns, 'roi'));
        end

        function testGetUMTSchemaRejectsUnsupportedVersion(testCase)
            %TESTGETUMTSCHEMAREJECTSUNSUPPORTEDVERSION Reject unknown schema versions.

            testCase.verifyError(@() getUMTSchema(999), ...
                'Umitoolbox:getUMTSchema:unsupportedVersion');
        end

        function testValidateUMTStructRejectsUnsupportedVersion(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSUNSUPPORTEDVERSION Reject UMT with unsupported version.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);

            S = genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels);

            S.version = 999;

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:getUMTSchema:unsupportedVersion');
        end

        function testGenUMTStructAcceptsMixedCaseImageInputs(testCase)
            %TESTGENUMTSTRUCTACCEPTSMIXEDCASEIMAGEINPUTS Accept valid mixed-case image inputs.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X','E'}, testCase.BaseData);
            labelsLower = struct('e', {labels.E});

            S = genUMTStruct(value, ...
                'kind', 'ImAgE', ...
                'entryName', 'mixedCaseImage', ...
                'dimNames', {'y','x','e'}, ...
                'labels', labelsLower);

            testCase.verifyEqual(S.kind, 'image');
            testCase.verifyEqual(S.data.mixedCaseImage.dimNames, {'Y','X','E'});
            testCase.verifyTrue(isfield(S.labels, 'E'));
            testCase.verifyWarningFree(@() validateUMTStruct(S, 'requireEventInfo', false));

            S = TestIOHelpers.appendDefaultEventInfoIfNeeded(S);
            testCase.verifyWarningFree(@() validateUMTStruct(S));
        end

        function testGenUMTStructAcceptsMixedCaseROIInputs(testCase)
            %TESTGENUMTSTRUCTACCEPTSMIXEDCASEROIINPUTS Accept valid mixed-case roi inputs.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI','Measure'}, testCase.BaseData);
            labelsLower = struct('roi', {labels.ROI}, 'measure', {labels.Measure});

            S = genUMTStruct(value, ...
                'kind', 'RoI', ...
                'entryName', 'mixedCaseROI', ...
                'dimNames', {'roi','measure'}, ...
                'labels', labelsLower);

            testCase.verifyEqual(S.kind, 'roi');
            testCase.verifyEqual(S.data.mixedCaseROI.dimNames, {'ROI','Measure'});
            testCase.verifyTrue(isfield(S.labels, 'ROI'));
            testCase.verifyTrue(isfield(S.labels, 'Measure'));
            testCase.verifyWarningFree(@() validateUMTStruct(S));
        end

        function testValidateUMTStructAcceptsAllImagePatterns(testCase)
            %TESTVALIDATEUMTSTRUCTACCEPTSALLIMAGEPATTERNS Accept every valid image layout.

            schema = getUMTSchema(1);
            imagePatterns = [{ {} }, schema.allowedPatterns.image];

            for iPattern = 1:numel(imagePatterns)
                dimNames = imagePatterns{iPattern};
                [value, labels] = TestIOHelpers.makeUMTPayload('image', dimNames, testCase.BaseData);

                S = genUMTStruct(value, ...
                    'kind', 'image', ...
                    'entryName', sprintf('imageEntry_%d', iPattern), ...
                    'dimNames', dimNames, ...
                    'labels', labels);

                testCase.verifyEqual(S.version, 1);
                testCase.verifyEqual(S.kind, 'image');
                testCase.verifyTrue(isfield(S.data, sprintf('imageEntry_%d', iPattern)));

                entry = S.data.(sprintf('imageEntry_%d', iPattern));
                testCase.verifyEqual(entry.dimNames, dimNames);
                testCase.verifyWarningFree(@() validateUMTStruct(S, 'requireEventInfo', false));

                S = TestIOHelpers.appendDefaultEventInfoIfNeeded(S);
                testCase.verifyWarningFree(@() validateUMTStruct(S));
            end
        end

        function testValidateUMTStructAcceptsAllROIPatterns(testCase)
            %TESTVALIDATEUMTSTRUCTACCEPTSALLROIPATTERNS Accept every valid roi layout.

            schema = getUMTSchema(1);
            roiPatterns = [{ {} }, schema.allowedPatterns.roi];

            for iPattern = 1:numel(roiPatterns)
                dimNames = roiPatterns{iPattern};
                [value, labels] = TestIOHelpers.makeUMTPayload('roi', dimNames, testCase.BaseData);

                S = genUMTStruct(value, ...
                    'kind', 'roi', ...
                    'entryName', sprintf('roiEntry_%d', iPattern), ...
                    'dimNames', dimNames, ...
                    'labels', labels);

                testCase.verifyEqual(S.version, 1);
                testCase.verifyEqual(S.kind, 'roi');
                testCase.verifyTrue(isfield(S.data, sprintf('roiEntry_%d', iPattern)));

                entry = S.data.(sprintf('roiEntry_%d', iPattern));
                testCase.verifyEqual(entry.dimNames, dimNames);
                testCase.verifyWarningFree(@() validateUMTStruct(S, 'requireEventInfo', false));

                S = TestIOHelpers.appendDefaultEventInfoIfNeeded(S);
                testCase.verifyWarningFree(@() validateUMTStruct(S));
            end
        end

        function testGenUMTStructAppendsImageMeasurements(testCase)
            %TESTGENUMTSTRUCTAPPENDSIMAGEMEASUREMENTS Append a second image entry.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);
            [value2, labels2] = TestIOHelpers.makeUMTPayload('image', {'Y','X','E'}, testCase.BaseData);

            S = genUMTStruct(value1, ...
                'kind', 'image', ...
                'entryName', 'map', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels1);

            S = genUMTStruct(S, ...
                'value', value2, ...
                'entryName', 'eventMap', ...
                'dimNames', {'Y','X','E'}, ...
                'labels', labels2);

            testCase.verifyTrue(isfield(S.data, 'map'));
            testCase.verifyTrue(isfield(S.data, 'eventMap'));
            testCase.verifyEqual(S.data.map.value, value1);
            testCase.verifyEqual(S.data.eventMap.value, value2);
            testCase.verifyWarningFree(@() validateUMTStruct(S, 'requireEventInfo', false));

            S = TestIOHelpers.appendDefaultEventInfoIfNeeded(S);
            testCase.verifyWarningFree(@() validateUMTStruct(S));
        end

        function testGenUMTStructAppendsROIMeasurements(testCase)
            %TESTGENUMTSTRUCTAPPENDSROIMEASUREMENTS Append multiple roi entries.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('roi', {'ROI','T'}, testCase.BaseData);
            [value2, labels2] = TestIOHelpers.makeUMTPayload('roi', {'ROI','Measure'}, testCase.BaseData);

            S = genUMTStruct(value1, ...
                'kind', 'roi', ...
                'entryName', 'trace', ...
                'dimNames', {'ROI','T'}, ...
                'labels', labels1);

            S = genUMTStruct(S, ...
                'value', value2, ...
                'entryName', 'stats', ...
                'dimNames', {'ROI','Measure'}, ...
                'labels', labels2);

            testCase.verifyTrue(isfield(S.data, 'trace'));
            testCase.verifyTrue(isfield(S.data, 'stats'));
            testCase.verifyEqual(S.data.trace.value, value1);
            testCase.verifyEqual(S.data.stats.value, value2);
            testCase.verifyWarningFree(@() validateUMTStruct(S));
        end

        function testGenUMTStructRejectsDuplicateEntryWithoutOverwrite(testCase)
            %TESTGENUMTSTRUCTREJECTSDUPLICATEENTRYWITHOUTOVERWRITE Reject duplicate entry names.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('roi', {'ROI'}, testCase.BaseData);
            [value2, ~] = TestIOHelpers.makeUMTPayload('roi', {'ROI'}, testCase.BaseData);

            S = genUMTStruct(value1, ...
                'kind', 'roi', ...
                'entryName', 'metric', ...
                'dimNames', {'ROI'}, ...
                'labels', labels1);

            testCase.verifyError(@() genUMTStruct(S, ...
                'value', value2 + 10, ...
                'entryName', 'metric', ...
                'dimNames', {'ROI'}), ...
                'Umitoolbox:genUMTStruct:invalidInput');
        end

        function testGenUMTStructOverwriteExistingEntry(testCase)
            %TESTGENUMTSTRUCTOVERWRITEEXISTINGENTRY Overwrite an existing entry.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('roi', {'ROI'}, testCase.BaseData);
            [value2, ~] = TestIOHelpers.makeUMTPayload('roi', {'ROI'}, testCase.BaseData);
            value2 = value2 + 100;

            S = genUMTStruct(value1, ...
                'kind', 'roi', ...
                'entryName', 'metric', ...
                'dimNames', {'ROI'}, ...
                'labels', labels1);

            S = genUMTStruct(S, ...
                'value', value2, ...
                'entryName', 'metric', ...
                'dimNames', {'ROI'}, ...
                'overwrite', true);

            testCase.verifyEqual(S.data.metric.value, value2);
            testCase.verifyWarningFree(@() validateUMTStruct(S));
        end

        function testGenUMTStructOverwriteLabelsSucceeds(testCase)
            %TESTGENUMTSTRUCTOVERWRITELABELSSUCCEEDS Overwrite conflicting shared labels.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('roi', {'ROI','T'}, testCase.BaseData);
            [value2, labels2] = TestIOHelpers.makeUMTPayload('roi', {'ROI','Measure'}, testCase.BaseData);

            labels2.ROI = {'R_A', 'R_B', 'R_C', 'R_D'};

            S = genUMTStruct(value1, ...
                'kind', 'roi', ...
                'entryName', 'trace', ...
                'dimNames', {'ROI','T'}, ...
                'labels', labels1);

            S = genUMTStruct(S, ...
                'value', value2, ...
                'entryName', 'stats', ...
                'dimNames', {'ROI','Measure'}, ...
                'labels', labels2, ...
                'overwrite', true);

            testCase.verifyEqual(S.labels.ROI, labels2.ROI);
            testCase.verifyEqual(S.labels.Measure, labels2.Measure);
            testCase.verifyWarningFree(@() validateUMTStruct(S));
        end

        function testGenUMTStructRejectsInvalidImagePattern(testCase)
            %TESTGENUMTSTRUCTREJECTSINVALIDIMAGEPATTERN Reject non-image dims for image kind.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI','T'}, testCase.BaseData);

            testCase.verifyError(@() genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'badImage', ...
                'dimNames', {'ROI','T'}, ...
                'labels', labels), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testGenUMTStructRejectsInvalidROIPattern(testCase)
            %TESTGENUMTSTRUCTREJECTSINVALIDROIPATTERN Reject image dims for roi kind.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);

            testCase.verifyError(@() genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'badROI', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsImageRowVector(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSIMAGEROWVECTOR Reject 1D row-vector data.

            value = single(1:7);

            testCase.verifyError(@() genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'badRow', ...
                'dimNames', {'T'}), ...
                'Umitoolbox:genUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsROIRowVector(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSROIROWVECTOR Reject roi row-vector data.

            value = single(1:4);

            testCase.verifyError(@() genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'badRow', ...
                'dimNames', {'ROI'}), ...
                'Umitoolbox:genUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsMissingTopField(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSMISSINGTOPFIELD Reject missing required field.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels);

            S = rmfield(S, 'data');

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsUnsupportedTopField(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSUNSUPPORTEDTOPFIELD Reject unknown fields.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI','Measure'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'stats', ...
                'dimNames', {'ROI','Measure'}, ...
                'labels', labels);

            S.unexpectedField = 1;

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsMissingEntryValue(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSMISSINGENTRYVALUE Reject entry missing value.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'main', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels);

            S.data.main = rmfield(S.data.main, 'value');

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsMissingEntryDimNames(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSMISSINGENTRYDIMNAMES Reject entry missing dimNames.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI','T'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'trace', ...
                'dimNames', {'ROI','T'}, ...
                'labels', labels);

            S.data.trace = rmfield(S.data.trace, 'dimNames');

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsNonScalarEntryStruct(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSNONSCALARENTRYSTRUCT Reject non-scalar entry structs.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'metric', ...
                'dimNames', {'ROI'}, ...
                'labels', labels);

            S.data.metric = repmat(S.data.metric, 1, 2);

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsUnsupportedEntryField(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSUNSUPPORTEDEXTRYFIELD Reject unknown entry-level fields.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X','T'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'movie', ...
                'dimNames', {'Y','X','T'}, ...
                'labels', labels);

            S.data.movie.foo = 1;

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsInvalidPayloadClassChar(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSINVALIDPAYLOADCLASSCHAR Reject char payloads.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'metric', ...
                'dimNames', {'ROI'}, ...
                'labels', labels);

            S.data.metric.value = 'abc';

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsInvalidPayloadClassString(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSINVALIDPAYLOADCLASSSTRING Reject string payloads.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'metric', ...
                'dimNames', {'ROI'}, ...
                'labels', labels);

            S.data.metric.value = "abc";

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsInvalidPayloadClassCell(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSINVALIDPAYLOADCLASSCELL Reject cell payloads.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'map', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels);

            S.data.map.value = {1, 2, 3};

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsInvalidPayloadClassStruct(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSINVALIDPAYLOADCLASSSTRUCT Reject struct payloads.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'map', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels);

            S.data.map.value = struct('bad', 1);

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsWrongLabelLength(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSWRONGLABELLENGTH Reject mismatched label lengths.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI','Measure'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'stats', ...
                'dimNames', {'ROI','Measure'}, ...
                'labels', labels);

            S.labels.ROI = {'ROI_01', 'ROI_02'}; % Wrong length.

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsUnusedLabelDimension(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSUNUSEDLABELDIMENSION Reject labels for unused dims.

            [value, labels] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'image', ...
                'entryName', 'map', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels);

            S.labels.E = {'Event_1', 'Event_2'};

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsDuplicateLabelFieldsAfterCaseNormalization(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSDUPLICATELABELFIELDSAFTERCASENORMALIZATION Reject duplicate label dims.

            [value, labels] = TestIOHelpers.makeUMTPayload('roi', {'ROI'}, testCase.BaseData);
            S = genUMTStruct(value, ...
                'kind', 'roi', ...
                'entryName', 'metric', ...
                'dimNames', {'ROI'}, ...
                'labels', labels);

            S.labels.roi = S.labels.ROI;

            testCase.verifyError(@() validateUMTStruct(S), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsCrossEntrySharedDimensionSizeMismatch(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSCROSSENTRYSHAREDDIMENSIONSIZEMISMATCH Reject inconsistent ROI size across entries.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('roi', {'ROI','T'}, testCase.BaseData);

            value2 = reshape(single(1:6), [3, 2]); % 3 ROIs instead of 4

            S = genUMTStruct(value1, ...
                'kind', 'roi', ...
                'entryName', 'trace', ...
                'dimNames', {'ROI','T'}, ...
                'labels', labels1);

            testCase.verifyError(@() genUMTStruct(S, ...
                'value', value2, ...
                'entryName', 'stats', ...
                'dimNames', {'ROI','Measure'}), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testValidateUMTStructRejectsConflictingSharedLabels(testCase)
            %TESTVALIDATEUMTSTRUCTREJECTSCONFLICTINGSHAREDLABELS Reject inconsistent top-level labels.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('roi', {'ROI','T'}, testCase.BaseData);
            [value2, labels2] = TestIOHelpers.makeUMTPayload('roi', {'ROI','Measure'}, testCase.BaseData);

            S = genUMTStruct(value1, ...
                'kind', 'roi', ...
                'entryName', 'trace', ...
                'dimNames', {'ROI','T'}, ...
                'labels', labels1);

            labels2.ROI = {'Other_1', 'Other_2', 'Other_3', 'Other_4'};

            testCase.verifyError(@() genUMTStruct(S, ...
                'value', value2, ...
                'entryName', 'stats', ...
                'dimNames', {'ROI','Measure'}, ...
                'labels', labels2), ...
                'Umitoolbox:genUMTStruct:invalidInput');
        end

        function testSaveAndLoadUMTImageRoundTrip(testCase)
            %TESTSAVEANDLOADUMTIMAGEROUNDTRIP Round-trip save/load for image-kind UMT.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('image', {'Y','X'}, testCase.BaseData);
            [value2, labels2] = TestIOHelpers.makeUMTPayload('image', {'Y','X','E'}, testCase.BaseData);

            S = genUMTStruct(value1, ...
                'kind', 'image', ...
                'entryName', 'map', ...
                'dimNames', {'Y','X'}, ...
                'labels', labels1);

            S = genUMTStruct(S, ...
                'value', value2, ...
                'entryName', 'eventMap', ...
                'dimNames', {'Y','X','E'}, ...
                'labels', labels2);
            S = TestIOHelpers.appendDefaultEventInfoIfNeeded(S);

            filePath = fullfile(testCase.TempFolder, 'image_processed.randomExt');
            saveData(filePath, S);

            savedFile = fullfile(testCase.TempFolder, 'image_processed.umt');
            testCase.verifyTrue(isfile(savedFile));

            [outFile, Info] = loadData(savedFile);

            testCase.verifyEqual(Info.FileType, '.umt');
            testCase.verifyEqual(outFile, S);
        end

        function testSaveAndLoadUMTROIRoundTrip(testCase)
            %TESTSAVEANDLOADUMTROIROUNDTRIP Round-trip save/load for roi-kind UMT.

            [value1, labels1] = TestIOHelpers.makeUMTPayload('roi', {'ROI','T'}, testCase.BaseData);
            [value2, labels2] = TestIOHelpers.makeUMTPayload('roi', {'ROI','Measure'}, testCase.BaseData);

            S = genUMTStruct(value1, ...
                'kind', 'roi', ...
                'entryName', 'trace', ...
                'dimNames', {'ROI','T'}, ...
                'labels', labels1);

            S = genUMTStruct(S, ...
                'value', value2, ...
                'entryName', 'stats', ...
                'dimNames', {'ROI','Measure'}, ...
                'labels', labels2);

            filePath = fullfile(testCase.TempFolder, 'roi_processed.ext');
            saveData(filePath, S);

            savedFile = fullfile(testCase.TempFolder, 'roi_processed.umt');
            testCase.verifyTrue(isfile(savedFile));

            [outFile, Info] = loadData(savedFile);

            testCase.verifyEqual(Info.FileType, '.umt');
            testCase.verifyEqual(outFile, S);
        end

        function testAppendImportedChannelInfoAddsDefaultCamIdx(testCase)
            %TESTAPPENDIMPORTEDCHANNELINFOADDSDEFAULTCAMIDX Add CamIdx for non-multicamera data.

            data = testCase.BaseData;
            AcqInfoStream = TestIOHelpers.makeAcqInfoStreamForData(data, 10);

            channelInfo = TestIOHelpers.makeImportedChannelInfo( ...
                'red.dat', ...
                size(data,3), ...
                10, ...
                5);

            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, channelInfo);

            testCase.verifyTrue(isfield(AcqInfoStream, 'ImportedChannels'));
            testCase.verifyEqual(numel(AcqInfoStream.ImportedChannels), 1);
            testCase.verifyEqual(AcqInfoStream.ImportedChannels(1).DatFile, 'red.dat');
            testCase.verifyEqual(AcqInfoStream.ImportedChannels(1).CamIdx, 1);
        end

        function testSaveDATAndLoadDATExplicitRate(testCase)
            %TESTSAVEDATANDLOADDATEXPLICITRATE Round-trip .dat with explicit axes and rate.
            %
            %   saveData takes the axes and frame rate from its arguments and
            %   never writes AcqInfos.mat (.dat header Phase 6a).

            data = testCase.BaseData;
            [Ny, Nx, Nt] = size(data);

            filePath = fullfile(testCase.TempFolder, 'rawdata.whatever');
            saveData(filePath, data, 'DimNames', {'Y','X','T'}, 'FrameRateHz', 12.5);

            savedFile = fullfile(testCase.TempFolder, 'rawdata.dat');
            testCase.verifyTrue(isfile(savedFile));
            testCase.verifyFalse(isfile(fullfile(testCase.TempFolder, 'AcqInfos.mat')));

            [outFile, Info] = loadData(savedFile);

            testCase.verifyEqual(size(outFile), size(data));
            testCase.verifyEqual(outFile, data);

            % .dat Info schema (old names removed in .dat header Phase 8a).
            testCase.verifyEqual(Info.dimSizes, [Ny, Nx, Nt]);
            testCase.verifyEqual(Info.frameRateHz, 12.5);
            testCase.verifyEqual(Info.dataClass, 'single');
            testCase.verifyEqual(Info.dimNames, {'Y','X','T'});
            testCase.verifyEqual(Info.filePath, savedFile);

            % loadMetaData should return file-facing metadata, not the full
            % acquisition structure.
            testCase.verifyFalse(isfield(Info, 'fileName'));
            testCase.verifyFalse(isfield(Info, 'OriginalLength'));
            testCase.verifyFalse(isfield(Info, 'ImportedChannels'));
            testCase.verifyFalse(isfield(Info, 'AISampleRate'));
            testCase.verifyFalse(isfield(Info, 'Illumination1'));
            testCase.verifyFalse(isfield(Info, 'Tag'));
            testCase.verifyFalse(isfield(Info, 'Color'));
        end

        function testSaveDATIgnoresAcqInfosTimeline(testCase)
            %TESTSAVEDATIGNORESACQINFOSTIMELINE No frame-rate fallback to AcqInfos.mat.
            %
            %   AcqInfos.mat describes the raw acquisition; temporal binning
            %   at import changes the imported rate, so saveData does not
            %   guess it from there (.dat header Phase 6a): without
            %   'FrameRateHz' or Info it errors, even with a matching
            %   ImportedChannels entry.

            baseData = testCase.BaseData;
            [Ny, Nx, Nt] = size(baseData);

            redData = cat(3, baseData, baseData + max(baseData(:)));
            redLength = size(redData, 3);

            AcqInfoStream = TestIOHelpers.makeAcqInfoStreamForData(baseData, 10);
            AcqInfoStream = appendImportedChannelInfo(AcqInfoStream, [ ...
                TestIOHelpers.makeImportedChannelInfo('red.dat', redLength, 20, 4), ...
                TestIOHelpers.makeImportedChannelInfo('green.dat', Nt, 10, 6)]);
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');

            filePath = fullfile(testCase.TempFolder, 'red.randomExt');
            savedFile = fullfile(testCase.TempFolder, 'red.dat');
            testCase.verifyError(@() saveData(filePath, redData, 'DimNames', {'Y','X','T'}), ...
                'Umitoolbox:saveData:missingFrameRate');
            testCase.verifyFalse(isfile(savedFile));

            saveData(filePath, redData, 'DimNames', {'Y','X','T'}, 'FrameRateHz', 20, ...
                'Info', struct('exposureMsec', 4));
            [outFile, Info] = loadData(savedFile);

            testCase.verifyEqual(size(outFile), [Ny, Nx, redLength]);
            testCase.verifyEqual(outFile, redData);
            testCase.verifyEqual(datAxisSize(Info, 'T'), redLength);

            % saveData writes a header with the given rate and exposure.
            testCase.verifyEqual(Info.format, 'header');
            testCase.verifyEqual(Info.frameRateHz, 20);
            testCase.verifyEqual(Info.exposureMsec, 4);

            testCase.verifyFalse(isfield(Info, 'fileName'));
            testCase.verifyFalse(isfield(Info, 'OriginalLength'));
            testCase.verifyFalse(isfield(Info, 'ImportedChannels'));
            testCase.verifyFalse(isfield(Info, 'Tag'));
            testCase.verifyFalse(isfield(Info, 'Color'));
        end

        function testSaveDATAcceptsUnknownTimeline(testCase)
            %TESTSAVEDATACCEPTSUNKNOWNTIMELINE A headered .dat may have any T.
            %
            %   Before the .dat header, saveData rejected a T that matched no
            %   AcqInfos timeline; the header now describes the file itself.

            data = testCase.BaseData;

            otherData = data(:,:,1:end-1);
            filePath = fullfile(testCase.TempFolder, 'otherTimeline');
            saveData(filePath, otherData, 'DimNames', {'Y','X','T'}, 'FrameRateHz', 10);

            savedFile = fullfile(testCase.TempFolder, 'otherTimeline.dat');
            testCase.verifyTrue(isDatWithHeader(savedFile));
            [outData, Info] = loadData(savedFile);
            testCase.verifyEqual(outData, otherData);
            testCase.verifyEqual(Info.frameRateHz, 10);
        end

        function testLoadDATRejectsAcqInfosBoundFile(testCase)
            %TESTLOADDATREJECTSACQINFOSBOUNDFILE Headerless .dat without sidecar is rejected.
            %
            %   Even when AcqInfos.mat fully describes it (matching size and
            %   timeline), a headerless file without a legacy sidecar is no
            %   longer readable (.dat header Phase 5b).

            data = testCase.BaseData;
            AcqInfoStream = TestIOHelpers.makeAcqInfoStreamForData(data, 10);
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');

            boundFile = fullfile(testCase.TempFolder, 'acqInfosBound.dat');
            fid = fopen(boundFile, 'w');
            fwrite(fid, data, 'single');
            fclose(fid);

            id = 'Umitoolbox:loadMetaData:acqInfosBoundUnsupported';
            testCase.verifyError(@() loadMetaData(boundFile), id);
            testCase.verifyError(@() loadData(boundFile), id);
        end

        function testLoadDATUsesLegacySidecarMetadata(testCase)
            %TESTLOADDATUSESLEGACYSIDECARMETADATA Preserve Astrocyte-style sidecar compatibility.

            data = testCase.BaseData(:,:,1:end-1);
            [Ny, Nx, Nt] = size(data);

            datFile = fullfile(testCase.TempFolder, 'legacy_red.dat');
            fid = fopen(datFile, 'w');
            fwrite(fid, data, 'single');
            fclose(fid);

            datSize = [Ny, Nx]; %#ok<NASGU>
            datLength = Nt; %#ok<NASGU>
            Freq = 7.5; %#ok<NASGU>
            Datatype = 'single'; %#ok<NASGU>
            dim_names = {'Y','X','T'}; %#ok<NASGU>
            save(fullfile(testCase.TempFolder, 'legacy_red.mat'), ...
                'datSize', 'datLength', 'Freq', 'Datatype', 'dim_names');

            [outFile, Info] = loadData(datFile);

            testCase.verifyEqual(outFile, data);
            testCase.verifyEqual(Info.dimSizes, [Ny, Nx, Nt]);
            testCase.verifyEqual(Info.frameRateHz, Freq);
            testCase.verifyEqual(Info.format, 'legacySidecar');
        end

        function testLoadDATDetectsLegacyEventSplitAsUnsupported(testCase)
            %TESTLOADDATDETECTSLEGACYEVENTSPLITASUNSUPPORTED Legacy 2-and-2 split raises
            % a clear, specific warning-style error instead of the generic
            % invalidDatSize error. See task-legacy-split-error-messaging.

            datFile = fullfile(testCase.TempFolder, 'legacy_events.dat');
            fid = fopen(datFile, 'w');
            fwrite(fid, zeros(1, 10, 'single'), 'single');
            fclose(fid);

            % Astrocyte's pre-fix split_data_by_event.m convention: dim_names =
            % {E,Y,X,T}, with datSize/datLength splitting the 4 dims 2-and-2
            % (datSize = [nTrials, Height], datLength = [Width, trialLen]).
            dim_names = {'E','Y','X','T'}; %#ok<NASGU>
            datSize = [3, 112]; %#ok<NASGU>
            datLength = [112, 320]; %#ok<NASGU>
            Freq = 10; %#ok<NASGU>
            Datatype = 'single'; %#ok<NASGU>
            save(fullfile(testCase.TempFolder, 'legacy_events.mat'), ...
                'dim_names', 'datSize', 'datLength', 'Freq', 'Datatype');

            testCase.verifyError(@() loadData(datFile), ...
                'Umitoolbox:loadMetaData:legacyEventSplitUnsupported');

            try
                loadData(datFile);
                messageText = '';
            catch ME
                messageText = ME.message;
            end
            testCase.verifySubstring(messageText, 'not currently supported');
            testCase.verifySubstring(messageText, 'no automated way to convert');
        end

        function testLoadDATDetectsLegacyEventSplitFixture(testCase)
            %TESTLOADDATDETECTSLEGACYEVENTSPLITFIXTURE Regression check against the
            % preserved Astrocyte-produced fixture from
            % task-legacy-compat-fixture-testing.

            fixtureFolder = fullfile('D:', 'UMIT-DEV', 'LegacyCompatFixtures', ...
                'SaveFolder_Astrocyte_EventsSplit_20260819_074413');
            testCase.assumeTrue(isfolder(fixtureFolder), ...
                'Legacy compat fixture is not available on this machine.');

            eventsFile = fullfile(fixtureFolder, 'DATABYEVENTS.dat');
            testCase.verifyError(@() loadData(eventsFile), ...
                'Umitoolbox:loadMetaData:legacyEventSplitUnsupported');

            % The same fixture's plain continuous channel must still open normally.
            [~, Info] = loadData(fullfile(fixtureFolder, 'green.dat'));
            testCase.verifyEqual(Info.dimNames, {'Y','X','T'});
        end

        function testLoadDATCurrentFourDimSplitStillLoads(testCase)
            %TESTLOADDATCURRENTFOURDIMSPLITSTILLLOADS dev's own {Y,X,T,E} split
            % convention (datSize spans all 4 dims) must be unaffected by the new
            % legacy-split detection.

            data = testCase.BaseData(:,:,1:end-1);
            [Ny, Nx, Nt] = size(data);
            nEvents = 2;
            eventData = repmat(data, [1, 1, 1, nEvents]);

            datFile = fullfile(testCase.TempFolder, 'current_events.dat');
            fid = fopen(datFile, 'w');
            fwrite(fid, eventData, 'single');
            fclose(fid);

            % Convention 1: datSize spans all 4 dims, no separate datLength, so
            % Height/Width/T are resolved directly from datSize + dim_names.
            dim_names = {'Y','X','T','E'}; %#ok<NASGU>
            datSize = [Ny, Nx, Nt, nEvents]; %#ok<NASGU>
            Freq = 10; %#ok<NASGU>
            Datatype = 'single'; %#ok<NASGU>
            save(fullfile(testCase.TempFolder, 'current_events.mat'), ...
                'dim_names', 'datSize', 'Freq', 'Datatype');

            [~, Info] = loadData(datFile);
            testCase.verifyEqual(Info.dimNames, {'Y','X','T','E'});
            testCase.verifyEqual(Info.dimSizes, [Ny, Nx, Nt, nEvents]);
        end

        function testLoadDATRejectsUnrelatedMalformedDatSize(testCase)
            %TESTLOADDATREJECTSUNRELATEDMALFORMEDDATSIZE A non-2-and-2, non-YXT-only
            % mismatch must still fall through to the original generic error, so
            % the new legacy-split check doesn't over-match unrelated corruption.

            datFile = fullfile(testCase.TempFolder, 'malformed.dat');
            fid = fopen(datFile, 'w');
            fwrite(fid, zeros(1, 10, 'single'), 'single');
            fclose(fid);

            % 4 dim_names incl. non-YXT 'E', but datSize/datLength split 1-and-3
            % (not the legacy 2-and-2 pattern).
            dim_names = {'E','Y','X','T'}; %#ok<NASGU>
            datSize = 3; %#ok<NASGU>
            datLength = [10, 20, 5]; %#ok<NASGU>
            Freq = 10; %#ok<NASGU>
            Datatype = 'single'; %#ok<NASGU>
            save(fullfile(testCase.TempFolder, 'malformed.mat'), ...
                'dim_names', 'datSize', 'datLength', 'Freq', 'Datatype');

            testCase.verifyError(@() loadData(datFile), ...
                'Umitoolbox:loadMetaData:invalidDatSize');
        end

        function testSaveDATAcceptsOtherXY(testCase)
            %TESTSAVEDATACCEPTSOTHERXY A headered .dat may have any Y and X.
            %
            %   Before the .dat header, saveData rejected frames whose size
            %   differed from AcqInfos.mat; the header now describes the file.

            data = testCase.BaseData;
            [Ny, Nx, ~] = size(data);

            AcqInfoStream = TestIOHelpers.makeAcqInfoStreamForData(data, 10);
            AcqInfoStream.Height = Ny + 1;
            AcqInfoStream.Width = Nx;
            save(fullfile(testCase.TempFolder, 'AcqInfos.mat'), 'AcqInfoStream');

            filePath = fullfile(testCase.TempFolder, 'rawdata_other');
            saveData(filePath, data, 'DimNames', {'Y','X','T'}, 'FrameRateHz', 10);

            [outData, Info] = loadData(fullfile(testCase.TempFolder, 'rawdata_other.dat'));
            testCase.verifyEqual(outData, data);
            testCase.verifyEqual(Info.dimSizes, size(data));
        end

        function testLoadUMTRejectsMissingUMTContainer(testCase)
            %TESTLOADUMTREJECTSMISSINGUMTCONTAINER Reject files without a UMT container.

            badFile = fullfile(testCase.TempFolder, 'bad.umt');
            badStruct = struct('foo', 1);
            save(badFile, '-struct', 'badStruct', '-mat');

            testCase.verifyError(@() loadData(badFile), ...
                'Umitoolbox:loadMetaData:invalidUMT');
        end

        function testLoadUMTRejectsInvalidUMTStruct(testCase)
            %TESTLOADUMTREJECTSINVALIDUMTSTRUCT Reject malformed scalar UMT structures.

            badFile = fullfile(testCase.TempFolder, 'bad_scalar.umt');
            badUMT = struct();
            badUMT.version = 1;
            badUMT.kind = 'image';
            badUMT.data = struct();
            badUMT.data.main = struct('dimNames', {{'Y','X'}}); % Missing value.
            save(badFile, '-struct', 'badUMT', '-mat');

            testCase.verifyError(@() loadData(badFile), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end

        function testSaveUMTRejectsInvalidStruct(testCase)
            %TESTSAVEUMTREJECTSINVALIDSTRUCT Reject invalid struct when saving.

            badStruct = struct('foo', 1);
            badFile = fullfile(testCase.TempFolder, 'bad_save');

            testCase.verifyError(@() saveData(badFile, badStruct), ...
                'Umitoolbox:validateUMTStruct:invalidInput');
        end
    end
end
