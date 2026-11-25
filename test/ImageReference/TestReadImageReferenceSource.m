classdef (SharedTestFixtures = {UMITTemporaryRootFixture}) TestReadImageReferenceSource < matlab.unittest.TestCase
    %TESTREADIMAGEREFERENCESOURCE MAT-file sources of readImageReferenceSource.
    %
    %   A reference exported from the Image Reference Manager is a MAT file with one
    %   struct variable, ImageReference. "Create from file" used to refuse it ("no 2D
    %   numeric array variables"); it now reads ImageReference.image. Plain 2D numeric
    %   MAT files keep working. Images are non-square (6 x 8).

    properties
        Folder = ''
    end

    methods (TestMethodSetup)
        function createFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fixture.Folder;
        end
    end

    methods (Test)
        function anExportedImageReferenceFileIsReadAsItsImage(testCase)
            image = single(reshape(1:48, [6 8]));
            reference = genImageReferenceStruct(image, 'Name', 'Exported reference');
            file = exportImageReferenceFile(reference, fullfile(testCase.Folder, 'ref.mat'));

            [out, info] = readImageReferenceSource(file, 'Interactive', false);

            testCase.verifyEqual(out, image);
            testCase.verifyEqual(info.imageSizeYX, [6 8]);
            testCase.verifyEqual(info.selectedVariable, 'ImageReference');
            testCase.verifySubstring(info.processingNotes, 'exported Image Reference');
        end

        function aPlainTwoDimensionalMatFileStillWorks(testCase)
            plain = single(reshape(1:48, [6 8]));
            file = fullfile(testCase.Folder, 'plain.mat');
            save(file, 'plain');

            [out, info] = readImageReferenceSource(file, 'Interactive', false);

            testCase.verifyEqual(out, single(reshape(1:48, [6 8])));
            testCase.verifyEqual(info.selectedVariable, 'plain');
        end

        function aStructWithoutAnImageIsStillRefused(testCase)
            ImageReference = struct('name', 'no image here');
            file = fullfile(testCase.Folder, 'noimage.mat');
            save(file, 'ImageReference');

            testCase.verifyError(@() readImageReferenceSource(file, 'Interactive', false), ...
                'readImageReferenceSource:No2DNumericVariables');
        end

        function aNamedOtherVariableOverridesTheExportedReference(testCase)
            reference = genImageReferenceStruct(single(ones(6, 8)), 'Name', 'ref');
            other = single(reshape(1:48, [6 8]));
            ImageReference = reference;
            file = fullfile(testCase.Folder, 'both.mat');
            save(file, 'ImageReference', 'other');

            out = readImageReferenceSource(file, 'Variable', 'other', 'Interactive', false);

            testCase.verifyEqual(out, other);
        end
    end
end
