classdef TestRetinotopyFFTWindow < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPath(testCase)
            thisFile = mfilename('fullpath');
            projectRoot = extractBefore(fileparts(thisFile), [filesep 'test']);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                projectRoot, 'IncludingSubfolders', true));
        end
    end

    methods (Test)
        function testNonAverageModeConcatenatesOnlyStimulusEpochs(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            folder = fixture.Folder;

            writeMinimalAcqInfosMat(folder, ...
                'FrameRateHz', 1, 'AISampleRate', 1000);

            eventID = uint16([1; 1; 2; 2; 1; 1; 2; 2]);
            state = logical([1; 0; 1; 0; 1; 0; 1; 0]);
            timestamps = single([1; 3; 5; 7; 9; 11; 13; 15]);
            eventNameList = {'0','180'};
            repetitionID = uint16([1; 1; 1; 1; 2; 2; 2; 2]);
            selectedEvents = true(size(eventID));
            baselinePeriod = single(1);
            save(fullfile(folder, 'events.mat'), ...
                'eventID', 'state', 'timestamps', 'eventNameList', ...
                'repetitionID', 'selectedEvents', 'baselinePeriod');

            data = 100 * ones(1, 1, 16, 'single');
            data(1,1,[1 2 9 10]) = single([1 -1 1 -1]);
            data(1,1,[5 6 13 14]) = single([-1 1 -1 1]);

            out = genRetinotopyMaps(data, folder, ...
                'FrameRateHz', 1, ...
                'Direction', 'Azimuth_only', ...
                'b_useAverageMovie', false);

            testCase.verifyEqual( ...
                out.data.AzimuthMap.value(1,1,1), single(2), ...
                'AbsTol', single(10 * eps('single')));
        end
    end
end
