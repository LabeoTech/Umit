classdef TestEventsManagerPlot < matlab.unittest.TestCase
    %TESTEVENTSMANAGERPLOT Smoke tests for plotting behavior.

    properties
        TempFolder char
        OldFigureVisible
    end

    methods (TestMethodSetup)
        function createTempFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fx = testCase.applyFixture(TemporaryFolderFixture);
            testCase.TempFolder = fx.Folder;
            addpath(fullfile(fileparts(mfilename('fullpath')), 'EMTestHelpers'));

            testCase.OldFigureVisible = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            close all force
        end
    end

    methods (TestMethodTeardown)
        function closeFigures(testCase)
            close all force
            set(groot, 'DefaultFigureVisible', testCase.OldFigureVisible);
        end
    end

    methods (Test)
        function testPlotAnalogSmoke(testCase)
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000 15000], ...
                'PulseWidth', 80, ...
                'StimName', 'Main');

            obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            ax = obj.plot();

            % plot() returns the per-channel axes, not the parent figure --
            % reach the figure through the axes to check its Tag.
            testCase.verifyTrue(all(ishghandle(ax)));
            f = ancestor(ax(1), 'figure');
            testCase.verifyEqual(get(f, 'Tag'), 'EventsManager_AnalogINPlot');
            testCase.verifyNumElements(ax, numel(obj.AIChanList));
        end

        function testPlotReusesSingleTaggedFigure(testCase)
            buildAnalogAcquisitionFolder(testCase.TempFolder, ...
                'PulseStarts', [5000 10000], ...
                'PulseWidth', 80, ...
                'StimName', 'Main');

            obj = EventsManager(testCase.TempFolder, testCase.TempFolder, 'csv');
            ax1 = obj.plot();
            testCase.verifyTrue(all(ishghandle(ax1)));
            ax2 = obj.plot();
            testCase.verifyTrue(all(ishghandle(ax2)));

            figs = findall(groot, 'Type', 'figure', 'Tag', 'EventsManager_AnalogINPlot');
            testCase.verifyNumElements(figs, 1);
            testCase.verifyEqual(figs, ancestor(ax2(1), 'figure'));
        end

        function testPlotExternalSignalSmoke(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.TempFolder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000 11000], 60, 'Amplitude', 5);
            obj.getTriggersFromSignal(signal, 10000, false);

            ax = obj.plot();

            testCase.verifyTrue(all(ishghandle(ax)));
            testCase.verifyNumElements(ax, 1);
            ttl = get(get(ax, 'Title'), 'String');
            testCase.verifyEqual(ttl, 'extSignal');
        end

        function testPlotExternalSignalIgnoresChannelSelection(testCase)
            writeMinimalAcqInfosMat(testCase.TempFolder, 'FrameRateHz', 60, 'AISampleRate', 10000);
            obj = EventsManager(testCase.TempFolder, '', 'csv');
            signal = makePulseSignal(22000, [3000 7000], 60, 'Amplitude', 5);
            obj.getTriggersFromSignal(signal, 10000, false);

            warnState = warning('off', 'all');
            cleanupObj = onCleanup(@() warning(warnState));
            lastwarn('');
            ax = obj.plot({'StimAna1'});
            [warnMsg, ~] = lastwarn;

            testCase.verifyTrue(all(ishghandle(ax)));
            testCase.verifyNotEmpty(warnMsg);
        end
    end
end
