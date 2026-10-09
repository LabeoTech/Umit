classdef TestGenImageTimeSeriesUMT < matlab.unittest.TestCase
    methods (Test)
        function testPipelineInfo(testCase)
            info = genImageTimeSeriesUMT('pipelineInfo'); %#ok<NASGU>
            testCase.verifyEqual(info.name, 'genImageTimeSeriesUMT');
            testCase.verifyEqual(info.outputs(1).type, {'ProcessedData'});
        end

        function testPackagesYXTData(testCase)
            data = rand(8, 7, 11, 'single');
            out = genImageTimeSeriesUMT(data);

            testCase.verifyTrue(isstruct(out) && isscalar(out));
            testCase.verifyEqual(out.kind, 'image');
            testCase.verifyTrue(isfield(out.data, 'main'));
            testCase.verifyEqual(cellstr(string(out.data.main.dimNames)), {'Y','X','T'});
            testCase.verifyEqual(out.data.main.value, data);
        end

        function testCustomEntryName(testCase)
            data = rand(5, 6, 9, 'single');
            out = genImageTimeSeriesUMT(data, 'EntryName', 'fluo');
            testCase.verifyTrue(isfield(out.data, 'fluo'));
            testCase.verifyFalse(isfield(out.data, 'main'));
        end

        function testOptionalLabels(testCase)
            data = rand(4, 5, 6, 'single');
            labels = struct();
            labels.T = {'t1','t2','t3','t4','t5','t6'};

            out = genImageTimeSeriesUMT(data, 'Labels', labels);
            testCase.verifyTrue(isfield(out, 'labels'));
            testCase.verifyEqual(out.labels.T, labels.T);
        end

        function testRejectsNon3DInput(testCase)
            data = rand(5, 6, 'single');
            testCase.verifyError(@() genImageTimeSeriesUMT(data), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
        end
    end
end
