%% Main script for applying a Pipeline to Labeo Datasets
% Add necessary folders to Matlab Path
toolboxFolder = 'D:\Academico\PostDoc_UdeM\LabeoTech\IsaToolbox';
IOI_ANAfolder = 'D:\Academico\PostDoc_UdeM\LabeoTech\ioi_ana';
addpath(genpath(toolboxFolder));
addpath(genpath(IOI_ANAfolder));
addpath(genpath('G:\RestingState_and_EyeTracking_SampleData'))
clearvars
%% Create Protocol Object
maindir = 'G:\RestingState_and_EyeTracking_SampleData';
savedir = 'G:\RS_data';
protocol = Protocol('TestProtocol', maindir, savedir, @protoFunc_RS, []);
protocol.generateList;
protocol.generateSaveFolders;
save([protocol.SaveDir protocol.Name '.mat'], 'protocol');
%% Query filter
% Clear previously saved Filter structure:
protocol.clearFilterStruct
% Query subjects
protocol.FilterStruct.Subject.PropName = '';
protocol.FilterStruct.Subject.Expression = '';
protocol.FilterStruct.Subject.LogicalOperator ='NOT';
    % Excludes a subject from query:
protocol.FilterStruct.Subject(2).PropName = 'ID';
protocol.FilterStruct.Subject(2).Expression = 'M0005844209';
% Query Acquisition and Modality:
protocol.FilterStruct.Acquisition.PropName = 'ID';
protocol.FilterStruct.Acquisition.Expression = 'RS_o'; % Leave empty to select all
protocol.FilterStruct.Modality.PropName = ''; % Leave empty to select all. This works as well.
% protocol.FilterStruct.Modality.Expression = '';

% Choose query method:
protocol.FilterStruct.FilterMethod = 'contains'; % Options: 'contains', 'regexp', 'strcmp';
% Perform query:
protocol.queryFilter;
% Display indices of selected branches from the Protocol hierarchy:
idx = protocol.Idx_Filtered

%% Preprocessing Pipeline
% Create pipeline
pipe = PipelineManager([], protocol);
% Delete intermediate files (optional). It will keep the first and last files from the Pipeline:
pipe.EraseIntermediate = true;
% TESTING:
pipe.addTask('FluorescenceImaging', 'dummyMethodForTesting')
pipe.addTask('FluorescenceImaging', 'dummyMethodForTesting_WithError')
pipe.addTask('FluorescenceImaging', 'dummyMethodForTesting3')

% Select Optional Parameters:
opts = pipe.setOpts('FluorescenceImaging', 'run_ImagesClassification');
pipe.addTask('FluorescenceImaging', 'run_ImagesClassification', 'opts', opts);
% pipe.addTask('FluorescenceImaging', 'run_ImagesClassification', 'output', 'Fluo_channel_file', 'opts', opts);
pipe.addTask('FluorescenceImaging', 'run_GSR');
% pipe.addTask('FluorescenceImaging', 'run_TemporalFilter');
pipe.addTask('FluorescenceImaging', 'run_SeedPixelCorrelation');

pipe.run_pipeline
% Save Pipeline
pipe.savePipe('testPipeline')
% Save Protocol
save(fullfile(protocol.SaveDir,[protocol.Name '.mat']), 'protocol')
%% Load Pipeline and run
pipe = PipelineManager([], protocol);
pipe.loadPipe('testPipeline.mat')
pipe.run_pipeline

%% Data Visualization
% SeedPixelCorrelation Visualization:
tmpM = protocol.Array.ObjList(2).Array.ObjList(1).Array.ObjList(1);
tmpM.viz_SeedPxCorr;






