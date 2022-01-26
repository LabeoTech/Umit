%% Main script for applying a Pipeline to Labeo Datasets
% Add necessary folders to Matlab Path
clearvars
toolboxFolder = 'D:\Academico\PostDoc_UdeM\LabeoTech\IsaToolbox';
IOI_ANAfolder = 'D:\Academico\PostDoc_UdeM\LabeoTech\ioi_ana';
addpath(genpath(toolboxFolder));
addpath(genpath(IOI_ANAfolder));
%% Create Protocol Object
maindir = 'F:\SolennRawData';
savedir = 'E:\TestSaveDir';
protocol = Protocol('TestProtocol2', maindir, savedir, @TestDataProtoFunc, []);
protocol.generateList;
protocol.generateSaveFolders;
save(fullfile(protocol.SaveDir, [protocol.Name '.mat']), 'protocol');
%% Load existing Protocol Object:
addpath(genpath('G:\DummyDataSet_4_testing'));
maindir = 'F:\Solenn';
savedir = 'E:\Solenn_Data\ProjectSaveDir';
ProjectName = 'Solenn_backupData_2.mat';
load(fullfile(savedir, ProjectName));
%% Query filter
% Clear previously saved Filter structure:
protocol.clearFilterStruct
% Query subjects
protocol.FilterStruct.Subject.PropName = 'ID';
protocol.FilterStruct.Subject.Expression = 'M5';
% protocol.FilterStruct.Subject.LogicalOperator ='NOT';
    % Excludes a subject from query:
% protocol.FilterStruct.Subject(2).PropName = 'ID';
% protocol.FilterStruct.Subject(2).Expression = 'M0005844209';
% Query Acquisition and Modality:
protocol.FilterStruct.Acquisition.PropName = 'ID';
protocol.FilterStruct.Acquisition.Expression = 'OD_SFTF_21ST18'; % Leave empty to select all
% protocol.FilterStruct.Acquisition.LogicalOperator ='OR';
% protocol.FilterStruct.Acquisition(2).PropName = 'ID';
% protocol.FilterStruct.Acquisition(2).Expression = 'RS'; % Leave empty to select all

% protocol.FilterStruct.Modality.PropName = 'ID'; % Leave empty to select all. This works as well.
% protocol.FilterStruct.Modality.Expression = 'Ctx';

% Choose query method:
protocol.FilterStruct.FilterMethod = 'contains'; % Options: 'contains', 'regexp', 'strcmp';
% Perform query:
protocol.queryFilter;
% Display indices of selected branches from the Protocol hierarchy:
idx = protocol.Idx_Filtered

%% Preprocessing Pipeline
% Create pipeline
pipe = PipelineManager(protocol, 'FluorescenceImaging');
% Show list of Available analysis functions:
pipe.showFuncList;
% Set Optional Parameters for "run_ImagesClassification" function:
pipe.setOpts(15)
% Add "run_ImagesClassification" to the pipeline:
% Example of adding tasks to pipeline using functions indices:
pipe.addTask(3);
pipe.addTask(7);
pipe.addTask(4);
% Example of adding tasks and saving outputs:
pipe.addTask(1,true, 'testout');
pipe.addTask(7, true, 'tempOut');
% Example of pipeline construction with function names as input:
pipe.addTask('alignFrames');
pipe.addTask('calculateDF_F0');
pipe.addTask('getEventsFromSingleChannel')
pipe.addTask('SeedPixCorr')
%% Overview of pipeline
pipe.showPipeSummary
%% Run pipeline
pipe.run_pipeline
%% Save Pipeline
pipe.savePipe('testPipeline')
%% Load Pipeline and run
pipe = PipelineManager(protocol, 'FluorescenceImaging');
pipe.loadPipe('TestPipe1')
pipe.showPipeSummary
pipe.run_pipeline

