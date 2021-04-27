%% Main script for applying a Pipeline to Labeo Datasets
% Add necessary folders to Matlab Path
clearvars
toolboxFolder = 'D:\Academico\PostDoc_UdeM\LabeoTech\IsaToolbox';
IOI_ANAfolder = 'D:\Academico\PostDoc_UdeM\LabeoTech\ioi_ana';
addpath(genpath(toolboxFolder));
addpath(genpath(IOI_ANAfolder));
%% Create Protocol Object
addpath(genpath('E:\Solenn_Data'));
maindir = 'F:\Solenn';
savedir = 'E:\Solenn_Data\ProjectSaveDir';
protocol = Protocol('Solenn_backupData_2', maindir, savedir, @SolenDataProtoFunc, []);
protocol.generateList;
protocol.generateSaveFolders;
save(fullfile(protocol.SaveDir, [protocol.Name '.mat']), 'protocol');
%% Load existing Protocol Object:
addpath(genpath('G:\DummyDataSet_4_testing'));
maindir = 'G:\DummyDataSet_4_testing';
savedir = 'G:\DummyDataSet_4_testing\ProjectSaveDir';
ProjectName = 'dummyTest2';
load(fullfile(savedir, ProjectName));
%% Query filter
% Clear previously saved Filter structure:
protocol.clearFilterStruct
% Query subjects
protocol.FilterStruct.Subject.PropName = 'ID';
protocol.FilterStruct.Subject.Expression = 'M13';
% protocol.FilterStruct.Subject.LogicalOperator ='NOT';
    % Excludes a subject from query:
% protocol.FilterStruct.Subject(2).PropName = 'ID';
% protocol.FilterStruct.Subject(2).Expression = 'M0005844209';
% Query Acquisition and Modality:
protocol.FilterStruct.Acquisition.PropName = 'ID';
protocol.FilterStruct.Acquisition.Expression = 'SFTF'; % Leave empty to select all
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
pipe = PipelineManager([], protocol, 'D:\Academico\PostDoc_UdeM\LabeoTech\IsaToolbox\Analysis');
% Delete intermediate files (optional). It will keep the first and last files from the Pipeline:
% pipe.EraseIntermediate = true;
pipe.IgnoreLoggedFiles = true;
% Select Optional Parameters:
pipe.addTask('FluorescenceImaging', 'dummyFunc4Testing');
pipe.addTask('FluorescenceImaging', 'dummyFunc4Testing_3');
pipe.addTask('FluorescenceImaging', 'dummyFunc4Testing_4');
pipe.addTask('FluorescenceImaging', 'dummyFunc4Testing_5');

    
pipe.addTask('FluorescenceImaging', 'alignFrames');
pipe.addTask('FluorescenceImaging', 'calculateDF_F0');
pipe.addTask('FluorescenceImaging', 'getEventsFromSingleChannel')
pipe.addTask('FluorescenceImaging', 'SeedPixCorr')


opts = pipe.setOpts('run_ImagesClassification');
pipe.addTask('FluorescenceImaging', 'run_ImagesClassification', opts);
opts = pipe.setOpts('getEventsFromSingleChannel');
pipe.addTask('FluorescenceImaging', 'getEventsFromSingleChannel', opts);
opts = pipe.setOpts('alignFrames');
pipe.addTask('FluorescenceImaging', 'alignFrames', opts);
pipe.addTask('FluorescenceImaging', 'calculateDF_F0');
opts = pipe.setOpts('event_triggered_average');
pipe.addTask('FluorescenceImaging', 'event_triggered_average', opts);
% pipe.addTask('FluorescenceImaging', 'tempFiltNormalize');
% pipe.addTask('FluorescenceImaging', 'GSR')
pipe.addTask('FluorescenceImaging', 'MV_getCentroidValues');
pipe.addTask('FluorescenceImaging', 'MVcalculate_SF_TF_average');
pipe.addTask('FluorescenceImaging', 'SeedPixCorr');
pipe.addTask('FluorescenceImaging','lookForMetaDataFile');
%% Overview of pipeline
pipe.showPipeSummary
%% Run pipeline
pipe.run_pipeline
%% Save Pipeline
pipe.savePipe('testPipeline')
% Save Protocol
save(fullfile(protocol.SaveDir,[protocol.Name '.mat']), 'protocol')
%% Load Pipeline and run
pipe = PipelineManager([], protocol);
pipe.loadPipe('testPipeline.mat')
pipe.run_pipeline

%% Data Visualization

selec = protocol.extractFilteredObjects(3);
retino_maps(selec,'fChan.dat')





