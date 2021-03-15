%% Main script for applying a Pipeline to Labeo Datasets
% Add necessary folders to Matlab Path
clearvars
toolboxFolder = 'D:\Academico\PostDoc_UdeM\LabeoTech\IsaToolbox';
IOI_ANAfolder = 'D:\Academico\PostDoc_UdeM\LabeoTech\ioi_ana';
addpath(genpath(toolboxFolder));
addpath(genpath(IOI_ANAfolder));
%% Create Protocol Object
maindir = 'G:\Project_RS_EyeAndBodyTracking';
addpath(genpath(maindir))
savedir = 'G:\Testing_new';
protocol = Protocol('TestProtocol', maindir, savedir, @protocolFunction_RS, []);
protocol.generateList;
protocol.generateSaveFolders;
save([protocol.SaveDir protocol.Name '.mat'], 'protocol');
%% Query filter
% Clear previously saved Filter structure:
protocol.clearFilterStruct
% Query subjects
protocol.FilterStruct.Subject.PropName = 'ID';
protocol.FilterStruct.Subject.Expression = 'M0004855451';
% protocol.FilterStruct.Subject.LogicalOperator ='NOT';
    % Excludes a subject from query:
% protocol.FilterStruct.Subject(2).PropName = 'ID';
% protocol.FilterStruct.Subject(2).Expression = 'M0005844209';
% Query Acquisition and Modality:
% protocol.FilterStruct.Acquisition.PropName = 'ID';
% protocol.FilterStruct.Acquisition.Expression = 'RS_o'; % Leave empty to select all
protocol.FilterStruct.Modality.PropName = 'ID'; % Leave empty to select all. This works as well.
protocol.FilterStruct.Modality.Expression = 'Ctx';

% Choose query method:
protocol.FilterStruct.FilterMethod = 'contains'; % Options: 'contains', 'regexp', 'strcmp';
% Perform query:
protocol.queryFilter;
% Display indices of selected branches from the Protocol hierarchy:
idx = protocol.Idx_Filtered;

%% Preprocessing Pipeline
% Create pipeline
pipe = PipelineManager([], protocol, 'D:\Academico\PostDoc_UdeM\LabeoTech\IsaToolbox\Analysis');
% Delete intermediate files (optional). It will keep the first and last files from the Pipeline:
pipe.EraseIntermediate = true;
pipe.IgnoreLoggedFiles = true;
% Select Optional Parameters:
opts = pipe.setOpts('run_ImagesClassification');
pipe.addTask('FluorescenceImaging', 'run_ImagesClassification', opts);
% pipe.addTask('FluorescenceImaging', 'run_ImagesClassification')
pipe.addTask('FluorescenceImaging', 'dummyFunc4Testing');
pipe.addTask('FluorescenceImaging', 'dummyFunc4Testing_2');
pipe.addTask('FluorescenceImaging', 'dummyFunc4Testing_3');

pipe.addTask('FluorescenceImaging', 'tempFiltNormalize');
pipe.addTask('FluorescenceImaging', 'GSR')
pipe.addTask('FluorescenceImaging', 'SeedPixCorr');
% 
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
% SeedPixelCorrelation Visualization:
tmpM = protocol.Array.ObjList(1).Array.ObjList(2).Array.ObjList(1);
tmpM.viz_SeedPxCorr;






