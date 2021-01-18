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
protocol = Protocol('RSandET_project', maindir, savedir, @protoFunc_RS, []);
protocol.generateList;
protocol.generateSaveFolders;
save([savedir filesep 'ProtocolFile.mat'], 'protocol');
%% Query filter
protocol.FilterStruct.Subject.PropName = '';
protocol.FilterStruct.Subject.Expression = '';
protocol.FilterStruct.Acquisition.PropName = 'ID';
protocol.FilterStruct.Acquisition.Expression = 'RS_and_ET'; % Leave empty to select all
protocol.FilterStruct.Modality.PropName = ''; % Leave empty to select all. This works as well.
% protocol.FilterStruct.Modality.Expression = '';
protocol.FilterStruct.FilterMethod = 'strcmp'; % Options: 'contains', 'regexp', 'strcmp';
protocol.queryFilter;
idx = protocol.Idx_Filtered;
%% Preprocessing Pipeline
% Create pipeline
pipe = PipelineManager([], protocol);
% Select Optional Parameters:
opts = pipe.setOpts('CalciumImaging', 'run_ImagesClassification');
pipe.addTask('CalciumImaging', 'run_ImagesClassification', 'output', 'Fluo_channel_file', 'opts', opts);
pipe.addTask('CalciumImaging', 'run_GSR');
% pipe.addTask('CalciumImaging', 'run_TemporalFilter');
pipe.addTask('CalciumImaging', 'run_SeedPixelCorrelation');

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






