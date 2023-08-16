%% Main script for applying a Pipeline to Labeo Datasets
% Add necessary folders to Matlab Path
clearvars
toolboxFolder = 'C:\PATH_TO_TOOLBOX_FOLDER';
addpath(genpath(toolboxFolder));
%% Set Directories Paths and Project Name:
maindir = 'C:\PATH_TO_RAW_FOLDER';
savedir = 'C:\PATH_TO_SAVE_FOLDER';
ProjectName = 'TYPEprojectNAME';

%% Create Protocol Object
protocol = Protocol(ProjectName, maindir, savedir, @protocolFcn_template, []);
protocol.generateList;
protocol.generateSaveFolders;
save(fullfile(protocol.SaveDir, [protocol.Name '.mat']), 'protocol');
%% Load existing Protocol Object:
load(fullfile(savedir, ProjectName));
%% Query filter
% Clear previously saved Filter structure:
protocol.clearFilterStruct
% Query subjects
% protocol.FilterStruct.Subject.PropName = 'ID';
% protocol.FilterStruct.Subject.Expression = 'M5';
% protocol.FilterStruct.Subject.LogicalOperator ='NOT';
    % Excludes a subject from query:
% protocol.FilterStruct.Subject(2).PropName = 'ID';
% protocol.FilterStruct.Subject(2).Expression = 'M0005844209';
% Query Acquisition and Modality:
protocol.FilterStruct.Acquisition.PropName = 'ID';
protocol.FilterStruct.Acquisition.Expression = 'RS'; % Leave empty to select all
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
idx = protocol.Idx_Filtered;

%% Preprocessing Pipeline
% Create pipeline
pipe = PipelineManager(protocol, 'FluorescenceImaging');
% Show list of Available analysis functions:
pipe.showFuncList;
% Set Optional Parameters for "run_ImagesClassification" function:
pipe.setOpts('run_ImagesClassification')
% Add "run_ImagesClassification" to the pipeline:
% Example of adding tasks to pipeline using functions names:
pipe.addTask('run_ImagesClassification');
pipe.addTask('GSR');
pipe.addTask('spatialGaussFilt');
% Example of adding tasks and saving outputs:
pipe.addTask('GSR',true, 'testout');
pipe.addTask('alignFrames', true, 'tempOut');
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

%% Data analysis
% In this section, we put together the data created using one of the
% "Analysis" functions (e.g. transformDat2Mat, genCorrelationMatrix). 
% These files contain the data from observations (i.e. ROIs). Here, we will
% use the "StatsManager" class to group the data into one or more
% experimental groups and gather the data for plotting.

% Select .MAT file containing the data:
fileName = 'corrMatrix_avg_vs_avg.mat';
% Get the list of objects containing the data
list_of_objs = protocol.extractFilteredObjects(3);
% Get the list of observations contained in all files:
obs_list = {'A_R','AL_R','AM_R','V1_L', 'M1_L', 'M2_L'};
% Set a list of experimental groups for each item in the "list_of_obj":
list_of_groups = repmat({'Test'},size(list_of_objs));
% Instantiate the "StatsManager" object:
statMngr = StatsManager(list_of_objs, obs_list,list_of_groups, fileName);
%% Plot grouped data
PlotLongData(statMngr); % Plotting tool for scalar and time-series data types.
%% Plot Correlation Matrices
PlotCorrMatrix(statMngr) % Plotting tool for correlation matrices.
%% Perform statistical analysis for scalar data
% Average acquisitions
statMngr.getAcqIndexList('original')
statMngr.resetAvgIndex; % Reset to original Acquisitions 
% Merge all
statMngr.setAcquisitionRange([1 4])
% Merge consecutive groups
statMngr.setAcquisitionRange([1 2])
statMngr.setAcquisitionRange([3 4])
% Set statistical independent and grouping variables:
statMngr.setStatsVariables({'ROI','Group'},'Acquisition')
% Run stats:
[a,b,c] = statMngr.runStatsOnScalar;
% Generate stats report:
statMngr.exportReportToTXT('reportMatrix.txt')
%% Perform statistical analysis for correlation matrices
[qMatrix, statsReport] = statMngr.runStatsOnMatrix;
% Plot matrix:
figure; 
s1=subplot(121); 
imagesc(s1,qMatrix,[0,1]);
title(s1,'qValues');
colormap(s1,'jet'); axis square;
colorbar(s1,'eastoutside')
s2=subplot(122);
imagesc(s2,qMatrix<=statMngr.pAlpha);
title(s2,['qValues \leq' num2str(statMngr.pAlpha)]); 
colormap(s2,'gray'); axis square; colorbar(s2,'eastoutside')
sgtitle('FDR-corrected p Values (q Values)'); 







