%% Main script for applying a Pipeline to Labeo Datasets
% Add necessary folders to Matlab Path
paths = fullfile('D:\Academico\PostDoc_UdeM\LabeoTech', {'IsaToolbox', 'ioi_ana', 'OOP_scripts'});
for i = 1:length(paths)
addpath(genpath(paths{i}));
end
clearvars
%% Create Protocol Object
maindir = 'D:\Academico\PostDoc_UdeM\LabeoTech\OOP_scripts\myDummyDataSet';
myProtocol = Protocol(maindir, maindir, @ProtoFunc4, []);
myProtocol.generateList;

%% Query filter
query_filter = myProtocol.createFilterStruct;
query_filter.Subject.PropName = 'ID';
query_filter.Subject.Expression = 'Mouse_00\d*';
query_filter.Acquisition.PropName = 'ID';
query_filter.Acquisition.Expression = 'Acq_00\d';
query_filter.Modality.PropName = 'ID';
query_filter.Modality.Expression = 'Labeo';

%% 


