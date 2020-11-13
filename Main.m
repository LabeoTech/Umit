%% Main script for applying a Pipeline to Labeo Datasets
% Add necessary folders to Matlab Path
addpath(genpath('D:\Academico\PostDoc_UdeM\LabeoTech\IsaToolbox'));
addpath(genpath('D:\Academico\PostDoc_UdeM\LabeoTech\ioi_ana'));
addpath(genpath('D:\Academico\PostDoc_UdeM\LabeoTech\OOP_scripts'));
addpath(genpath('G:\SampleTestFolder'));
clearvars

%% Create Protocol Object
maindir = 'G:\RetinoSampleData';
myProtocol = Protocol(maindir, maindir, @ProtoFunc4, []);
myProtocol.generateList;

%% Query filter
query_filter = myProtocol.createFilterStruct;
query_filter.Subject.PropName = 'ID';
query_filter.Subject.Expression = 'Mouse_\d*';
query_filter.Acquisition.PropName = 'ID';
query_filter.Acquisition.Expression = 'Acq_\d*';
query_filter.Modality.PropName = 'ID';
query_filter.Modality.Expression = 'Labeo';

%% Preprocessing Pipeline
ppPipe.a.ClassName = 'Labeo';
ppPipe.a.FuncName = 'ImagesClassification';
ppPipe.a.FuncParams = 'FOLDER, 1, 1, 0, 0';

ppPipe.b.ClassName = 'Labeo';
ppPipe.b.FuncName = 'Ana_IOI_FullFrame';
ppPipe.b.FuncParams = 'FOLDER, 0, 1, 1, []';

%% Run PreProcessing Pipeline
myProtocol.runPreProcessingPipeline(ppPipe,[])






test = myProtocol.Array.ObjList.Array.ObjList(1);
test.Array.ObjList(1).run_Ana_IOI_FullFrame(test.Array.ObjList(1).Folder, 0, 1, 1, [])
log = test.Array.ObjList(1).LastLog;
test.save2LogBook(log);










