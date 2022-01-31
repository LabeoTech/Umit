function out = genDataHistory(fcnInfo, fcnStr, optsStruct, outFileList)
% This function creates a structure containing information about an
% analysis function. This structure will be added to the metaData of the
% data analysed. 
% This function is called by the "updateHistory" methods in
% "PipelineManager" and "DataViewer_pipelineMngr" classes, as well as by
% "manual_alignFrames" function.

% Inputs:
%   fcnInfo (struct): structure containing the function's basic informations with
%       fields:
%           -name (char): name of the analysis function.
%           -folder (char): path where the analysis function file is located.
%           -creationDatetime(datetime): timestamp of the creation of the
%               analysis function file.
%   fcnStr (char): string of the analysis function called by an EVAL function.
%   optsStruct (struct): structure containing optional parameters of the
%       analysis function.
%   outFileList (cell) : list of fileNames created by the analysis
%       function.
% Output:
%   out (struct): structure with the information necessary for the
%       dataHistory variable in the data's metaData.

out = struct('runDatetime', datetime('now'), 'name', {fcnInfo.name},...
                'folder', {fcnInfo.folder}, 'creationDatetime', datetime(fcnInfo.date),...
                'opts', optsStruct, 'funcStr', {fcnStr}, 'outputFile_list', outFileList);
end