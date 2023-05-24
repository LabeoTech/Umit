function out = genDataHistory(fcnInfo, optsStruct, inputFileName, outFileName)
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
%   optsStruct (struct): structure containing optional parameters of the
%       analysis function.
%   inputFileName(cell|char): name of the input file(s) to the function. 
%   outFileName(cell|char): name of the output file(s) from the function.
%   This field is used just by functions that create files already. Just to
%   keep track of the files that were created.
% Output:
%   out (struct): structure with the information necessary for the
%       dataHistory variable in the data's metaData.

out = struct('runDatetime', datetime('now'), 'name', {fcnInfo.name},...
                'folder', {fcnInfo.folder}, 'creationDatetime',...
                datetime(fcnInfo.datenum, 'ConvertFrom', 'datenum'),...
                'opts', optsStruct, 'inputFileName',inputFileName, 'outFileName',outFileName);
end