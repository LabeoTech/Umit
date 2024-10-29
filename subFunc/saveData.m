function saveData(filename,SaveFolder,data,varargin)
% SAVE2DAT saves data to binary .dat files.
% Inputs:
%   filename (str): name of the file to save the data.
%   SaveFolder (str): path to folder where the file will be saved.
%   data (numerical array | struct): non-empty "single" numeric 3D matrix
%   OR data packaged in a structure (see genDataStructure.m, for details, basic explanation below).
%
%   (Optional): !Options available only for numerical data!
%       AcqInfoStream(struct): structure with acquisition information. See "AcqInfos.mat"
%       file.
%       Append(bool): If true, appends new data to existing .dat file.
%
% Note #1: If the input data is a 3D numerical array, it will be saved with a ".dat" file extension. 
%          It must have the same dimensions stated in the "AcqInfos.mat" file in the save folder.
% Note #2: If there is no "AcqInfos.mat" file in the saveFolder, it is necessary to
%          provide the "AcqInfoStream" structure!
% Note #3: If the input data is a structure, it will be saved with a
%          ".datstat" file extension. This structure must contain the
%          following fields to be considered valid:
%               data (struct): structure containing the data.
%               obsID (cell): cell array with the same length as data
%                             containing the observation IDs.
%               b_hasEvents (bool): Boolean flag indicating if the data is
%                                   split by events.
%               b_hasMultipleMeasures(bool): Boolean flag indicating if th
%                                           data structure contain multiple measures (structure with
%                                           more than a single field).
%               dataCategory (struct): structure containing the data
%                                      category for each field of the "data" structure.
%                                      It must be one of the following:
%                           1) scalar
%                           2) time-vector
%                           3) map (2D) 
%                           4) image-time-series (3D)
%
%

% Arguments validation
p = inputParser;
addRequired(p,'filename', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p,'SaveFolder',@isfolder);
addRequired(p, 'data', @(x) ( isnumeric(x) & strcmpi(class(x),'single') ) | isstruct(x));
addOptional(p, 'AcqInfoStream', struct.empty(0,1),@isstruct);
addParameter(p,'Append',false,@islogical);
parse(p, filename, SaveFolder, data,varargin{:});
%%%%%%
filename = convertStringsToChars(filename);
SaveFolder = convertStringsToChars(SaveFolder);
[~,file,~] = fileparts(filename);
if isstruct(data)
    % Save to .mat file:
    filename = [file, '.datstat']; % Add extension for data in structure
    save2mat(fullfile(SaveFolder,filename),data);
else
     % Save to .dat file:   
     filename = [file '.dat']; % Add extension for data as numeric array.
     save2dat(filename,SaveFolder,data,p.Results.AcqInfoStream, p.Results.Append);
end
disp(['Data saved as "' filename '" in folder "' SaveFolder '"']);

end
% Local functions:

function save2dat(file, saveFolder,data,AcqInfoStream, b_append)
% Saves numeric data array to a .dat file.
% Check if there is an "AcqInfos.mat" file in the SaveFolder:

% 
if ~isfile(fullfile(saveFolder,'AcqInfos.mat'))    
    if isempty(AcqInfoStream)
        % Raise error. An "AcqInfos.mat" is necessary!
        error(['Failed to save "' file  '"! No "AcqInfos.mat" file found in folder "' saveFolder '"!']);
    else
        % Save provided AcqInfoStream
        save(fullfile(saveFolder,'AcqInfos.mat'),'AcqInfoStream');
    end
end

% Check if the data has the same dimensions as in AcqInfoStream:
load(fullfile(saveFolder,'AcqInfos.mat'));%#ok
assert(isequaln(size(data),[AcqInfoStream.Width, AcqInfoStream.Height, AcqInfoStream.Length]),...
    'Umitoolbox:saveData:invalidInput',...
    'Operation aborted! Data does not have the same dimensions saved in the "AcqInfos.mat" file!');
clear AcqInfoStream

% Set permission type:
permission = 'w';
if b_append
    % Set data writing permission to "append":
    permission = 'a';
    % Update Data Length in AcqInfos.mat file:
    load(fullfile(saveFolder,'AcqInfos.mat'));%#ok      
    AcqInfoStream.Length = AcqInfoStream.Length + size(data,3);
    save(fullfile(saveFolder,'AcqInfos.mat'),'AcqInfoStream');
    clear AcqInfoStream;
end
% Save the data to a .dat file:
disp('Writing data to .DAT file ...');
fid = fopen(fullfile(saveFolder,file), permission);
fwrite(fid, data, 'single'); % Write data as single precision.
fclose(fid);
end

function save2mat(filename,data)
% Saves data structure to a .datstat file.

% Validate data structure before saving:
fn = fieldnames(data);
errID = 'Umitoolbox:saveData:invalidInput';
errMsg = 'Operation aborted! Invalid data structure! See genDataStructure.m documentation for details.';
assert(all(ismember({'data','obsID','b_hasEvents','b_hasMultipleMeasures','dataCategory'},...
    fn)),errID, errMsg)
assert(isstruct(data.data) & ...
    iscell(data.obsID) & ...
    islogical(data.b_hasEvents) & ...
    islogical(data.b_hasMultipleMeasures) & ...
    isstruct(data.dataCategory), errID, errMsg);
% Save data to "datstat" file:
disp('Writing data to .DATSTAT file ...');
save(filename,'-struct','data','-mat')
end

