function saveDat(DatFileName, data,varargin)
% SAVE2DAT saves data to binary .dat files.
% Inputs:
%   DatFileName (str): fullpath of .DAT file.
%   data (numerical array): non-empty "single" numeric multi-dimensional matrix.
%   (Optional:
%       AcqInfoStream(struct): structure with acquisition information. See "AcqInfos.mat"
%       file.
%       Append(bool):If true, appends new data to existing .dat file.
% Note:
% If there is no "AcqInfos.mat" file in the saveFolder, it is necessary to
% provide the "AcqInfoStream" structure!

% Arguments validation
p = inputParser;
addRequired(p,'DatFileName', @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(p, 'data', @(x) validateattributes(x, {'single'}, {'nonempty'}));
addOptional(p, 'AcqInfoStream', struct.empty(0,1),@isstruct);
addParameter(p,'Append',false,@islogical);
parse(p, DatFileName, data,varargin{:});
DatFileName = p.Results.DatFileName;
data = p.Results.data;
AcqInfoStream = p.Results.AcqInfoStream;
%%%%%%

% Check if there is an "AcqInfos.mat" file in the SaveFolder:
[saveFolder,file,ext] = fileparts(DatFileName);
if isempty(saveFolder)
    saveFolder = pwd;
end
% 
if ~isfile(fullfile(saveFolder,'AcqInfos.mat'))    
    if isempty(AcqInfoStream)
        % Raise error. An "AcqInfos.mat" is necessary!
        error(['Failed to save "' file ext '"! No "AcqInfos.mat" file found in folder "' saveFolder '"!']);
    else
        % Save provided AcqInfoStream
        save(fullfile(saveFolder,'AcqInfos.mat'),'AcqInfoStream');
    end
end
clear AcqInfoStream
% Force filename extension:
if ~endsWith(DatFileName, '.dat')
    DatFileName = fullfile(saveFolder,[file, '.dat']);
end
% Set permission type:
permission = 'w';
if p.Results.Append
    % Set data writing permission to "append":
    permission = 'a';
    % Update Data Length in AcqInfos.mat file:
    load(fullfile(saveFolder,'AcqInfos.mat'));      
    AcqInfoStream.Length = AcqInfoStream.Length + size(data,3);
    save(fullfile(saveFolder,'AcqInfos.mat'),'AcqInfoStream');
    clear AcqInfoStream;
end

disp('Writing data to .DAT file ...');
fid = fopen(DatFileName, permission);
fwrite(fid, data, 'single'); % Write data as single precision.
fclose(fid);
disp(['Data saved in : "' DatFileName '"']);
end