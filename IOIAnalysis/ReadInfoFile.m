function out = ReadInfoFile(FolderPath, varargin)
% This function parses the "info.txt" file and saves the data to a
% structure.
% Inputs:
% FolderPath (char): Path to folder containing the "info.txt" file.
% infoFile (char): Optional Parameter. Name of the .TXT file containing the
% acquisition information. Use this parameter to read a file with a
% different name as "info".

% Read the info.txt file:
if nargin == 1
    txt = readcell(fullfile(FolderPath, 'info.txt'), 'Delimiter', ':', 'NumHeaderLines',1);
else
    [~,infoFile,ext] = fileparts(varargin{:});
    if isempty(ext)
        ext = '.txt';
    end        
    txt = readcell(fullfile(FolderPath, [infoFile, ext]), 'Delimiter', ':', 'NumHeaderLines',1);
end
    
% Rebuild strings that were split by the delimiter:
b_hasMissingVals = cellfun(@(x) isa(x,'missing'), txt);
if size(txt,2)>2
    for i = 1:size(txt,1)
        if all(~b_hasMissingVals(i,2:end))            
            txt{i,2} = strjoin(cellfun(@num2str,txt(i,2:end), 'UniformOutput',false),': ');
        end
    end
end
% Check if an "Events" table exists in "info.txt" file:
hasEvntTable = find(strcmpi(txt(:,1), 'events description'));
if any(hasEvntTable)
    % Read table located beneath the string "Events Description:":
    evntTable = readcell(fullfile(FolderPath, 'info.txt'), 'Delimiter', '\t', 'NumHeaderLines',hasEvntTable+2);
end
% Remove Parameters with missing values:
txt(b_hasMissingVals(:,2),:) = [];
% Replace white spaces in Parameters column by underscores:
txt(:,1) = cellfun(@(x) strrep(x, ' ', '_'),txt(:,1), 'UniformOutput', false);
% Create structure from text:
out = struct();
bCamIndex = 0;
for i = 1:size(txt,1)
    Param = txt{i,1};
    Value = txt{i,2};
    % Parse "Illumination" string:
    if( startsWith(Param, 'Illumination') )
        if( endsWith(Param, 'CameraIdx') )
            bCamIndex = 1;
            out.(Param(1:13)).CamIdx = Value;
        elseif( endsWith(Param, 'FrameIdx') )
            out.(Param(1:13)).FrameIdx = Value;
            bCamIndex = 1;
        elseif( regexp(Param, '[0-9]') )
            out.(Param(1:13)).ID = str2double(Param(13:end));
            out.(Param(1:13)).Color = Value;
        end
        continue
    end
    % Parse "Stimulation[1-2]" string:
%     if regexp(Param, 'Stimulation[1-2]')
    % Transform string with array of digits to num. array:
    if ischar(Value)
        if ( all(isstrprop(erase(Value, ' '), 'digit')) )
            Value = str2num(Value);%#ok. "str2double" would return NaN in this case.
        end
    end
    % Save parameter value to structure:
    out.(Param) = Value;
end
% Save MultiCam index:
if( bCamIndex )
    out.MultiCam = 1;
else
    out.MultiCam = 0;
end
% Save Events info:
if any(hasEvntTable)
    header = matlab.lang.makeValidName(evntTable(1,:)); % Create valid fieldnames for structure.
    for i = 2:size(evntTable,1) % Skip header
        for j = 1:size(evntTable,2)
            out.(['Stim' num2str(i-1)]).(header{j}) = evntTable{i,j};
        end
    end
end