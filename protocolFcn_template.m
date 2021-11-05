function out = protocolFcn_template(obj)
% User-defined function. This Function is used inside "Protocol" object to generate the
% list of Subjects, Acquisitions and Modalities.
% This function will be used to update the content of the obj.MainDir (Folder where the Data is).
% Input:
% obj (Protocol class) = self referring handle of Protocol class.
% Output:
% out (Subject class) = array of objects from Subject class containing Acquisition and
% Modalities.

% This is a template for creating a protocol function.
% This template reads mesoscale imaging datasets created with the imaging
% systems from Labeo Technologies Inc. The imaging raw data are stored in a
% series of binary files (.bin). This template function scans a directory
% (obj.MainDir) containing folders with the imaging data and a text file
% with the recording information.
% The folder tree is organized as :
%      MainDir --
%                |
%                |--RECORDINGFOLDER-- (ex. "M2D_RS_21ST0401")
%                                                   |
%                                                   |-xxxx.bin
%                                                   |-xxxx.bin
%                                                   |-xxxx.bin
%                                                      ...
%                                                   |-xxxx.bin
%                                                   |-info.txt

% The names of subjects (Subject ID) and acquisitions (Acquisition ID) are encoded in the
% folder names separated by an underscore as SUBJID "_" ACQID. For example:
% Lets take the example folder above (M2D_RS_21ST0401), here we have:
% Subject ID: "M2D"
% Acquisition ID: "RS_21ST0401

% Minimal information required:
%   In order to make this function work, we need to set the following properties
%   for each object created in this function:
%       Subject:
%               -ID (str)
%       Acquisition:
%               -ID (str)
%               -Start_datetime (datetime): Timestamp of the beginning of
%               the recording.
%       Modality:
%               -ID (str)
%               -RawFolder (str): Full path for an existing folder
%                containing raw data from a given modality.
%               -RawFiles (cell): list of raw file names with the file extensions.
% * Other properties can be set as well. Please, check the documentation
% for the classes "Subject", "Acquisition", "Modality" and those derived
% from it.

% Regular expression to find the raw data files:
expLabeo = '(ai|img)_\d*.bin'; % For Labeo Data
% Generate list of all Recording folders:
folder_list = dir(obj.MainDir); folder_list = folder_list(~ismember({folder_list.name}, {'.','..'}) & [folder_list.isdir] == 1);
% Keep folders that start with upper letters and/or numbers followed by an underscore:
idx = cellfun(@(x) ~isempty(x), regexp({folder_list.name}, '[A-Z0-9]*_.*')); folder_list = folder_list(idx);
% Parse recording folder names to get subject and acquisition IDs:
IDs = regexp({folder_list.name}', '_', 'split', 'once'); IDs = vertcat(IDs{:});
% Here, we create an array of subjects based on the IDs extracted from the
% folder names in which all the acquisition and modality objects will be
% stored in:
uniqS = unique(IDs(:,1)); % Get the list of all subjects.
% Create an array of emtpy Subject objects:
for i = 1:numel(uniqS)
    out(i) = Subject();%#ok
end
% Loop across all folders and populate Subject objects with with Acquisition and modality objects:
cnt = 1;
for i = 1:length(folder_list)
    % Try to find existing Subject in the array. If not, add its ID to it.
    idx_S = strcmp(IDs(i,1),{out.ID});
    if ~any(idx_S)
        % Add Subject ID to object (other info can be also added here):
        idx_S = cnt;
        out(idx_S).ID = IDs{i,1};
        cnt = cnt+1;
    end
    % Create an Acquisition object:
    tmpA = Acquisition();
    tmpA.ID = IDs{i,2}; % Set the Aquisition ID.
    % Add Cortical Imaging Object to the Acquisition object:
    tmpA.Array.addObj(FluorescenceImaging()); % Adds a FluorescenceImaging object to the current Acquisition.
    tmpA.Array.ObjList(end).ID = 'CortexImaging'; % Set the FluorescenceImaging object ID.
    tmpA.Array.ObjList(end).RawFolder = fullfile(folder_list(i).folder, folder_list(i).name); % Set FluorescenceImaging RawFolder property.
    % Get a list of raw data files from the raw folder using regular expression:
    files = getNamesFromDir(fullfile(folder_list(i).folder, folder_list(i).name), expLabeo, 0, 'match');
    tmpA.Array.ObjList(end).RawFiles = files; % Set RawFiles list property.
    % Find text file to get Aquisition Start_datetime and other info:
    infoFile = fullfile(folder_list(i).folder, folder_list(i).name, 'info.txt'); % File full path
    txt = fileread(infoFile); % Read text file.
    DateTime = regexp(txt, '(?<=DateTime:\s*)\S+','match', 'once'); % Get date and time as string.
    DateTime = datetime(DateTime, 'InputFormat', 'yyyyMMdd_HHmmss'); % Transform string to datetime object.
    tmpA.Start_datetime = DateTime; % Set Start_datetime property.
    % Get number of recording channels (optional):
    numChan = numel(regexp(txt, 'Illumination', 'match'));
    % Get sample rate (optional)
    sr = str2double(regexp(txt, '(?<=FrameRateHz:\s+)\d+\.*\d*', 'match'));
    % Add optional parameters to Modality object:
    tmpA.Array.ObjList(end).SampleRateHz = sr;
    tmpA.Array.ObjList(end).NumberOfChannels = numChan;
    tmpA.Array.ObjList(end).RecordingSystem = 'LabeoTechSystem';
    % Finally, add the Acquisition Object to the current Subject object:
    out(idx_S).Array.addObj(tmpA);
end

end