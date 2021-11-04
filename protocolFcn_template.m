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

% Change current directory to the "main directory":
cd(obj.MainDir);

% Generate list of all Recording folders:
FolderNames = dir('**\21ST*\M*'); % *All "Day folders" starts with "21ST" and all "recording folders" starts with the letter "M". 
FolderNames = FolderNames([FolderNames.isdir] == 1); % Get only folders. 
FolderNames = cellfun(@(x,y) fullfile(x, y), {FolderNames.folder}, {FolderNames.name}, 'UniformOutput', false); % Extract folder names:

% Parse recording folder names:
idx_underscore = cellfun(@(x) strfind(x,'_'), FolderNames, 'UniformOutput', false); % Find underscores. This will be used later.
% Get Sunject IDs:
subjS = regexp(FolderNames, 'M\w+?(?=\_)', 'match', 'once'); % Get characters before the first underscore.
uniqS= unique(subjS); % Get the list of all subjects.
% Create emtpy array of Subject objects:
out = repmat(Subject(),1, numel(uniqS)); 
% Loop across all Subjects and populate each one with Acquisition and
% modality objects:
for i = 1:length(uniqS)
    out(i).ID = uniqS{i}; % Set Subjet "ID";
    out(i).MyParent = obj; % Set protocol object (obj) as the Subject Parent. !! Important: This step is essential for the protocol object to work. Do not change or erase this line.
    % Find all aquisitions from a Subject:
    idxAcq = contains(FolderNames, [uniqS{i} '_']);
    AcqList = cellfun(@(x,y) x(y(1)+1:end), FolderNames(idxAcq), idx_underscore(idxAcq), 'UniformOutput', false); % Get the Acquisition names
    subFolders = FolderNames(idxAcq); % Get list of folders containing the recordings for a Subject.
    % Loop across the list of acquisitions:
    for j = 1:length(AcqList)
        % Instantiate an Acquisition object:
        tmpA = Acquisition();
        tmpA.MyParent = out(i); % Set the current Subject as the Acquisition Parent.
        tmpA.ID = AcqList{j}; % Set the Aquisition ID.
        % Add Cortical Imaging Object to the Acquisition object:
        tmpA.Array.addObj(FluorescenceImaging()); % Adds a FluorescenceImaging object to the current Acquisition.
        tmpA.Array.ObjList(end).MyParent = tmpA; % Set the current Acquisition as the FluorescenceImaging object Parent.
        tmpA.Array.ObjList(end).ID = 'CortexImagingData'; % Set the FluorescenceImaging object ID.
        tmpA.Array.ObjList(end).RawFolder = subFolders{j}; % Set FluorescenceImaging RawFolder property.
        % Get a list of raw data files from the raw folder using regular expression: 
        files = getNamesFromDir(subFolders{j}, expLabeo, 0, 'match'); 
        tmpA.Array.ObjList(end).RawFiles = files; % Set RawFiles list property.
        % Find text file to get Aquisition Start_datetime and other info:
        infoFile = fullfile(subFolders{j}, 'info.txt'); % File full path
        txt = fileread(infoFile); % Read text file.
        DateTime = regexp(txt, '(?<=DateTime:\s*)\S+','match', 'once'); % Get date and time as string.
        DateTime = datetime(DateTime, 'InputFormat', 'yyyyMMdd_HHmmss'); % Transform string to datetime object.
        tmpA.Start_datetime = DateTime; % Set Start_datetime property.
        
        % Get number of recording channels (optional):
        numChan = numel(regexp(txt, 'Illumination', 'match'));
        
        % Get sample rate
        sr = str2double(regexp(txt, '(?<=FrameRateHz:\s+)\d+\.*\d*', 'match'));
        
        % Set optional parameters:
        tmpA.Array.ObjList(end).SampleRateHz = sr;
        tmpA.Array.ObjList(end).NumberOfChannels = numChan;       
        tmpA.Array.ObjList(end).RecordingSystem = 'LabeoTechSystem';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % If the recordings include event-triggered datasets, the toolbox
        % needs the full path to a "meta data file" containing the stimulus
        % information:
        matFileList = dir([subFolders{j} '\*.mat']); matFileList = {matFileList.name}; % Find a .MAT file containing the stimulus data.
        idx = contains(matFileList, uniqS{i});
        if sum(idx) > 0
            indx = find(idx);
            matFile = matFileList{indx(1)}; % Picks the first one, but I assume that there is only one file per folder.
            tmpA.Array.ObjList(end).MetaDataFileName = matFile; % Set the MetaDataFileName property of FluorescenceImaging object.
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Finally, add the Acquisition Object to the current Subject object:
        out(i).Array.addObj(tmpA);
    end
end









