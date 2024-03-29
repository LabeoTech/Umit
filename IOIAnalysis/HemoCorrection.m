function [fData, fMetaData] = HemoCorrection(Folder, fData, bSave, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Infos:
%
% This function is used to remove the hemodynamic fluctuations from any
% fluorescence signal.
%
% Inputs:
% 1. Folder: Folder contaning the dataset to work with.
% 2. fData (3D num. array | char): Image time series (with dimensions "Y","X","T")
%       with the fluorescence channel to be corrected OR name of the .DAT file
%       containing the fluorescence signal.
% 3. bSave (logical): If TRUE, this function saves the corrected
%       fluorescence signal to a .DAT file named "fluo_corrected.dat".
% Optional inputs:
%
% 4. ChannelList(cell): list of reflectance channel names to be used by the
%       Hemodynamic correction algorithm. It must be two or more from the
%       list {"Red","Amber","Green"}. If not provided (default), a dialog box
%       will appear to select the channels.
% 5. fMetaData (struct): Name-value pair argument. Meta data structure
%       associated with the fluorescence signal. This parameter must be provided
%       if the input "fData" is a 3D array.
% 6. LowPassFreq(num scalar): Name-value pair argument. Cut-off frequency
%       in Hertz for a low pass filter applied to the reflectance signals prior
%       to the hemodynamic correction.
% Ouput:
%   - fData(3D num. array):  the result of the correction will be given back
%       through this output.
%   - fMetaData (struct): meta data associated with the signal stored in
%       fData.

% Examples:
%
%   1. Using the function with minimal parameters, prompting the user to
%      select channels and saving the corrected data:
%
%      Folder = 'DataFolder';
%      fData = 'fluodata.dat';
%      bSave = true;
%      HemoCorrection(Folder, fData, bSave);
%
%   2. Using the function with a specific list of channels:
%
%      Folder = 'DataFolder';
%      fData = 'fluodata.dat';
%      bSave = true;
%      ChannelList = {'Red', 'Green'};
%      HemoCorrection(Folder, fData, bSave, ChannelList);
%
%   3. Using the function with a low-pass filter and saving the corrected
%      data:
%
%      Folder = 'DataFolder';
%      fData = 'fluodata.dat';
%      bSave = true;
%      LowPassFreq = 1.0; % 1 Hz cutoff frequency for the low-pass filter
%      HemoCorrection(Folder, fData, bSave, 'LowPassFreq', LowPassFreq);
%
%   4. Using the function to perform the correction and receive the
%      corrected data as output:
%
%      Folder = 'DataFolder';
%      fData = 'fluodata.dat';
%      bSave = false; % Do not save the corrected data to a file
%      [fData, fMetaData] = HemoCorrection(Folder, fData, bSave);
%      % The corrected data is now available in the "fData" variable
%
%   5. Using the function to perform the correction using the fluorescence
%       data as input and output the corrected data.
%      [fData, fMetaData] = loadDatFile('fluodata.dat');
%        bSave = false;
%      [fData,fMetaData] = HemoCorrection(Folder,fData,bSave,'fMetaData',fMetaData);


%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'Folder',@isfolder)
addRequired(p,'fData',@(x) ischar(x) || isnumeric(x))
addRequired(p,'bSave',@isscalar)
addOptional(p,'ChannelList',{''},@(x) iscell(x) && ischar([x{:}]));
addParameter(p,'fMetaData',[],@(x) isa(x,'matlab.io.MatFile') | isstruct(x));
addParameter(p,'LowPassFreq',[],@isscalar)
% Parse inputs:
parse(p,Folder,fData,bSave,varargin{:});
% Instantiate input parameters:
Folder = p.Results.Folder;
fData = p.Results.fData;
bSave = p.Results.bSave;
cList = p.Results.ChannelList;
fMetaData = p.Results.fMetaData;
sFreq = p.Results.LowPassFreq;
clear p
%

if( ~strcmp(Folder(end),filesep) )
    Folder = strcat(Folder, filesep);
end
% Get list of channels to be used:
if( isempty([cList{:}]) )
    cList = dir([Folder '*.dat']);
    fn = {};
    for ind = 1:size(cList,1)
        if( ~strcmp(cList(ind).name(1),'f') )
            fn{end+1} = cList(ind).name;%#ok
        end
    end
    
    [idx, tf] = listdlg('PromptString',{'Select channels to be used to',...
        'compute hemodynamic correction.',''},...
        'ListString',fn);
    
    if( tf == 0 )
        return;
    end
    
    fn = fn(idx);
    clear cList idx ind tf;
else
    fn = {};
    for ind = 1:length(cList)
        tag = lower(cList{ind});
        switch tag
            case 'red'
                if( exist([Folder 'rChan.dat'], 'file') )
                    fn{end+1} = 'rChan.dat'; %#ok
                else
                    fn{end+1} = 'red.dat';%#ok
                end
            case {'amber', 'yellow'}
                if( exist([Folder 'yChan.dat'], 'file') )
                    fn{end+1} = 'yChan.dat';%#ok
                else
                    fn{end+1} = 'yellow.dat';%#ok
                end
            case 'green'
                if( exist([Folder 'gChan.dat'], 'file') )
                    fn{end+1} = 'gChan.dat';%#ok
                else
                    fn{end+1} = 'green.dat';%#ok
                end
        end
    end
end
% Open fluo data file:
if ischar(fData)
    fprintf('Loading non-corrected fluo file "%s"...\n',fData);
    [~,file,~] = fileparts(fData);
    fDataFile= [file,'.dat']; % Enforce file extension to .DAT
    fMetaData = load([Folder, file, '.mat']);
    % Load fluo file:    
    eval(['fid = fopen(''' Folder fDataFile ''');']);
    fData = fread(fid, inf, '*single');%#ok
    fData = reshape(fData, fMetaData.datSize(1,1)*fMetaData.datSize(1,2), []);
    fData = reshape(fData, fMetaData.datSize(1,1), fMetaData.datSize(1,2), []);
    fclose(fid);
    fprintf('Done.\n')
else
    if isempty(fMetaData)
        error('MetaData .MAT file must be provided when the data is provided as input!')
    end
end

NbFrames = fMetaData.datLength;
HemoData = zeros(size(fn,2), prod(fMetaData.datSize), NbFrames, 'single');

if( isempty(sFreq) )
    bFilt = false;
    sFreq = fMetaData.Freq/2;
else
    bFilt = true;
    if( sFreq >= fMetaData.Freq/2 )
        bFilt = false;
        sFreq = fMetaData.Freq/2;
    else
        fprintf('Low pass temporal filter cut-off set at %0.1f Hz.\n',sFreq);
    end
end
% Reading Hemo file:
for ind = 1:size(fn,2)
    fprintf('Opening: %s \n', fn{ind});
    eval(['fid = fopen(''' Folder fn{ind} ''');']);
    tmp = fread(fid, inf, '*single');
    tmp = reshape(tmp, fMetaData.datSize(1,1)*fMetaData.datSize(1,2), []);
    if( bFilt )
        fprintf('Time Filtering...\n');
        f = fdesign.lowpass('N,F3dB', 4, sFreq, fMetaData.Freq);
        lpass = design(f,'butter');
        tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
    end
    tmp = reshape(tmp, fMetaData.datSize(1,1), fMetaData.datSize(1,2), []);
    fprintf('Spatial Filtering...\n');
    tmp = imgaussfilt(tmp,1, 'Padding', 'symmetric');
    tmp = reshape(tmp, [], size(tmp,3));
    tmp = (tmp - mean(tmp,2))./mean(tmp,2);
    HemoData(ind, :, :) = tmp;
    fprintf('Done.\n');
    fclose(fid);
end
if( size(HemoData, 2) ~= fMetaData.datLength )
    sz = size(HemoData);
    dimCanaux = find(sz == size(fn,2));
    dimPix = find(sz == prod(fMetaData.datSize));
    dimTime = find(sz == prod(fMetaData.datLength));
    
    HemoData = permute(HemoData,[dimCanaux, dimPix, dimTime]);
end
clear tmp fn fid ind NbPts

% Correction:
fData = reshape(fData,prod(fMetaData.datSize),[]);
m_fData = mean(fData,2);
fData = (fData - m_fData)./m_fData;
fprintf('Hemodynamic Correction: ');
warning('off', 'MATLAB:rankDeficientMatrix');
h = waitbar(0, 'Fitting Hemodyn on Fluorescence');
for indF = 1:size(fData,1)
    X = [ones(1, size(fData,2)); linspace(0,1,size(fData,2)); reshape(HemoData(:,indF,:),sz(dimCanaux), sz(dimTime))];
    B = X'\fData(indF,:)';
    fData(indF,:) = fData(indF,:) - (X'*B)';
    waitbar(indF/size(fData,1), h);
end
close(h);
warning('on', 'MATLAB:rankDeficientMatrix');
clear B X;
fData = bsxfun(@times, fData, m_fData) + m_fData;
fData = reshape(fData, fMetaData.datSize(1,1), fMetaData.datSize(1,2), []);
fprintf('Finished Hemodyn on Fluorescence.\n')
% Update output meta data file:
[~,filename,~] = fileparts(fMetaData.datFile);
if isempty(filename)
    filename = 'fluo';
end
fMetaData.datFile = [filename '_corrected.dat']; % Update datfile.
% Save corrected data to file:
if ( bSave )
    fprintf('Saving corrected fluo file...\n');      
    % Save .DAT file:
    fid = fopen([Folder fMetaData.datFile],'W');
    fwrite(fid, fData, '*single');
    fclose(fid);
    % Save .MAT file:    
    save([Folder fMetaData.datFile], '-struct','fMetaData')
    fprintf('Corrected fluo data saved in "%s" as "%s"\n',Folder, fMetaData.datFile);
end

end
