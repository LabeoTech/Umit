function varargout = HemoCorrection(Folder, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Infos:
%
% This function is used to remove the hemodynamic fluctuations from any
% fluorescence signal. 
%
% Inputs:
% 1. Folder: Folder contaning the dataset to work with.
% 2. Varargin -> if empty: a dialog box will be prompt to ask user which
%                   channels to use to do the correction.
%             -> cell array of string: to specify which channels to use.
%                   Ex: HemoCorrection(pwd, {'Red', 'Green'});
% 3. Optional: lowpass filter 
%           -> Value of the cutoff frequency to use on a lowpass filter
%           applied to intrinsic signals. THIS PARAMETER IS OPTIONAL
%           To keep the data as is, just use the function with 2 parameters
% Ouput:
% - If an output is set, the result of the correction will be given back
% through this output. All the data in the folder will remain unchanged.
% - If no output is specified, the fChan.dat files in Folder will be overwritten
% with the corrected data.
% 
% Exemples:
%
% 1- HemoCorrection(pwd); 
% The fluorescence dat files in the folder will be overwriten and a dialog
% box will be used in order to select which channels must be used to
% compute the correction.
% 2- NewFluo = HemoCorrection(pwd, {'Green'});
% The dat files in the folder won't be overwriten. The new corrected
% fluorescence data will be in NewFluo. Only the green channel will be used
% to compute the correction.
% 3- HemoCorrection(pwd, {'Red, 'Green', 'Amber'});
% fChan_475.dat will be overwriten with the corrected fluorescence data.
% All three hemodynamic channels will be used to compute the correction.

if( ~strcmp(Folder(end),filesep) )
    Folder = strcat(Folder, filesep);
end

if( nargin <= 1 )
    cList = dir([Folder '*.dat']);
    fn = {};
    for ind = 1:size(cList,1)
        if( ~strcmp(cList(ind).name(1),'f') )
            fn{end+1} = cList(ind).name;
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
    tmp = varargin{1};
    for ind = 1:size(tmp,2)
        tag = lower(tmp{ind});
        switch tag
            case 'red'
                if( exist([Folder 'rChan.dat'], 'file') )
                    fn{end+1} = 'rChan.dat';
                else
                    fn{end+1} = 'red.dat';
                end
            case {'amber', 'yellow'}
                if( exist([Folder 'yChan.dat'], 'file') )
                    fn{end+1} = 'yChan.dat';
                else
                    fn{end+1} = 'yellow.dat';
                end
            case 'green'
                if( exist([Folder 'gChan.dat'], 'file') )
                    fn{end+1} = 'gChan.dat';
                else
                    fn{end+1} = 'green.dat';
                end
        end
    end
end
fList = dir([Folder 'fluo*.mat']);
if( isempty( fList ) )
    fList = dir([Folder 'Data_Fluo*.mat']);
end
idx = arrayfun(@(x) contains(fList(x).name,'splitAndNorm'), 1:size(fList,1));
fList(idx) = [];

Infos = matfile([Folder fList(1).name]);
NbFrames = Infos.datLength; 
HemoData = zeros(size(fn,2), prod(Infos.datSize), NbFrames, 'single');

if( nargin <= 2 )
    bFilt = false;
    sFreq = Infos.Freq/2;
else
    sFreq = varargin{2};
    bFilt = true;
    if( sFreq >= Infos.Freq/2 )
        bFilt = false;
        sFreq = Infos.Freq/2;
    end
end
%Reading Hemo file:
for ind = 1:size(fn,2)
    fprintf('Opening: %s \n', fn{ind});
    eval(['fid = fopen(''' Folder fn{ind} ''');']);
    tmp = fread(fid, inf, '*single');
    tmp = reshape(tmp, Infos.datSize(1,1)*Infos.datSize(1,2), []);
    if( bFilt )
        fprintf('Time Filtering...\n');
        f = fdesign.lowpass('N,F3dB', 4, sFreq, Infos.Freq);
        lpass = design(f,'butter');
        tmp = single(filtfilt(lpass.sosMatrix, lpass.ScaleValues, double(tmp')))';
    end
    tmp = reshape(tmp, Infos.datSize(1,1), Infos.datSize(1,2), []);
    fprintf('Spatial Filtering...\n');
    tmp = imgaussfilt(tmp,1, 'Padding', 'symmetric');
    tmp = reshape(tmp, [], size(tmp,3));
    tmp = (tmp - mean(tmp,2))./mean(tmp,2);
    HemoData(ind, :, :) = tmp;
    fprintf('Done.\n');
    fclose(fid);
end
if( size(HemoData, 2) ~= Infos.datLength )
    sz = size(HemoData);
    dimCanaux = find(sz == size(fn,2));
    dimPix = find(sz == prod(Infos.datSize));
    dimTime = find(sz == prod(Infos.datLength));
    
    HemoData = permute(HemoData,[dimCanaux, dimPix, dimTime]);
end
clear tmp fn fid ind NbPts

%Correction:
fList = dir([Folder 'fluo*.dat']);
if( isempty( fList ) )
    fList = dir([Folder 'fChan*.dat']);
end
idx = arrayfun(@(x) contains(fList(x).name,'splitAndNorm'), 1:size(fList,1));
fList(idx) = [];
for ind = 1:size(fList,1)
    fprintf('Opening fluo data: %s \n', fList(ind).name);
    eval(['fid = fopen(''' Folder fList(ind).name ''');']);
    fData = fread(fid, inf, '*single');
    fclose(fid);
    fData = reshape(fData, prod(Infos.datSize), []);
    m_fData = mean(fData,2);
    fData = (fData - m_fData)./m_fData;
    
    fprintf('Hemodynamic Correction: ');
    warning('off', 'MATLAB:rankDeficientMatrix');
    h = waitbar(0, 'Fitting Hemodyn on Fluorescence');

    for indF = 1:size(fData,1)       
        X = [ones(1, size(fData,2)); linspace(0,1,size(fData,2)); reshape(HemoData(:,indF,:),sz(dimCanaux), sz(dimTime))];      
        B = X'\fData(indF,:)';
        A(indF,:) = B;
        fData(indF,:) = fData(indF,:) - (X'*B)';
        waitbar(indF/size(fData,1), h);
    end
    close(h);
    warning('on', 'MATLAB:rankDeficientMatrix');
    clear B X locHemoDyn;
    
    fData = bsxfun(@times, fData, m_fData) + m_fData;
    
    fData = reshape(fData, Infos.datSize(1,1), Infos.datSize(1,2), []);
    if( nargout == 0 )
        eval(['fid = fopen(''' Folder fList(ind).name ''', ''w'');']);
        fwrite(fid, fData, 'single');
        fclose(fid);
    else
        varargout{ind} = fData;
    end 
end

end