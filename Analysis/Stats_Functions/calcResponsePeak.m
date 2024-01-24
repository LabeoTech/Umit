function outData = calcResponsePeak(data, varargin)
% CALCRESPONSEPEAK calculates the amplitude, peak and onset latencies of a
% response. Here, a "response" means a signal that crosses the "baseline"
% level by a threshold (1 std by default). 
% The input data for this function should be the output of the function
% "getDataFromROI". Thus, the input data is a time vector split by events
% representing the temporal response profile of each ROI to a given stimulus.

% Inputs:
%   data: structure containing the data and meta data created by the
%       "getDataFromROI" function.
%   Parameters:
%       - !!TDB!!
% Output:
%   outData: structure containing stats-ready data extracted from ROIs.

% Defaults:
default_Output = 'PeakStats.mat'; %#ok This line is here just for Pipeline management.
default_opts = struct('STD_threshold', 1, 'ResponsePolarity','positive');
opts_values = struct('STD_threshold',[eps Inf], 'ResponsePolarity', {{'positive','negative'}});%#ok This is here only as a reference for PIPELINEMANAGER.m.
default_object = ''; % This line is here just for Pipeline management to be able to detect this input.
%%% Arguments parsing and validation %%%
p = inputParser;
addRequired(p,'data',@(x) isstruct(x)); % Validate if the input is a structure:
% Optional Parameters:
addOptional(p, 'opts', default_opts,@(x) isstruct(x));
addOptional(p, 'object', default_object, @(x) isempty(x) || isa(x,'Acquisition') || isa(x,'Modality')); % Used by the umIToobox app ONLY.

% Parse inputs:
parse(p,data, varargin{:});
% Initialize Variables:
outData = p.Results.data;
opts = p.Results.opts;
object = p.Results.object;
clear p
%%%%%%%%%%%%%%%%
% Check if the input data is "time vector split by events":
assert(all(ismember(upper(outData.dim_names), {'O','E','T'})),...
    'Wrong input data type. Data must be a time vector split by event(s).');
% For each ROI, calculate the average response and use the average time
% vector to calculate the peak stats:
ROIdata = cell(1,length(outData.data));
evntList = unique(outData.eventID);
for ii = 1:length(outData.data)
        
    % Instantiate output arrays:
    PeakAmp_arr = nan(size(evntList),'single');
    PeakLat_arr = PeakAmp_arr; onsetAmp_arr = PeakAmp_arr; onsetLat_arr = PeakAmp_arr;
    
    for jj = 1:length(evntList)
        idx = ( outData.eventID == evntList(jj) );
        % Calculate average response to the current event:
        avgResp = squeeze(mean(outData.data{ii}(:,idx,:),2,'omitnan'));
        if strcmpi(opts.ResponsePolarity, 'negative')
            % Flip the data around its mean to make the response peak
            % positive:
            avgResp = -1.*avgResp + 2*mean(avgResp,1);
        end
        % Calculate threshold as a multiple of the standard deviation of the 
        %   baseline period from the average response:
        evntFr = round(outData.preEventTime_sec*outData.Freq);
        avgBsln = mean(avgResp(1:evntFr),'omitnan'); % Average baseline amplitude        
        thr = avgBsln + opts.STD_threshold*std(avgResp(1:evntFr),0,1, 'omitnan');
        % Find the maximum ('peak') response value that crosses the
        % threshold:
        [PeakValue,indxPeak] = max(avgResp(evntFr+1:end),[],1);
        indxPeak = evntFr + indxPeak;
        if PeakValue < thr
            % Skip step when no peak is found
            continue
        end
        % Calculate Peak amplitude:
        PeakAmp_arr(jj) = PeakValue - avgBsln;
        % Calculate the onset(threshold crossing point) amplitude and latencies in seconds:
        onsetIndx = evntFr + find(avgResp(evntFr+1:end) > thr,1,'first');
        % Calculate onset amplitude:
        onsetAmp_arr(jj) = avgResp(onsetIndx) - avgBsln;
        % Calculate the Peak latency:
        PeakLat_arr(jj) = (indxPeak - evntFr)/outData.Freq;
        % Calculate onset latency:        
        onsetLat_arr(jj) = (onsetIndx - evntFr)/outData.Freq;
    end
    % Put data inside "ROIdata" structure:
    ROIdata{ii} = PeakAmp_arr;    
end
% Update meta data:
outData.eventID = evntList;
outData.label = outData.eventNameList;
outData.data = ROIdata;

end