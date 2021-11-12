function outData = run_HemoCorrection(SaveFolder, varargin)
% RUN_HEMOCORRECTION calls the function
% HEMOCORRECTION from the IOI library (LabeoTech).
% In brief, this function removes hemodynamic fluctuations from any
% fluorescence signal using the information of two or more channels.

% The algorithm used here is described in:
% Valley, Matthew & Moore, Michael & Zhuang, Jun & Mesa, Natalia & Castelli, Dan & Sullivan,
% David & Reimers, Mark & Waters, Jack. (2019). Separation of hemodynamic signals from
% GCaMP fluorescence measured with widefield imaging. Journal of Neurophysiology. 123.
% 10.1152/jn.00304.2019.


% Defaults:
default_Output = 'hemoCorr_fluo.dat'; %#ok. This line is here just for Pipeline management.

%%% Arguments parsing and validation %%%
p = inputParser;
% Save folder:
addRequired(p, 'SaveFolder', @isfolder);
% Optional Parameters:
% opts structure:
default_opts = struct('Red', true, 'Green', true, 'Amber', true);
addOptional(p, 'opts', default_opts,@(x) isstruct(x) && ~isempty(x));
% Parse inputs:
parse(p,SaveFolder, varargin{:});

% Translate opts to char cell array:
fields = fieldnames(p.Results.opts);
idx = cellfun(@(x) p.Results.opts.(x), fields);
list = fields(idx)';

% Run HemoCorrection function from IOI library:
disp('Performing hemodynamic correction in fluo channel...')
outData = HemoCorrection(p.Results.SaveFolder, list);
disp('Finished hemodynamic correction.')
end