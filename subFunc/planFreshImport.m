function plan = planFreshImport(SaveFolder)
%PLANFRESHIMPORT List what a fresh import of SaveFolder would delete and keep.
%
%   plan = planFreshImport(SaveFolder)
%
%   A re-import into an existing SaveFolder is a FRESH IMPORT: every
%   artefact derived from the previous data is deleted, so the folder is
%   equivalent to a first import (brief section 1.6). Read-only: nothing is
%   written or deleted here (see applyFreshImport).
%
%   Deletion scope (brief section 4 Q2, Bruno 2026-10-15): only files UMIT
%   KNOWS it created, top level only:
%       - import output: *.dat, AcqInfos.mat, ROI.mat, tformDualCam.mat and
%         legacy sidecars <channel>.mat / <channel>_info.mat;
%       - folder state and derived files: *.umt, *.roi, DataParams.mat,
%         events.mat, dataHistory.mat, pipeLog.mat,
%         manualAlignmentReport_*.mat, *MotionCorrection.mat,
%         registrationQC_*.fig/.png, *MotionCorrectionQC.png,
%         ImageAlignmentFailedReport_*.txt, pipeErrorLog_*.csv, DataViewer
%         ROI exports *_ROI_*.csv, std_speckle.tiff, ConvertToTiff exports
%         of the folder's channels (<channel>.tif, <channel>_C*_R*.tif,
%         <channel>_info.txt), legacy Text_events.mat, LogBook.mat,
%         ImagingReferenceFrame.mat;
%       - every file recorded as an output in dataHistory.mat;
%       - temporary leftovers: PMTMP_*, *.dat.tmp, *.aligning.tmp,
%         *.aligning.bak, saveMatAtomic/ledger temps tp<uuid>[.bak].mat|json.
%   Always kept:
%       - raw acquisition files: *.bin, info.txt / *info.txt (except the
%         ConvertToTiff exports above), Comments.txt, Snapshot*.png, input
%         *.tif / *.tiff (except the exports above), info.json;
%       - ledgers: UMITRigLedger.json, UMITProjectLedger.json, their lock
%         folders, a legacy UMITProjectBinding.umitlink, and quarantined
%         invalid rig ledgers (UMITRigLedger.invalid_<timestamp>.json);
%       - backup zips (*.zip) of earlier imports;
%       - anything else (unknown files and every subfolder), reported in
%         plan.unknown so a confirmation can list it.
%
%   Output PLAN (scalar struct):
%       saveFolder       - absolute folder
%       hasPreviousImport - true when the folder holds an import or any
%                          deletable artefact or a rig/project ledger: the
%                          next import is then a fresh import
%       toDelete         - cellstr of file names (relative) to delete
%       keptRaw, keptLedgers, keptBackups, unknown - cellstr of names kept
%       rigLedgerExists, projectLedgerBound - ledgers that get a
%                          freshImport event
%
%   See also applyFreshImport, ImagesClassification, importFromTif.

SaveFolder = char(string(SaveFolder));
plan = struct('saveFolder', SaveFolder, 'hasPreviousImport', false, ...
    'toDelete', {{}}, 'keptRaw', {{}}, 'keptLedgers', {{}}, 'keptBackups', {{}}, ...
    'unknown', {{}}, 'rigLedgerExists', false, 'projectLedgerBound', false);
if ~isfolder(SaveFolder)
    return
end

listing = dir(SaveFolder);
listing = listing(~ismember({listing.name}, {'.', '..'}));
names = {listing.name};
isDir = [listing.isdir];

datNames = names(~isDir & endsWith(names, '.dat', 'IgnoreCase', true));
channelBases = cellfun(@(n) n(1:end-4), datNames, 'UniformOutput', false);
historyOutputs = iHistoryOutputs(SaveFolder);

ledgerNames = {'UMITRigLedger.json', 'UMITProjectLedger.json', ...
    'UMITProjectBinding.umitlink', 'UMITRigLedger.lock', 'UMITProjectLedger.lock'};
for iEntry = 1:numel(names)
    name = names{iEntry};
    if ismember(name, ledgerNames) || ...
            ~isempty(regexp(name, '^UMITRigLedger\.invalid_\d{8}T\d{6}(_\d+)?\.json$', 'once'))
        plan.keptLedgers{end+1} = name;
    elseif isDir(iEntry)
        plan.unknown{end+1} = [name filesep];
    elseif iIsDeletable(name, channelBases, historyOutputs)
        plan.toDelete{end+1} = name;
    elseif iIsRaw(name)
        plan.keptRaw{end+1} = name;
    elseif endsWith(name, '.zip', 'IgnoreCase', true)
        plan.keptBackups{end+1} = name;
    else
        plan.unknown{end+1} = name;
    end
end

plan.rigLedgerExists = isfile(fullfile(SaveFolder, 'UMITRigLedger.json'));
plan.projectLedgerBound = UMITProjectLedger.isBound(SaveFolder);
plan.hasPreviousImport = ~isempty(plan.toDelete) || plan.rigLedgerExists || ...
    UMITProjectLedger.exists(SaveFolder);
end

function tf = iIsDeletable(name, channelBases, historyOutputs)
lowerName = lower(name);
[~, base, ext] = fileparts(name);
ext = lower(ext);

fixedNames = lower({'AcqInfos.mat', 'ROI.mat', 'tformDualCam.mat', 'DataParams.mat', ...
    'events.mat', 'dataHistory.mat', 'pipeLog.mat', 'Text_events.mat', 'LogBook.mat', ...
    'ImagingReferenceFrame.mat', 'std_speckle.tiff'});
if ismember(lowerName, fixedNames) || ismember(ext, {'.dat', '.umt', '.roi'}) || ...
        any(strcmpi(name, historyOutputs))
    tf = true;
    return
end

patterns = { ...
    '^manualalignmentreport_.*\.mat$', '^.*motioncorrection\.mat$', ...
    '^registrationqc_.*\.(fig|png)$', '^.*motioncorrectionqc\.png$', ...
    '^imagealignmentfailedreport_.*\.txt$', '^pipeerrorlog_.*\.csv$', ...
    '^.*_roi_.*\.csv$', '^pmtmp_.*', '^.*\.dat\.tmp$', '^.*\.aligning\.(tmp|bak)$', ...
    '^tp[0-9a-f_]{36}(\.bak)?\.(mat|json)$'};
if any(~cellfun(@isempty, regexp(lowerName, patterns, 'once')))
    tf = true;
    return
end

% Legacy sidecars and ConvertToTiff exports of the folder's own channels.
tf = false;
for iBase = 1:numel(channelBases)
    channel = lower(channelBases{iBase});
    if ismember(lowerName, {[channel '.mat'], [channel '_info.mat'], [channel '.tif'], ...
            [channel '_info.txt']}) || ...
            ~isempty(regexp(lowerName, ['^' regexptranslate('escape', channel) ...
            '_c.*_r.*\.tif$'], 'once'))
        tf = true;
        return
    end
end
if strcmp(ext, '.mat') && any(strcmpi([base '.dat'], channelBases))
    tf = true;
end
end

function tf = iIsRaw(name)
lowerName = lower(name);
tf = endsWith(lowerName, '.bin') || endsWith(lowerName, 'info.txt') || ...
    strcmp(lowerName, 'comments.txt') || strcmp(lowerName, 'info.json') || ...
    (startsWith(lowerName, 'snapshot') && endsWith(lowerName, '.png')) || ...
    endsWith(lowerName, '.tif') || endsWith(lowerName, '.tiff');
end

function outputs = iHistoryOutputs(SaveFolder)
% Output file names recorded in dataHistory.mat (best effort, read-only).
outputs = {};
historyPath = fullfile(SaveFolder, 'dataHistory.mat');
if ~isfile(historyPath)
    return
end
try
    loaded = load(historyPath);
catch
    return
end
values = struct2cell(loaded);
for iValue = 1:numel(values)
    outputs = [outputs, iCollectFileNames(values{iValue})]; %#ok<AGROW>
end
outputs = unique(outputs);
end

function names = iCollectFileNames(value)
names = {};
if isstruct(value)
    fields = fieldnames(value);
    for iItem = 1:numel(value)
        for iField = 1:numel(fields)
            % PipelineManager history entries store the output in 'file'.
            if any(strcmpi(fields{iField}, {'file', 'outFile', 'outputFile', 'FileName', ...
                    'outputFileName', 'savedFile'}))
                names = [names, iTextNames(value(iItem).(fields{iField}))]; %#ok<AGROW>
            end
        end
    end
elseif istable(value)
    for name = intersect(value.Properties.VariableNames, {'FileName', 'fileName', 'outFile'})
        names = [names, iTextNames(value.(name{1}))]; %#ok<AGROW>
    end
end
end

function names = iTextNames(value)
names = {};
if ischar(value) || isstring(value) || iscellstr(value)
    text = cellstr(string(value));
    for iText = 1:numel(text)
        [~, base, ext] = fileparts(text{iText});
        if ~isempty(base)
            names{end+1} = [base ext]; %#ok<AGROW>
        end
    end
end
end
