function summary = describeFolderRegistration(saveFolder)
%DESCRIBEFOLDERREGISTRATION DataViewer summary of a folder's image alignment.
%
%   summary = describeFolderRegistration(saveFolder)
%
%   Read-only: reads DataParams.registration (the only record of a folder's
%   reference and alignment) and the project ledger, and writes nothing.
%   Never throws for a missing or invalid DataParams: it is reported as a
%   state. SUMMARY fields:
%       folder
%       state            'NoDataParams', 'Invalid', 'NoReference',
%                        'ReferenceSelected', 'Estimated', 'EstimateStale'
%                        (the estimate was made against another reference)
%                        or 'Aligned'
%       issues           load/validation messages (cellstr; 'Invalid' only)
%       isAligned        the folder's data were transformed (committed)
%       referenceName    name of the selected or committed reference
%       sourceKind       'project', 'umitDat', 'file', 'legacy' or ''
%       sourcePath       recorded source path (human aid only)
%       shortHash        committed hash when aligned, else the selected one
%       alignedSourceWarning  aligned-source warning recorded at selection
%       appliedOn, appliedBy, isReviewed
%       isProjectBound   the folder's project ledger is in state Bound
%       displayText      short text for the folder context label
%       tooltipLines     detailed lines (cellstr); in project-less mode they
%                        include a note that projects let a reference be
%                        reused across recordings
%
%   See also describeRegistrationReference, describeFolderRigLedger,
%   describeFolderProjectLedger.

saveFolder = char(string(saveFolder));
summary = struct('folder', saveFolder, 'state', 'NoDataParams', 'issues', {{}}, ...
    'isAligned', false, 'referenceName', '', 'sourceKind', '', 'sourcePath', '', ...
    'shortHash', '', 'alignedSourceWarning', '', 'appliedOn', '', 'appliedBy', '', ...
    'isReviewed', false, 'isProjectBound', false, 'displayText', '', ...
    'tooltipLines', {{}});

if isempty(saveFolder) || ~isfolder(saveFolder)
    summary.displayText = '<none>';
    summary.tooltipLines = {'Alignment: no accessible save folder.'};
    return
end
summary.isProjectBound = UMITProjectLedger.isBound(saveFolder);

if ~isfile(fullfile(saveFolder, 'DataParams.mat'))
    summary.displayText = 'not aligned';
    summary.tooltipLines = [{'Alignment: no DataParams.mat (no reference selected).'}, ...
        iProjectNote(summary)];
    return
end
try
    DataParams = loadDataParams(saveFolder);
catch ME
    summary.state = 'Invalid';
    summary.issues = {ME.message};
    summary.displayText = 'invalid DataParams';
    summary.tooltipLines = {'Alignment: DataParams.mat could not be read:', ['  ' ME.message]};
    return
end

registration = DataParams.registration;
reference = registration.reference;
estimate = registration.estimate;
apply = registration.apply;
summary.isAligned = logical(apply.isRegistered);
if ~isempty(reference)
    summary.referenceName = char(string(reference.name));
    summary.sourceKind = reference.source.kind;
    summary.sourcePath = reference.source.path;
    summary.alignedSourceWarning = reference.source.warning;
    summary.shortHash = UMITLedgerCore.shortHash(reference.contentHash);
end

if summary.isAligned
    summary.state = 'Aligned';
    summary.appliedOn = apply.appliedOn;
    summary.appliedBy = apply.appliedBy;
    summary.isReviewed = logical(apply.isReviewed);
    if ~isempty(apply.appliedReferenceHash)
        summary.shortHash = UMITLedgerCore.shortHash(apply.appliedReferenceHash);
    end
    summary.displayText = sprintf('aligned to "%s" (%s)', summary.referenceName, summary.shortHash);
    lines = {sprintf('Alignment: ALIGNED (committed) to %s.', ...
        describeRegistrationReference(reference, apply.appliedReferenceHash))};
    lines{end+1} = sprintf('Applied on %s by %s.', apply.appliedOn, apply.appliedBy);
    if ~summary.isReviewed
        lines{end+1} = 'The alignment was not marked as reviewed.';
    end
    lines{end+1} = 'A different reference requires a fresh import of the raw data.';
elseif isempty(reference)
    summary.state = 'NoReference';
    summary.displayText = 'not aligned';
    lines = {'Alignment: not aligned; no reference selected.'};
elseif isempty(estimate)
    summary.state = 'ReferenceSelected';
    summary.displayText = sprintf('not aligned; reference "%s"', summary.referenceName);
    lines = {sprintf('Alignment: not aligned; %s selected, no transform estimated.', ...
        describeRegistrationReference(reference))};
elseif strcmp(estimate.referenceHash, reference.contentHash)
    summary.state = 'Estimated';
    summary.displayText = sprintf('not aligned; estimated against "%s"', summary.referenceName);
    lines = {sprintf('Alignment: not aligned; transform estimated against %s, not applied.', ...
        describeRegistrationReference(reference))};
else
    summary.state = 'EstimateStale';
    summary.displayText = sprintf('not aligned; reference "%s" changed since the estimate', ...
        summary.referenceName);
    lines = {sprintf(['Alignment: not aligned; %s selected, but the transform was ' ...
        'estimated against %s. Estimate again before applying.'], ...
        describeRegistrationReference(reference), ...
        describeRegistrationReference([], estimate.referenceHash))};
end

if ~isempty(reference)
    if ~isempty(summary.sourcePath)
        lines{end+1} = sprintf('Source (%s, as recorded): %s', summary.sourceKind, summary.sourcePath);
    end
    if ~isempty(summary.alignedSourceWarning)
        lines{end+1} = ['Warning recorded at selection: ' summary.alignedSourceWarning];
    end
end
summary.tooltipLines = [lines, iProjectNote(summary)];
end

function lines = iProjectNote(summary)
% Project-less mode: point out what projects add (brief Phase 5).
lines = {};
if ~summary.isProjectBound
    lines = {['This folder is not bound to a project. Projects let one ' ...
        'reference be reused across recordings (subject reference library).']};
end
end
