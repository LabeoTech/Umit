function severity = setStatusMessage(statusLabel, message, varargin)
%SETSTATUSMESSAGE Display a color-coded status message in a MATLAB UI label.
%
%   severity = setStatusMessage(statusLabel, message)
%   severity = setStatusMessage(statusLabel, message, 'Severity', severity)
%   severity = setStatusMessage(statusLabel, message, 'Background', background)
%   severity = setStatusMessage(statusLabel, ME)
%
%   Updates a MATLAB App Designer status label using a severity-dependent
%   font color. The severity can be provided explicitly or inferred from the
%   message text.
%
%   Inputs:
%       statusLabel - MATLAB UI label handle, usually app.StatusLabel.
%       message     - Status text, string/char/cellstr, or MException.
%
%   Name-Value options:
%       Severity        - Message severity.
%                         Options:
%                             'auto'      Infer from message text.
%                             'neutral'   Normal status text.
%                             'info'      Informational message.
%                             'busy'      Running/loading/saving message.
%                             'success'   Successful operation.
%                             'warning'   Warning message.
%                             'error'     Error/failure message.
%                         Default: 'auto'
%
%       DefaultSeverity - Severity used when Severity='auto' and no keyword
%                         match is found.
%                         Default: 'info'
%
%       Background      - App background polarity used to choose readable
%                         status colors.
%                         Options:
%                             'auto'  Use the initialized figure theme,
%                                     otherwise light.
%                             'light'
%                             'dark'
%                         Default: 'auto'
%
%       MaxLength       - Maximum displayed text length. If 0, the full text
%                         is shown.
%                         Default: 0
%
%       UseTooltip      - If true, store the full message in the label
%                         tooltip when supported.
%                         Default: true
%
%       Refresh         - If true, call drawnow limitrate nocallbacks after
%                         updating the label.
%                         Default: true
%
%   Output:
%       severity        - Final severity used to color the label.
%
%   Examples:
%       setStatusMessage(app.StatusLabel, "Ready.");
%
%       setStatusMessage(app.StatusLabel, ...
%           "Warning: events.mat was not found.");
%
%       setStatusMessage(app.StatusLabel, ...
%           "Data exported successfully.", ...
%           "Severity", "success", ...
%           "Background", "light");
%
%       setStatusMessage(app.StatusLabel, ...
%           "Processing movie...", ...
%           "Severity", "busy", ...
%           "Background", "dark");
%
%       try
%           load(filePath, '-mat');
%       catch ME
%           setStatusMessage(app.StatusLabel, ME, ...
%               "Severity", "error", ...
%               "Background", "dark");
%       end
%
%   Notes:
%       - Explicit severity is preferred when the caller already knows the
%         message type.
%       - Text-based severity inference is intended as a fallback for legacy
%         code paths.
%       - The returned severity is useful for unit tests and debugging.
%       - Colors are RGB triplets in the [0, 1] range, compatible with
%         uilabel FontColor.

p = inputParser;
p.FunctionName = 'setStatusMessage';

addRequired(p, 'statusLabel', @(x) isobject(x) || isgraphics(x));
addRequired(p, 'message');

addParameter(p, 'Severity', 'auto', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'DefaultSeverity', 'info', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'Background', 'auto', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'MaxLength', 0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);

addParameter(p, 'UseTooltip', true, ...
    @(x) islogical(x) && isscalar(x));

addParameter(p, 'Refresh', true, ...
    @(x) islogical(x) && isscalar(x));

parse(p, statusLabel, message, varargin{:});

severityIn = lower(strtrim(char(string(p.Results.Severity))));
defaultSeverity = lower(strtrim(char(string(p.Results.DefaultSeverity))));
background = lower(strtrim(char(string(p.Results.Background))));

validSeverities = {'auto', 'neutral', 'info', 'busy', ...
    'success', 'warning', 'error'};

if ~ismember(severityIn, validSeverities)
    error('Umitoolbox:setStatusMessage:invalidSeverity', ...
        'Invalid severity "%s".', severityIn);
end

validDefaultSeverities = setdiff(validSeverities, {'auto'});
if ~ismember(defaultSeverity, validDefaultSeverities)
    error('Umitoolbox:setStatusMessage:invalidDefaultSeverity', ...
        'Invalid default severity "%s".', defaultSeverity);
end

validBackgrounds = {'auto', 'light', 'dark'};
if ~ismember(background, validBackgrounds)
    error('Umitoolbox:setStatusMessage:invalidBackground', ...
        'Invalid background "%s". Expected "auto", "light", or "dark".', ...
        background);
end

if strcmp(background, 'auto')
    background = iResolveFigureBackground(statusLabel);
end

[statusText, tooltipText, forcedSeverity] = iNormalizeStatusMessage(message);

if strcmp(severityIn, 'auto')
    if ~isempty(forcedSeverity)
        severity = forcedSeverity;
    else
        severity = iInferSeverity(statusText, defaultSeverity);
    end
else
    severity = severityIn;
end

displayText = statusText;

if p.Results.MaxLength > 0 && strlength(displayText) > p.Results.MaxLength
    displayText = extractBefore(displayText, p.Results.MaxLength);
    displayText = displayText + "...";
end

statusLabel.Text = char(displayText);
statusLabel.FontColor = iGetSeverityColor(severity, background);

if isprop(statusLabel, 'Tooltip')
    if p.Results.UseTooltip
        statusLabel.Tooltip = char(tooltipText);
    else
        statusLabel.Tooltip = '';
    end
end

if isprop(statusLabel, 'Interpreter')
    statusLabel.Interpreter = 'none';
end

if p.Results.Refresh
    drawnow limitrate nocallbacks
end

end

% =========================================================================
% Local helpers
% =========================================================================

function [statusText, tooltipText, forcedSeverity] = iNormalizeStatusMessage(message)
%INORMALIZESTATUSMESSAGE Convert supported message inputs to display text.

forcedSeverity = '';

if isa(message, 'MException')
    statusText = string(message.message);
    tooltipText = string(getReport(message, 'extended', 'hyperlinks', 'off'));
    forcedSeverity = 'error';
    return
end

if isstring(message)
    statusText = strjoin(message(:), newline);

elseif ischar(message)
    statusText = string(message);

elseif iscellstr(message)
    statusText = strjoin(string(message(:)), newline);

elseif iscell(message) && all(cellfun(@iIsTextScalar, message(:)))
    statusText = strjoin(string(message(:)), newline);

else
    try
        statusText = string(message);

        if ~isscalar(statusText)
            statusText = strjoin(statusText(:), " ");
        end

    catch
        statusText = string(evalc('disp(message)'));
    end
end

statusText = strtrim(statusText);

if strlength(statusText) == 0
    statusText = "Ready.";
end

tooltipText = statusText;

end

function tf = iIsTextScalar(x)
%IISTEXTSCALAR Return true for char vectors or scalar strings.

tf = ischar(x) || (isstring(x) && isscalar(x));

end

function severity = iInferSeverity(statusText, defaultSeverity)
%IINFERSEVERITY Infer status severity from message text.

txt = lower(char(statusText));

% Avoid false positives such as "No Errors" or "completed without warning".
if ~isempty(regexp(txt, '\b(no|without)\s+(error|errors|warning|warnings)\b', 'once'))
    severity = 'success';
    return
end

if ~isempty(regexp(txt, ...
        '\b(error|failed|failure|exception|invalid|aborted|unable|cannot|could not|unrecognized)\b', ...
        'once'))
    severity = 'error';
    return
end

if ~isempty(regexp(txt, ...
        '\b(warning|warn|skipped|deprecated|missing|not found)\b', ...
        'once'))
    severity = 'warning';
    return
end

if ~isempty(regexp(txt, ...
        '\b(done|completed|success|successful|saved|loaded|created|exported|imported)\b', ...
        'once'))
    severity = 'success';
    return
end

if ~isempty(regexp(txt, ...
        '\b(running|loading|saving|processing|validating|importing|exporting|please wait)\b', ...
        'once'))
    severity = 'busy';
    return
end

severity = defaultSeverity;

end

function color = iGetSeverityColor(severity, background)
%IGETSEVERITYCOLOR Return RGB color for a status severity and background.

switch background
    case 'dark'
        color = iGetDarkBackgroundColor(severity);

    otherwise
        color = iGetLightBackgroundColor(severity);
end

end

function background = iResolveFigureBackground(statusLabel)
%IRESOLVEFIGUREBACKGROUND Read the theme applied to the owning figure.

background = 'light';
figureHandle = ancestor(statusLabel, 'figure');

if isempty(figureHandle) || ~isappdata(figureHandle, 'UmitTheme')
    return
end

savedTheme = getappdata(figureHandle, 'UmitTheme');
if ischar(savedTheme) || (isstring(savedTheme) && isscalar(savedTheme))
    savedTheme = lower(strtrim(char(savedTheme)));
    if ismember(savedTheme, {'light', 'dark'})
        background = savedTheme;
    end
end

end

function color = iGetLightBackgroundColor(severity)
%IGETLIGHTBACKGROUNDCOLOR Return readable colors for light backgrounds.

switch severity
    case 'error'
        color = [0.80 0.00 0.00];   % red

    case 'warning'
        color = [0.70 0.45 0.00];   % dark yellow / amber

    case 'success'
        color = [0.00 0.45 0.15];   % dark green

    case 'busy'
        color = [0.00 0.25 0.70];   % dark blue

    case 'neutral'
        color = [0.25 0.25 0.25];   % gray

    case 'info'
        color = [0.10 0.10 0.10];   % near-black

    otherwise
        color = [0.10 0.10 0.10];   % fallback
end

end

function color = iGetDarkBackgroundColor(severity)
%IGETDARKBACKGROUNDCOLOR Return readable colors for dark backgrounds.

switch severity
    case 'error'
        color = [1.00 0.35 0.35];   % light red

    case 'warning'
        color = [1.00 0.75 0.20];   % light amber

    case 'success'
        color = [0.35 0.90 0.45];   % light green

    case 'busy'
        color = [0.45 0.70 1.00];   % light blue

    case 'neutral'
        color = [0.78 0.78 0.78];   % light gray

    case 'info'
        color = [0.92 0.92 0.92];   % near-white

    otherwise
        color = [0.92 0.92 0.92];   % fallback
end

end
