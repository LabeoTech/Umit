function varargout = plotTimeVector(dat,varargin)
% PLOTTIMEVECTOR strip plot for categorical variables.
% This type of plot is built-in in Maltab R2021b and above. Thus, this is
% an option for those with older Matlab versions.
% This function creates a jittered scatter plot around discrete (integer) X
% axis points to reduce too much overlap between the data points.
% Inputs:
%   data (int): matrix containing the data per observation (rows) and frames (columns).
%   Optional inputs:
%   varType (char | 'std'): name of variation measure for error bars. It
%       should be one of the following: "std" (standard deviation), "sem"
%       (standard error of the mean) or "ci" (95% confidence interval).
%   plotType ({'shaded','line'} | default = 'shaded'): Type of plot.
%   "Shaded" adds a shaded area around the average group line with the
%   data variation (defined by "varType"). "Line" shows individual lines
%   of "data".
%   group (int| default = 1): Column vector of group indices. The plot lines will be
%       separated by groups.
%   color (char or 3xN array | default = 'b'): color of each group.
%   axHandle (axis handle | default = [] ): handle to the axis where the
%       data will be plotted. If no handle is provided, a new figure will
%       be created with the plot.
%   dataTips (cell | default = {''}): cell array of characters containing
%       text associated with the data point in the scatter plot. This info
%       will be appended to the element's UserData structure and can be
%       used to customize a data tip or to be able to easily retrieve meta
%       data associated with a data point in the plot.
%   Freq (positive number| default = []): Sampling frequency. If provided,
%       the X axis will be transformed in seconds.
% Outputs:
%   Optional output:
%   axHandle : handle to axis containing the plotted data.

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
validateNumMat = @(x)isnumeric(x) & ismatrix(x);
addRequired(p,'data',validateNumMat); % Validate if the input is a numerical matrix.
addOptional(p,'varType','std',@(x) ismember(lower(x), {'sem', 'std','ci'}));
addOptional(p,'plotType','shaded',@(x) ismember(lower(x), {'shaded', 'line'}));
addParameter(p,'group', 1, validateNumMat)
addParameter(p,'color', 'b',@(x) ( ischar(x) || ismatrix(x))  || ( iscell(x) && ischar(x{:}) ));
addParameter(p,'axHandle',[],@(x) isa(x, 'matlab.graphics.axis.Axes')| isempty(x));
addParameter(p,'dataTip',{''},@iscell);
addParameter(p,'Freq',[], @(x) isempty(x) || ( isscalar(x) && x>0 ));
% Parse inputs:
parse(p,dat, varargin{:});
%Initialize Variables:
data = p.Results.data;
varType = lower(p.Results.varType);
plotType = p.Results.plotType;
group = p.Results.group;
color = p.Results.color;
axHandle = p.Results.axHandle;
dataTip = p.Results.dataTip;
Freq = p.Results.Freq;
clear p
%%%%
% Further input validation:
if isempty(axHandle)
    % Create new figure:
    figure; axHandle = gca;
end
[color,group,dataTip] = validateMatrix(data, color, group, dataTip);
ShadeAlpha = 0.2; % Patch transparency.
% Get number of groups per X point:
nG = numel(unique(group));
% Plotting:
hold(axHandle, 'on');

for ii = 1:nG
    hg = hggroup('Parent', axHandle);
    idxG = ( group == ii );
    dat = data(idxG,:);    
    if isempty(dat(~isnan(dat)))
        continue
    end
    myColor = unique(color(idxG,:),'rows');    
    % Create line plots:
    avg = mean(dat,1,'omitnan');
    switch varType
        case 'std'
            errDat = std(dat,0,1,'omitnan');
        case 'sem'
            errDat = std(dat,0,1,'omitnan')./sqrt(size(dat,1));
        case 'ci'
            errDat = (tinv(1-0.025,size(dat,1)-1)).*(std(dat,0,1,'omitnan')./sqrt(size(dat,1)));
    end
    % Plot Shaded area OR individual lines:    
     x = 1:size(dat,2);
     if ~isempty(Freq)
         x = x./Freq;
     end
    if strcmpi(plotType, 'shaded')
        % Plot shaded area:       
        p = patch(axHandle,[x(~isnan(avg)) fliplr(x(~isnan(avg)))],[avg(~isnan(avg)) - errDat(~isnan(errDat)),...
            fliplr(avg(~isnan(avg)) + errDat(~isnan(errDat)))], myColor,'Parent',hg);
        set(p,'FaceAlpha',ShadeAlpha,...
            'EdgeColor', myColor,...
            'Tag',['ShadeG' num2str(ii)]);        
    else
        % Plot individual lines:
        Tips = dataTip(idxG);
        for jj = 1:size(dat,1)
            p(jj) = plot(axHandle,x,dat(jj,:),'Color',myColor, 'Tag',['dataLn' num2str(jj), 'G', num2str(ii)], 'Parent',hg);
            p(jj).Color(4) = ShadeAlpha;
            p(jj).UserData.DataTip = Tips(jj);
        end
    end
    % Plot Average line on top of previous plot(s):
    avgLn(ii) = plot(x,avg, 'Color', myColor,'LineWidth',2,'Tag',['avgLnG', num2str(ii)], 'Parent', hg);%#ok. Otherwise, I cannot get a handle array.
    hg.Tag = ['timeVecG' num2str(ii)];  
    hg.UserData.Mean = avg;
    hg.UserData.std = std(dat,0,1,'omitnan');% Be sure to always have the standard deviation field.
    hg.UserData.Ndata = size(dat,1); % Store the number of data used to calculate the error bar. This can be used to recalculate the error later.
    if ~strcmpi(varType, 'std')
        hg.UserData.(varType) = errDat;
    end
end
hold(axHandle, 'off');
if nargout
    varargout = {axHandle};
end
end
%--------------------------------------------------------------------------
function [c,g,dt] = validateMatrix(x,c,g,dt)
% Check if all inputs have the same dimensions.

if isempty([dt{:}])
    dt = repmat({'null'}, size(x,1),1);
end
if isscalar(g)
    g = repmat(g, size(x,1),1);
end

if ( isnumeric(c) && numel(unique(g)) == size(c,1) )
    col = zeros(size(g,1),3);
    for ii = 1:size(col,1)
        col(ii,:) = c(g(ii),:);
    end
    c = col;
end
if strcmpi(c,'b') && numel(unique(g)) ~= 1
    % If no optional color was provided, create random colors for each
    % group:
    c = distinguishable_colors(numel(unique(g)), 'w');    
    col = zeros(size(g,1),3);
    for ii = 1:size(g,1)
        col(ii,:) = c(g(ii),:);
    end
    c = col; clear col        
end

if ischar(c) && numel(unique(g)) == 1
    c = {c};
end

if iscell(c)
    if numel(unique(c)) ~= numel(unique(g))
        ME = MException('stripchart:invalidInput',...
            '"Color" input must correspond to the number of groups.');
        throwAsCaller(ME)
    end
    col = zeros(size(x,1),3);
    for ii = 1:length(c)
        switch c{ii}
            case 'y'
                rgb = [1 1 0];
            case 'm'
                rgb = [1 0 1];
            case 'c'
                rgb = [0 1 1];
            case 'r'
                rgb = [1 0 0];
            case 'g'
                rgb = [0 1 0];
            case 'b'
                rgb = [0 0 1];
            case 'w'
                rgb = [1 1 1];
            case 'k'
                rgb = [0 0 0];
        end
        col(g == ii,:) = repmat(rgb, sum(g == ii),1);
    end
    c = col;
end
if ( ~isequaln(size(x,1), size(c,1), length(g), length(dt)) )
    ME = MException('stripchart:invalidInput', 'Input data must have the same size.');
    throwAsCaller(ME)
end
% Flip sime inputs as columns:

if size(g,1) < size(g,2)
    g = g';
end

if size(dt,1) < size(dt,2)
    dt = dt';
end
end