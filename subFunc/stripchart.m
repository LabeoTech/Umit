function stripchart(x,y,varargin)
% STRIPCHART strip plot for categorical variables.
% This type of plot is built-in in Maltab R2021b and above. Thus, this is
% an option for those with older Matlab versions.
% This function creates a jittered scatter plot around discrete (integer) X
% axis points to reduce too much overlap between the data points.
% Inputs:
%   x (int): column vector of x axis positions. This corresponds to the
%   index of the categorical variable.
%   y (float): column vector of data.
% Optional inputs:
%   group (int| default = 1): Column vector of extra group indices. The x axis will be split
%   between groups.
%   color (char or 3xN array | default = 'b'): color of scatter plot. If an array is
%   provided, it must have the same size of "y".
%   Xjitter (bool | default = false): Apply jitter in X axis. If TRUE, data points
%   will be jittered randomly with an uniform distribution. Otherwise, the
%   data points will be plotted over a single X axis point.
%   axHandle (axis handle | default = GCA ): handle to the axis where the
%   data will be plotted.

%%% Arguments parsing and validation %%%
p = inputParser;
% The input of the function must be a File , RawFolder or SaveFolder
validateNumVec = @(x)isnumeric(x) & isvector(x);
addRequired(p,'x',validateNumVec); % Validate if the input is a 3-D numerical matrix:
addRequired(p,'y',validateNumVec); % Validate if the input is a numerical vector:
addOptional(p,'varType','std',@(x) ismember(lower(x), {'sem', 'std','ci'}));
addParameter(p,'group', 1, validateNumVec)
addParameter(p,'color', 'b',@(x) ( ischar(x) || ismatrix(x))  || ( iscell(x) && ischar(x{:}) ));
addParameter(p,'Xjitter',false,@islogical)
addParameter(p,'Boxplot',false,@islogical);
addParameter(p,'axHandle',[],@(x) isa(x, 'matlab.graphics.axis.Axes')| isempty(x));
% Parse inputs:
parse(p,x, y, varargin{:});
%Initialize Variables:
x = p.Results.x;
y = p.Results.y;
varType = p.Results.varType;
group = p.Results.group;
color = p.Results.color;
b_Xjitter = p.Results.Xjitter;
b_showBox = p.Results.Boxplot;
axHandle = p.Results.axHandle;
clear p
%%%%
% Further input validation:
if isempty(axHandle)
    % Create new figure:
    figure; axHandle = gca;
end
[x,y,color,group] = validateVector(x,y,color,group);
% Extra Parameters:
lenSpread = 0.8; % Data spreads 80% of the space between X points.
MarkerAlpha = 0.5; % Dot transparency.
MarkerEdgeColor = 'k'; % Dot edge color.
MarkerSize = 20; % Dot size.
ErrorLineWidth = 3; % Width of error bars.
% Get number of groups per X point:
nG = numel(unique(group));
% Then we calculate the center of each data column around an X point.
xVec = unique(x);

xSpread = [xVec - lenSpread/2, xVec + lenSpread/2];
xSpread = unique(xSpread);xSpread = reshape(xSpread,2,[])';
xLim = zeros(size(xSpread,1),nG+1);
for ii = 1:size(xSpread,1)
    xLim(ii,:) = linspace(xSpread(ii,1), xSpread(ii,2), nG+1);
end
xCtr = (xLim(:,2:end) - xLim(:,1:end-1))./2 + xLim(:,1:end-1);
spacing = 0.2*min(diff(xLim(1,:))); % 20% of spacing between groups.

% Plotting:
hold(axHandle, 'on');
for ii = 1:xVec(end)
    idxX = ( x == xVec(ii) );
    for jj = 1:nG
        idxG = ( group == jj );
        data = y(idxX & idxG);
        if isempty(data)
            continue
        end
        myColor = unique(color(idxX & idxG,:),'rows');
        % Create Errorbar in X-centers:
        % Group lines in a single graphical object:
        if b_showBox
            errG = plotGroup(axHandle, data, xCtr(ii,jj), (lenSpread/nG) - spacing, 'box');
            % Customize the colors and line widths of boxplot elements:
            medLine = findobj(errG,'Tag','Median'); set(medLine, 'LineWidth',ErrorLineWidth,'Color',myColor);
            avgLine = findobj(errG, 'Tag','Mean'); set(avgLine, 'LineWidth',ErrorLineWidth, 'Color', myColor);
            boxPatch = findobj(errG,'Tag','Box'); set(boxPatch, 'FaceColor',myColor);                        
        else
            errG = plotGroup(axHandle, data,xCtr(ii,jj), (lenSpread/nG) - spacing, 'err', varType);
            arrayfun(@(x) set(x, 'Color',myColor,'LineWidth', 0.7*ErrorLineWidth), errG.Children);
            idxAvgLine = strcmp({errG.Children.Tag},'avgLine');
            set(errG.Children(idxAvgLine), 'LineWidth', ErrorLineWidth);
        end
        % Append the info structure with the original X axis position of
        % the data:
        errG.Tag = ['errG' num2str(jj) 'X' num2str(ii)];
        errG.UserData.Xpos = x(idxX & idxG); % Store data's initial X position in object's UserData.
        % Create scatter plot in selected positions:
        xPos = rand(size(data,1),1);
        if b_Xjitter
            xPos = rescale(xPos,xLim(ii,jj)+spacing/2,xLim(ii,jj+1)-spacing/2);
        else
            xPos = repmat(xCtr(ii,jj),length(data),1);
        end
        s = scatter(axHandle, xPos, data(:));
        s.MarkerFaceColor = myColor;
        s.MarkerFaceAlpha = MarkerAlpha;
        s.MarkerEdgeColor = MarkerEdgeColor;
        s.SizeData = MarkerSize;
        s.UserData.Xpos = x(idxX & idxG); % Store data's initial X position in UserData.
        s.Tag = ['scatG' num2str(jj) 'X' num2str(ii)];
    end
end
hold(axHandle, 'off');
end
%--------------------------------------------------------------------------
function [x,y,c,g] = validateVector(x,y,c,g)
% Check if all inputs have the same dimensions.

if isscalar(g)
    g = repmat(g, size(x));
end

if ischar(c)
    c = {c};
end
if iscell(c)
    if numel(unique(c)) ~= numel(unique(g))
        ME = MException('stripchart:invalidInput', '"Color" input must correspond to the number of groups.');
        throwAsCaller(ME)
    end
    col = zeros(length(x),3);
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
if ( ~isequaln(length(x), length(y), size(c,1), length(g)) )
    ME = MException('stripchart:invalidInput', 'Input vectors must have the same size.');
    throwAsCaller(ME)
end
% Flip all inputs as columns:
if size(x,1) < size(x,2)
    x = x';
end

if size(y,1) < size(y,2)
    y = y';
end

if size(c,1) < size(c,2)
    c = c';
end

if size(g,1) < size(g,2)
    g = g';
end
end

function hg = plotGroup(ax,data,xCtr, xRange,plotType, varargin)
errType = [varargin{:}];
xRange = xRange/2;
hg = hggroup('Parent', ax);
info = struct();
if strcmpi(plotType, 'err')
    avgDat = mean(data,'omitnan');
    line(ax,[xCtr - xRange xCtr + xRange],[avgDat avgDat],'Parent',hg,'Tag', 'avgLine');
    switch errType
        case 'std'
            errDat = std(data,0,'all','omitnan');
        case 'sem'
            errDat = std(data,0,'all','omitnan')./sqrt(length(data));
        case 'ci'
            errDat = (tinv(1-0.025,length(data)-1)).*(std(data,0,'all','omitnan')./sqrt(length(data)));
    end
    line(ax,[xCtr xCtr],[avgDat - errDat, avgDat + errDat],'Parent',hg);
    line(ax,[xCtr - 0.6*xRange, xCtr + 0.6*xRange],...
        [avgDat - errDat, avgDat - errDat; avgDat + errDat, avgDat + errDat]',...
        'Parent',hg,'Tag','capLine');
    info.Mean = avgDat;
    info.(errType) = errDat;
else
    data = sort(data);
    % Calculate boxplot elements:
    med = median(data,'omitnan');
    avg = mean(data,'omitnan');
    q1 = quantile(data,.25); q3 = quantile(data,.75);
    IQR = q3 - q1; % Interquartile range.
    wUb = 1.5*IQR + q3; % Upper Whisker position;
    wLb = q1 - 1.5*IQR; % Lowe Whisker position;
    ub = data(find(data <= wUb,1,'last'));% Whisker upper boundary.
    lb = data(find(data >= wLb,1,'first')); % Whisker lower boundary.
    % Draw boxplot:
    bx = patch(ax,[xCtr - xRange, xCtr + xRange, xCtr + xRange, xCtr - xRange], ...
        [q1, q1, q3, q3],'k', 'Parent',hg, 'Tag','Box');
    bx.EdgeColor = 'k';
    bx.FaceAlpha = .3;
    % Draw median line:
    line(ax,[xCtr - xRange, xCtr + xRange], [med med], 'Parent',hg, 'Tag', 'Median');
    % Draw mean line:
    line(ax,[xCtr - xRange, xCtr + xRange], [avg avg], 'Parent',hg, 'LineStyle',':','Tag', 'Mean');
    % Draw whiskers:
    line(ax,[xCtr, xCtr],[q1,lb],'Parent', hg, 'Tag','Lower Whisker');
    line(ax,[xCtr, xCtr],[q3,ub],'Parent', hg, 'Tag','Upper Whisker');
    % Draw whiskers caps:
    line(ax,[xCtr - xRange/2, xCtr + xRange/2],[lb lb],'Parent', hg, 'Tag','Lower Adjacent Value');
    line(ax,[xCtr - xRange/2, xCtr + xRange/2],[ub ub],'Parent', hg, 'Tag','Upper Adjacent Value');
    % Gather all boxplot values in a structure and store it in the Group's
    % UserData;
    info = struct();
    info.IQR = IQR;
    info.UpperAdjacent = ub;
    info.LowerAdjacent = lb;
    info.Mean = avg;
    info.Median = med;
    info.q1 = q1;
    info.q3 = q3;
end
% Save plot's data to Group UserData:
hg.UserData = info;
end

