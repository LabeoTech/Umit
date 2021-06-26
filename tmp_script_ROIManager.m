load('D:\Academico\PostDoc_UdeM\LabeoTech\IsaToolbox\GUI\DataViz\mouse_ctx_borders.mat')
% Splitting ROIs:
newCA = false([size(atlas.CA,1) size(atlas.CA,2) 2*size(atlas.CA,3)]);
for i = 1:size(atlas.CA,3)
    [r1,c1] = find(atlas.CA(:,:,i),1,'first');
    [r2,c2] = find(atlas.CA(:,:,i),1,'last');
%     middle = round(mean([c1,c2]));
    middle = round(mean([r1,r2]));
    newCA(:,middle+1:end,i) = atlas.CA(:,middle+1:end,i);
    newCA(:,1:middle-1,i+size(atlas.CA,3)) = atlas.CA(:,1:middle-1,i);
end

figure;
for i = 1:size(atlas.CA,3)/2
    imagesc(any(cat(3,atlas.CA(:,:,i), atlas.CA(:,:,(i+(size(atlas.CA,3)/2)))),3));
    axis image; colormap gray;
%     imagesc(any(cat(3, newCA(:,:,i), newCA(:,:,i+(size(atlas.CA,3)))),3)); 
%     axis image; colormap gray
    title(i);
    pause(.5)
end

figure;
perim_mask = boundarymask(newCA(:,:,2));

x_all = [];
y_all = [];
for i = 1:size(atlas.CA,3)
    B = bwboundaries(atlas.CA(:,:,i), 8, 'noholes');
    X = cellfun(@(x) x(:,1), B, 'UniformOutput', false);
    Y = cellfun(@(x) x(:,2), B, 'UniformOutput', false);
    x_all = [x_all; X];
    y_all = [y_all; Y];
end
flat_polyShape = polyshape(y_all, x_all);

figure;
hold on
for i = 1:length(polyArr)
    p = plot(polyArr);
end
hold off
axis image


%% Find bounding rectangle;
idxMap = atlas.flatAtlas ~=0;
bb = regionprops(idxMap, 'BoundingBox');

%% Test everything:

load('E:\Solenn_Data\ProjectSaveDir\M4X\ImagingReferenceFrame.mat');
%
close all
%
fig = figure;
ax = axes('Parent', fig);

img = imagesc(ax,reference_frame); colormap jet; axis image;

ax.NextPlot = 'replacechildren';
img.HandleVisibility = 'off';
% pHandle = plot(atlas.polyArr);
pHandle = plot(ax,atlas.polyShapes.polyArr);
% pHandle.FaceAlpha = .2;
% pHandle.FaceColor = [.3 .3 .3];
rect = drawrectangle(ax,'Position', atlas.BoundingBox, 'FixedAspectRatio', 1, ...
    'Rotatable', 1, 'DrawingArea','unlimited', 'Deletable', false);
addlistener(rect,'MovingROI',@(src,evt) updatePolyShapePos(src,evt,pHandle));
addlistener(rect,'ROIClicked',@(src,evt) closeRect(src,evt));

waitfor(rect);

polyShape = [pHandle.Shape];
% delete(pHandle);
delete(rect);


%% Crap code
% Create circle around point:
xCenter = atlas.BX;
yCenter = atlas.BY;
theta = 0 : 0.1 : 2*pi;
radius = 6;
x = radius * cos(theta) + xCenter;
y = radius * sin(theta) + yCenter;





% Highlight Pixels inside ROI (polyshape)
% Run ROI positioning tool before!
ptest = polyShape(end);

imag_test = reference_frame;
[r,c] = find(imag_test);
idx = ptest.isinterior(c,r);
vals = arrayfun(@(x,y) imag_test(y,x), c(idx), r(idx)) 

   
figure; imagesc(imag_test); axis image
hold on 
plot(c(idx), r(idx), 'r.');
hold off
%%
X = c(idx);
Y = r(idx);
for i = 1:length(X)
    x = X(i);
    y = Y(i);
    imag_test(y,x) = -10;
end

figure;
imagesc(imag_test) ; axis image



imag(i)

figure;
cx = 250.8562;
cy = 711.3073;
ax = gca;
ax.XLim = [100 400];
ax.YLim = [550 850];
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
pOrig = atlas.polyShapes.polyArr(end);
theta = 0 : 0.1 : 2*pi;
for i = 1:5:120
    new_x = i * cos(theta) + cx;
    new_y = i * sin(theta) + cy;
    shape = polyshape(new_x, new_y);
    pOut = intersect(shape,pOrig);
    p = plot(pOut);
    title(i);
    pause(.5)

end
    
%% Playing with Java and UITable
root = getenv('Umitoolbox')
fig = uifigure;
Data = table([1:10]',[15:24]', repmat('<body bgcolor = "#FF0000">    </body>',10,1))
uit = uitable('Parent', fig, 'Data', Data);
uit = uitable('Data',{'<body bgcolor="#FF0000">Hello</body>'})
cmap = jet(10);




% Initialize our custom cell renderer class object
javaaddpath(fullfile(root, 'GUI', 'DataViz', 'ColoredFieldCellRenderer.zip'));
cr = ColoredFieldCellRenderer(java.awt.Color.white);
cr.setDisabled(true);  % to bg-color the entire column
% Set specific cell colors (background and/or foreground)
for i = 1:10
    cr.setCellBgColor(2,i-1, java.awt.Color(cmap(i,1), cmap(i,2), cmap(i,3)));
end
% Replace Matlab's table model with something more renderer-friendly...
jTable.setModel(javax.swing.table.DefaultTableModel());
set(uit,'ColumnFormat',[]);
% Finally assign the renderer object to all the table columns
for colIdx = 1 : numel(uit.ColumnName)
  jTable.getColumnModel.getColumn(colIdx-1).setCellRenderer(cr);
end

%% Testing NaN region exclusion from ROI
close all
% Test image:
imag = reference_frame;
% add region with NaNs:
imag(50:70,40:60) = NaN;
imag(60,50) = 1;
imag(61,50) = 1;

% Plot:
figure;
ax =gca;
imH = imagesc(ax,imag); axis image;
% Create Polyshape by hand:
rect = drawpolygon(ax);
pS_orig = polyshape(rect.Position);
delete(rect);
hold on
p0 = plot(pS_orig);
p0.FaceColor = [0 1 0]; p0.EdgeColor = [0 0 0];

% Create non_nan_msk:
not_nan_msk = ~isnan(imag);

% Create mask from polyshape using POLY2MASK:
pS_orig_msk = poly2mask(pS_orig.Vertices(:,1), pS_orig.Vertices(:,2), ...
    size(imag,1), size(imag,2));
% Create intersect_msk:
intersect_msk = not_nan_msk & pS_orig_msk;
% Visualize:
hold(ax,'on')
imagesc(intersect_msk, 'AlphaData', .5*intersect_msk);
% Create new polyshape from intersect_msk:
% bwboundaries method:
B = bwboundaries(intersect_msk);
X = cellfun(@(x) x(:,1), B, 'UniformOutput', false);
Y = cellfun(@(x) x(:,2), B, 'UniformOutput', false);
pS_intersect = polyshape(Y,X);

%bwtraceboundary method:
% [x_init, y_init] = find(intersect_msk,1, 'first');
% tmp = bwtraceboundary(intersect_msk,[x_init y_init],'E');
% pS_intersect = polyshape(tmp(:,2), tmp(:,1));
% Visualize:
p1 = plot(pS_intersect);
hold(ax,'off');

% Chech if POLY2MASK of polyshape derived from maks degrades the ROI too
% much:

% POLY2MASK method:

% nReg = pS_intersect.NumRegions;
% test_msk = false([size(imag), nReg]);
% tmp_ps = regions(pS_intersect);
% for i = 1:nReg
%     test_msk(:,:,i) = poly2mask(tmp_ps(i).Vertices(:,1), tmp_ps(i).Vertices(:,2), ...
%     size(imag,1), size(imag,2));
% end
% test_msk = any(test_msk,3);

% Create isinterior method:
[r,c] = find(imag);
idx = pS_intersect.isinterior(c,r);
test_msk = false(size(imag));
X = c(idx);
Y = r(idx);
for i = 1:size(X,1)
    test_msk(Y(i),X(i)) = true;
end



figure;
imshowpair(intersect_msk, test_msk)



%%$$ OK!! But there is some erosion!

%% Controlling for ROIS partially out of image:

close all
figure;
ax =gca;
ax.Color = [0 0 0];
imH = imagesc(ax,imag); axis image;
% Create Polyshape by hand:
rect = drawpolygon(ax);
pS_orig = polyshape(rect.Position);
delete(rect);
hold on
p0 = plot(pS_orig);
p0.FaceColor = 'none'; p0.EdgeColor = [1 1 1];
% Create polyshape with image size:
[h,w] = size(imag);
pS_imag = polyshape((.5 + [0 0 w w]), (.5 + [h 0 0 h]));
p1 = plot(pS_imag);
p1.EdgeColor = [0 1 0];
p1.FaceColor = "none";
pS_intersect = pS_orig.intersect(pS_imag);

p2 = plot(pS_intersect);
p2.FaceColor = [1 0 0];
p2.EdgeColor = [0 0 1];

%%$$ OK!!  
%%








% Create NaN mask:
nan_msk = isnan(imag);
% Create a polyshape from nan_msk:
[x_init, y_init] = find(nan_msk,1, 'first');
tmp = bwtraceboundary(nan_msk,[x_init y_init],'E');
pS_nan= polyshape(tmp(:,2), tmp(:,1));
pS_subtract= pS_orig.subtract(pS_nan);

p1 = plot([pS_nan, pS_subtract]);

p1(1).EdgeColor = [0 0 0];
p1(2).EdgeColor = [0 1 0];
% Find pixels inside pS_subtract:
[r,c] = find(imag);
idx_in = pS_subtract.isinterior(c,r);
% Create not_nan_msk:
not_nan_msk = false(size(imag));
X = c(idx_in);
Y = r(idx_in);
for i = 1:size(X,1)
    not_nan_msk(Y(i),X(i)) = true;
end
% overlay not_nan_msk:
imagesc(not_nan_msk, 'AlphaData', .3.*not_nan_msk)


% Create Mask from polyshape:
[r,c] = find(imag);
idx = pS_orig.isinterior(c,r);
ps_msk = false(size(imag));
X = c(idx);
Y = r(idx);
for i = 1:size(X,1)
    ps_msk(Y(i),X(i)) = true;
end
% Plot ROI pixels:
imagesc(ps_msk, 'AlphaData', .3*ps_msk)
hold off
% Create mask from intersection:
intersect_msk = not_nan_msk.*ps_msk;
% Create polyshape from mask:
[x_init, y_init] = find(intersect_msk,1, 'first');
tmp = bwtraceboundary(intersect_msk,[x_init y_init],'E');
pS_rebuild = polyshape(tmp(:,2), tmp(:,1));
% Validate all:
close all
figure;
imagesc(imag); axis image;
hold on
imagesc(intersect_msk, 'AlphaData', .2.*intersect_msk);
p1 = plot(pS_orig); p1.FaceColor = 'none'
plot(pS_rebuild)
for i = 1:length(B)
    vtcs = B{i};
    p_all(i) = plot(polyshape(vtcs(:,2), vtcs(:,1)))
end


figure;
imagesc(intersect_msk); axis image
hold on
p0 = plot(tmp(:,2), tmp(:,1));
p1 = plot(pS_orig);
p1.FaceColor = 'none';
p1.EdgeColor = [1 1 1];
p2 = plot(pS_rebuild);

% Compare pixels inside:
[r,c] = find(imag);
b_isinOrig = pS_orig.isinterior(c,r);
b_isinRebuild = pS_rebuild.isinterior(c,r);

% Visualize differences:
msk_orig = false(size(intersect_msk));
msk_rebuild = msk_orig;
X_o = r(b_isinOrig);
Y_o = c(b_isinOrig);
X_r = r(b_isinRebuild);
Y_r = c(b_isinRebuild);
for i = 1:length(X_o)
    msk_orig(X_o(i),Y_o(i)) = true;
end
    
for i = 1:length(X_r)
    msk_rebuild(X_r(i),Y_r(i)) = true;
end
figure;
imagesc(msk_rebuild .*msk_orig)
imshowpair(msk_orig, msk_rebuild)
% Comparing ISINTERIOR vs poly2mask:
ps_rebuil_msk = poly2mask(pS_rebuild.Vertices(:,1), pS_rebuild.Vertices(:,2), ...
    size(imag,1), size(imag,2));

figure;
imshowpair(msk_rebuild, ps_rebuil_msk, 'falsecolor', 'Scaling', 'joint')

figure;
imshowpair(msk_orig, ps_msk, 'falsecolor', 'Scaling', 'joint')

%% Testing ROI shrinkage

close all
figure;
ax =gca;
ax.Color = [0 0 0];
imH = imagesc(ax,imag); axis image;
% Create Polyshape by hand:
rect = drawpolygon(ax);
pS_orig = polyshape(rect.Position);
delete(rect);
hold on
p0 = plot(pS_orig);
hold off
[r,c] = find(imag);
idx = pS_orig.isinterior(c,r);
ps_msk = false(size(imag));
X = c(idx);
Y = r(idx);
for i = 1:size(X,1)
    ps_msk(Y(i),X(i)) = true;
end


eroded_msk = bwmorph(ps_msk, 'shrink', 1);
imshowpair(ps_msk,eroded_msk, 'falsecolor', 'Scaling', 'joint');

%% Test DRAWPOINT

close all
figure;
ax =gca;
ax.Color = [0 0 0];
imH = imagesc(ax,imag); axis image;
% Create point:
pt= drawpoint(ax);

% New origin:
newOrigin = round(pt.Position);

% Shift X-Axis:
xLab = 0:20:ax.XLim(2);
xLab = unique([-xLab xLab]);
xLoc = xLab + newOrigin(1);
idx = xLoc >= 0;

set(ax, 'XTick', xLoc(idx), 'XTickLabel', arrayfun(@(x) num2str(x), xLab(idx),...
    'UniformOutput', false));
% Shift Y-Axis:
yLab = 0:20:ax.YLim(2);
yLab = unique([-yLab yLab]);
yLoc = yLab + newOrigin(2);
idx = yLoc >= 0;

set(ax, 'YTick', yLoc(idx), 'YTickLabel', arrayfun(@(x) num2str(x), yLab(idx),...
    'UniformOutput', false));
line([newOrigin(1) newOrigin(1)], ax.YLim, 'Color','k')
line(ax.XLim, [newOrigin(2) newOrigin(2)], 'Color','k')


xLab = unique([-xLab xLab]);
xPos = xLab + newOrigin(1);
idx = xPos>=0;


xLab = sort([-xLab xLab(2:end)]);
xLab = xLab + newOrigin(1);
xLab = xLab(xLab>=0);


newX = ax.XTick - newOrigin(1);
newY = ax.YTick - newOrigin(2);
newXLabel = arrayfun(@(x) num2str(x), newX, 'UniformOutput', false);
newYLabel = arrayfun(@(x) num2str(x), newY, 'UniformOutput', false);
ax.XTickLabel = newXLabel;
ax.YTickLabel = newYLabel;


















