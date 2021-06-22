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



