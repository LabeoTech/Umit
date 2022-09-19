function IOIFiguresGen(ExpFolderPath)

load(fullfile(ExpFolderPath, 'Hb_Concentrations.mat'));
load(fullfile(ExpFolderPath, 'ROIs.mat'));
load([ExpFolderPath filesep 'StimParameters.mat'],'PreStimLength');

mask = ROIs{1}.mask;
C1 = permute(C1,[4 1 2 3]);

minhbo = min(C1(1,:,mask(:)));
minhbo = min(minhbo(:));
maxhbo = max(C1(1,:,mask(:)));
maxhbo = max(maxhbo(:));

minhbr = min(C1(2,:,mask(:)));
minhbr = min(minhbr(:));
maxhbr = max(C1(2,:,mask(:)));
maxhbr = max(maxhbr(:));

myfilter = fspecial('gaussian',[5 5], 2);
vidObj = VideoWriter(fullfile(ExpFolderPath,'response_2d.avi'));
vidObj.FrameRate = 15;
open(vidObj);

h = figure('Visible','off');
for j=1:size(C1,2)
    subplot(1,2,1)
    pat_overlay_blend(Map,imfilter(squeeze(C1(1,j,:,:)),myfilter,'replicate'),mask,[minhbo*0.4,maxhbo*0.4],[minhbo,minhbo*0.2,maxhbo*0.05,maxhbo]);
    axis equal;
    axis off;
    title(['time: ',num2str(j/5-5),' sec']);
    subplot(1,2,2)
    pat_overlay_blend(Map,imfilter(squeeze(C1(2,j,:,:)),myfilter,'replicate'),mask,[minhbr*0.4,maxhbr*0.4],[minhbr,minhbr*0.05,maxhbr*0.2,maxhbr]);
    axis equal;
    axis off;
    title(['time: ',num2str(j/5-PreStimLength),' sec']);
    writeVideo(vidObj, getframe(gcf));
end
close(vidObj);
close(h);

h=figure('Visible','off');
t=(1:250)/5-PreStimLength;
plot(t,squeeze(mean(C1(1,:,mask(:)),3)),'r')
hold on
plot(t,squeeze(mean(C1(2,:,mask(:)),3)),'b')
xlabel('Time (sec)')
saveas(h, fullfile(ExpFolderPath, 'avg_response.jpg'))
close(h);
end