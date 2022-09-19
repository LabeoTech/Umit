Folder = uigetdir();
if( ~exist([Folder filesep 'Speckle_flow.mat'],'file') )
    disp('Folder does not contain speckle data...');
    return;
end
load([Folder filesep 'Speckle_flow.mat']);
load([Folder filesep 'StimParameters.mat']);
figure(2); colormap jet;
T = linspace(-PreStimLength, InterStim_min, size(speckle,1));  
for indF = 1:size(speckle,1)
     imagesc(squeeze(speckle(indF,:,:,1)),[0.8 1.2]); %scale is fixed between 80% and 120% if mean flow
     title(['Time(s): ' num2str(T(indF))]);
     pause(0.2);
end