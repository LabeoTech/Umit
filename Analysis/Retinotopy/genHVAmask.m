function outFile = genHVAmask(SaveFolder)
% GENHVAMASK uses the Azimuth and Elevation maps to segregate the visual
% areas.

% Load Azimuth and elevation maps:
az = load(fullfile(SaveFolder, 'AzimuthMap.dat'));
el = load(fullfile(SaveFolder, 'ElevationMap.dat'));
% Get phase component from maps:
az_phi = angle(az);
el_phi = angle(el);
% Create visual sign map:
VSM = genVSM(az_phi,el_phi);
% Filter VSM:
VSM = imgaussfilt(VSM,8);
% Threshold the VSM at 1.5 STD:
raw_bin = zeros(size(VSM), 'single');
thr = std(VSM(:));
raw_bin(VSM<-thr) = -1;
raw_bin(VSM>thr) = 1;
% Delineate visual cortex border:
SE = strel('disk',5);
vcb = imdilate(bwmorph(bwmorph(abs(raw_bin),'close',Inf),'open',Inf), SE, 'full');




end
% Local functions

function vsm = genVSM(phaseAz, phaseEl)
% This method calculates the visual sign maps of the Azimuth
% and Elevation maps.
disp('Calculating visual sign map...');

% Calculate phase map gradients
[gradAzx, gradAzy] = gradient(phaseAz);
[gradElx, gradEly] = gradient(phaseEl);
%
gradDirAz = atan2(gradAzy, gradAzx);
gradDirEl = atan2(gradEly, gradElx);
%
vsm = sin(angle(exp(1i.*gradDirAz).*exp(-1i.*gradDirEl)));
vsm(isnan(vsm)) = 0;
% Filter map
end