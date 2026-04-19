function eps_pathlength = ioi_epsilon_pathlength(whichCurve,...
    baseline_hbt, baseline_hbo, baseline_hbr, FilterSetName,CameraModel)
%	This function estimates epsilon * D, it takes into account the camera
%	response, the leds spectra and uses a pathlength factor either set from Kohl
%	or Dunn in the literature.
%   This module is dependent on this file which contains all hardware info for
%   the setup, needs to specify the leds and the camera response. we are still a
%   bit dependent on the RGY but we could abstract this (however the lambdas
%   would need to be registered to specific hardware still
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

dictSpectra = load('SysSpect.mat');
allFilterSets = load('FilterSets.mat');
cameraInfo = load('CameraSpect.mat');
% Get FilterSet spectra:
Filters.Excitation = 1;
Filters.Emission = 1;
if ~strcmpi(FilterSetName,'none')
    fn = fieldnames(allFilterSets);
    thisFilter = allFilterSets.(fn{strcmpi(fn,FilterSetName)});
    if ~strcmpi(thisFilter.Excitation,'none')
        Filters.Excitation = dictSpectra.(thisFilter.Excitation);
    end
    if ~strcmpi(thisFilter.Emission, 'none')
        Filters.Emission = dictSpectra.(thisFilter.Emission);
    end
end
allCameras = fieldnames(cameraInfo);
idxCam = strcmpi(CameraModel,allCameras);
if any(idxCam)
    c_camera = dictSpectra.(cameraInfo.(allCameras{idxCam}));
else
    c_camera = ones(1,301);
end

% Rough baseline concentrations (in uM) : 100 uM (in the brain)
c_tot = baseline_hbt*1e-6; %100e-6;
% Set Excitation:
c_led(1,:) = dictSpectra.Red.*Filters.Excitation;
c_led(2,:) = dictSpectra.Green.*Filters.Excitation;
c_led(3,:) = dictSpectra.Yellow.*Filters.Excitation;
% Set Emission
c_led(1,:) = c_led(1,:).*Filters.Emission;
c_led(2,:) = c_led(2,:).*Filters.Emission;
c_led(3,:) = c_led(3,:).*Filters.Emission;
%
c_pathlength = ioi_path_length_factor(400, 700, 301, c_tot*1000, whichCurve);
[c_ext_hbo,c_ext_hbr] = ioi_get_extinctions(400, 700, 301);

% Create vectors of values for the fits
CHbO = baseline_hbo/baseline_hbt*c_tot*(.85:.1:1.15);
CHbR = baseline_hbr/baseline_hbt*c_tot*(.85:.1:1.15);

% In this computation below we neglect the fact that pathlength changes
% with total concentration (it is fixed for a Ctot of 100e-6)
eps_pathlength=zeros(3,2);

IHbO = zeros(size(CHbO));
IHbR = zeros(size(CHbR));

%gMuS = abs(gMuS);
%Mus = gMuS(1)*(lambda_vec/500).^(-gMuS(2));
for iled=1:3
    for iconc = 1:length(CHbO)
        IHbO(iconc) = sum(c_camera.*c_led(iled,:).*exp(-c_pathlength*CHbO(iconc).*...
            (c_ext_hbo)),2); %	Measured intensity for different concentrations
        IHbR(iconc) = sum(c_camera.*c_led(iled,:).*exp(-c_pathlength*CHbR(iconc).*...
            (c_ext_hbr)),2);
    end
    IHbO = IHbO/median(IHbO);
    IHbR = IHbR/median(IHbR);
    
    % Compute effective eps
    p1 = polyfit(CHbO,-log(IHbO),1);
    p2 = polyfit(CHbR,-log(IHbR),1);
    HbRL = p2(1); %epsilon*D HbR effectif
    HbOL = p1(1);%epsilon*D HbO effectif
    eps_pathlength(iled,1)=HbOL;
    eps_pathlength(iled,2)=HbRL;
end


