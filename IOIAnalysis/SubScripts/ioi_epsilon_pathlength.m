function eps_pathlength = ioi_epsilon_pathlength(whichCurve,...
    baseline_hbt, baseline_hbo, baseline_hbr, Filters)
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

load('SysSpect.mat')

% Rough baseline concentrations (in uM) : 100 uM (in the brain)
c_tot = baseline_hbt*1e-6; %100e-6;

switch(Filters.Excitation)
    case 'GCaMP'
        c_led(1,:) = Red.*FF01496LP;
        c_led(2,:) = Green.*FF01496LP;
        c_led(3,:) = Yellow.*FF01496LP;
    case 'none'
        c_led(1,:) = Red;
        c_led(2,:) = Green;
        c_led(3,:) = Yellow;
    otherwise
        c_led(1,:) = Red;
        c_led(2,:) = Green;
        c_led(3,:) = Yellow;
end

switch(Filters.Emission)
    case 'GCaMP'
        c_led(1,:) = c_led(1,:).*FF01496LP;
        c_led(2,:) = c_led(2,:).*FF01496LP;
        c_led(3,:) = c_led(3,:).*FF01496LP;
    case 'jRGECO'
        c_led(1,:) = c_led(1,:).*FF01512630;
        c_led(2,:) = c_led(2,:).*FF01512630;
        c_led(3,:) = c_led(3,:).*FF01512630;
    otherwise
        c_led(1,:) = c_led(1,:);
        c_led(2,:) = c_led(2,:);
        c_led(3,:) = c_led(3,:);
end

switch(Filters.Camera)
    case 'CS2100M'
        c_camera = ThorQLux;
    case 'D1024'
        c_camera = PF1024;
    case 'D1312'
        c_camera = PF1312;
    otherwise
        c_camera = ones(1,301);
end
        
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


