function eps_pathlength = ioi_epsilon_pathlength(whichCurve, ...
    baseline_hbt, baseline_hbo, baseline_hbr, opticalInfo)
%IOI_EPSILON_PATHLENGTH Estimate effective extinction times pathlength.
%
%   epsPathlength = hemooptics.ioi_epsilon_pathlength(curve, HbT, HbO, HbR, opticalInfo)
%
%   Consumes spectra already resolved by UMITRigStore
%   (resolveOpticalConfiguration). The legacy form that loaded the IOI
%   optical MAT files by filter-set and camera name was removed when
%   HemoCompute absorbed IOIAnalysis (2026-10-14); those files are now only
%   the source of RigManagement/migrateLegacyOpticalLibrary.
%
%   opticalInfo is a scalar structure with these required fields:
%       wavelengthNm          301x1 wavelengths, exactly (400:700)'.
%       activeRows            3x1 logical mask. Rows are always ordered
%                             red, green, yellow; at least two are active.
%       illuminationResponse  3x301 illumination spectra in fixed row order.
%       cameraResponse        3x301 camera spectra selected per acquisition
%                             channel/CamIdx in the same fixed row order.
%       excitationResponse    301-sample excitation-filter transmission.
%       emissionResponse      301-sample emission-filter transmission.
%
%   Active spectral rows must be finite. Inactive illumination/camera rows
%   may contain NaN because they are not used. Responses use the canonical
%   400:700-nm grid and are normally normalized to [0,1] by the resolver.
%
%   epsPathlength is 3x2. Rows retain the fixed red/green/yellow order and
%   inactive rows are NaN. Column 1 is the effective HbO coefficient and
%   column 2 is the effective HbR coefficient.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

requiredFields = {'wavelengthNm', 'activeRows', 'illuminationResponse', ...
    'cameraResponse', 'excitationResponse', 'emissionResponse'};
assert(isstruct(opticalInfo) && isscalar(opticalInfo) && ...
    all(isfield(opticalInfo, requiredFields)), ...
    'Umitoolbox:ioi_epsilon_pathlength:InvalidOpticalInfo', ...
    'opticalInfo is missing required resolved-spectrum fields.');
assert(isequal(opticalInfo.wavelengthNm(:), (400:700)') && ...
    isequal(size(opticalInfo.illuminationResponse), [3 301]) && ...
    isequal(size(opticalInfo.cameraResponse), [3 301]) && ...
    numel(opticalInfo.excitationResponse) == 301 && ...
    numel(opticalInfo.emissionResponse) == 301 && ...
    numel(opticalInfo.activeRows) == 3, ...
    'Umitoolbox:ioi_epsilon_pathlength:InvalidOpticalInfo', ...
    'Resolved spectra must use three rows and the canonical 301-point grid.');

activeRows = logical(opticalInfo.activeRows(:));
assert(sum(activeRows) >= 2, ...
    'Umitoolbox:ioi_epsilon_pathlength:InvalidOpticalInfo', ...
    'At least two resolved optical rows are required.');
filters = opticalInfo.excitationResponse(:)' .* opticalInfo.emissionResponse(:)';
c_led = opticalInfo.illuminationResponse .* filters;
c_camera = opticalInfo.cameraResponse;

% Rough baseline concentrations (in uM) : 100 uM (in the brain)
c_tot = baseline_hbt*1e-6; %100e-6;
c_pathlength = ioi_path_length_factor(400, 700, 301, c_tot*1000, whichCurve);
[c_ext_hbo,c_ext_hbr] = ioi_get_extinctions(400, 700, 301);

% Create vectors of values for the fits
CHbO = baseline_hbo/baseline_hbt*c_tot*(.85:.1:1.15);
CHbR = baseline_hbr/baseline_hbt*c_tot*(.85:.1:1.15);

% In this computation below we neglect the fact that pathlength changes
% with total concentration (it is fixed for a Ctot of 100e-6)
eps_pathlength = nan(3,2);

IHbO = zeros(size(CHbO));
IHbR = zeros(size(CHbR));

%gMuS = abs(gMuS);
%Mus = gMuS(1)*(lambda_vec/500).^(-gMuS(2));
for iled = find(activeRows)'
    assert(all(isfinite(c_led(iled, :))) && all(isfinite(c_camera(iled, :))), ...
        'Umitoolbox:ioi_epsilon_pathlength:InvalidOpticalInfo', ...
        'Active optical rows must contain finite illumination and camera spectra.');
    for iconc = 1:length(CHbO)
        IHbO(iconc) = sum(c_camera(iled,:).*c_led(iled,:).*exp(-c_pathlength*CHbO(iconc).*...
            (c_ext_hbo)),2); %	Measured intensity for different concentrations
        IHbR(iconc) = sum(c_camera(iled,:).*c_led(iled,:).*exp(-c_pathlength*CHbR(iconc).*...
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
end
