function HBLookupTableGen(Wavelength)

Optics = ReadConfigFile();
[e_hbo, e_hbr] = ioi_get_extinctions(400, 700, 301);
e_hbo = single(e_hbo);
e_hbr = single(e_hbr);
pathlength = ioi_path_length_factor(400, 700, 301, 0.100, 'Silico');
pathlength = single(pathlength);
SP = load('SysSpect.mat');

switch Optics.Camera
    case 'PF1024'
        Camera = SP.PF1024;
    case 'PF1312'
        Camera = SP.PF1312;
    case 'BFly' 
        Camera = SP.BFly;
    case 'CS2100M'
        Camera = SP.ThorQLux;
    otherwise
        Camera = SP.PF1024;
end

switch Wavelength
    case '475'
        eval(['LED = SP.' Optics.GFPIllum ';']);
        eval(['iFilter = SP.' Optics.GFPFilter_ex ';']);
        eval(['dFilter = SP.' Optics.GFPFilter_em ';']);
        
        Illumination = single(LED.*iFilter);
        Detection =  single(Camera.*dFilter);
        Fluo_ex = SP.GFP_ex;
        Fluo_em = SP.GFP_em;
        File = ['FactCorrHB_' Wavelength '.mat'];
    case '567'
        eval(['LED = SP.' Optics.RFPIllum ';']);
        eval(['iFilter = SP.' Optics.RFPFilter_ex ';']);
        eval(['dFilter = SP.' Optics.RFPFilter_em ';']);
        
        Illumination = single(LED.*iFilter);
        Detection =  single(Camera.*dFilter);
        Fluo_ex = SP.RFP_ex;
        Fluo_em = SP.RFP_em;
        File = ['FactCorrHB_' Wavelength '.mat'];
    otherwise
        eval(['LED = SP.' Optics.GFPIllum ';']);
        eval(['iFilter = SP.' Optics.GFPFilter_ex ';']);
        eval(['dFilter = SP.' Optics.GFPFilter_em ';']);
        
        Illumination = single(LED.*iFilter);
        Detection =  single(Camera.*dFilter);
        Fluo_ex = SP.GFP_ex;
        Fluo_em = SP.GFP_em;
        File = ['FactCorrHB_' Wavelength '.mat'];
end

[HbO, HbR] = meshgrid(single(0:0.25:120), single(0:0.25:120));

Iillz = Illumination.*Fluo_ex;
Iillz = sum(Iillz(:));
Idetz = Detection.*Fluo_em;
Idetz = sum(Idetz(:));
IllumFact = Illumination.*Fluo_ex./Iillz;
DetectFact = Detection.*Fluo_em./Idetz;

extinct = bsxfun(@times, HbO, permute(e_hbo,[3 1 2])) + ...
    bsxfun(@times, HbR, permute(e_hbr,[3 1 2]));
extinct = extinct*1e-6;
extinct = bsxfun(@times, extinct, permute(pathlength,[3 1 2]));
extinct = exp(-extinct);

tmp = bsxfun(@times, extinct, permute(IllumFact,[3 1 2]));
CorrIllum = sum(tmp,3);

tmp = bsxfun(@times, extinct, permute(DetectFact,[3 1 2]));
CorrDetect = sum(tmp,3);

Corr = 1./(CorrIllum.*CorrDetect);

Folder = mfilename('fullpath');
Folder = Folder(1:strfind(Folder,'\HBLookup'));
save([Folder File], 'Corr');
end

