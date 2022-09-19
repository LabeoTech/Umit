%This function converts to raw binary files acquired by the imaging system
%to tiff (or big tiff files).
%No image classification to de-interlace the channels (colors) is done.
function ConvertBinToTiff(FolderPath)
%The folder path is passed as an argument. If not, a dialog will ask the
%user to select a path
if (nargin < 1)
    FolderPath = uigetdir(pwd,'Select a folder containing a set of image data.');
end

if ~ischar(FolderPath)
    error('Invalid folder.')
end
    
disp(['Starting Tiff conversion for: ' FolderPath]);

%Check how many data files there are
nFiles = 0;
totalSize = 0; %Used to know if we need to use big tif over 4 GB
while (exist([FolderPath filesep 'img_' sprintf('%05d',nFiles) '.bin'],'file'))
    file = dir([FolderPath filesep 'img_' sprintf('%05d',nFiles) '.bin']);
    totalSize = totalSize + file.bytes;
    nFiles = nFiles + 1;
end
disp([num2str(nFiles) ' files found.'])

outFName = strcat(FolderPath, filesep, 'img.tif'); %Output file name
if (totalSize < 3900000000)
    fTIF = Fast_Tiff_Write(outFName,1,0);
else
    fTIF = Fast_BigTiff_Write(outFName,1,0); %BigTiff over 4 GB (3.9 GB to be safe and account for headers)
    disp('Using BigTiff format')
end

for idxFile = 1:nFiles
    fid = fopen([FolderPath filesep 'img_' sprintf('%05d',idxFile-1) '.bin'],'r'); %Open the file
    header = fread(fid, 5, 'int32'); %Each file contains an header of 5 int32
    
    for idxImage = 0:header(5)-1 %Read all images in the file, header(5) is the number of images in the current fime
        %Seek the correct position in the file. Each image as an header of 3 uint64
        %header(2) is the width, header(3) is the height
        fseek(fid, 5*4+idxImage*(3*8+header(2)*header(3)*2)+3*8, -1);

        %Read the image
        theImage = uint16(fread(fid, header(2)*header(3), 'uint16'));

        if (~isempty(theImage))
            im = reshape(theImage,header(2),header(3));

            %Legacy Matlab tif write (slow)
%             imwrite(im, outFName, 'WriteMode', 'append',  'Compression','none');

            %Tif write using https://github.com/rharkes/Fast_Tiff_Write
            %Requires "Fast_Tiff_Write.m" and "Fast_BigTiff_Write.m"
            fTIF.WriteIMG(im);

            %OPTIONAL: Accumulate the image in an array (takes a lot of memory)
%             if (~exist('imArray'))
%                 imArray = [];
%                 imArray(:,:,1) = im;
%             else
%                 imArray(:,:,end+1) = im;
%             end

            %OPTIONAL: Display each image
%             imagesc(im);
%             title(num2str(idxImage + (idxFile-1)*header(5))) %The title is the image index
%             pause(0.01)            

            %OPTIONAL: Progress message
            if (mod(idxImage+1, 100) == 0)
                %The number of images per file is not equal to header(5)
                %for the last file, so we recalculate it
                file = dir([FolderPath filesep 'img_' sprintf('%05d',idxFile-1) '.bin']);
                nImagesPerFile = (file.bytes-5*4)/(3*8+header(2)*header(3)*2);
                disp([num2str(idxImage+1) '/' num2str(nImagesPerFile) ' images done on file ' num2str(idxFile) ' of ' num2str(nFiles)])
            end
        end
    end
    fclose(fid);
end
fTIF.close;
disp('Tiff conversion is done.')
