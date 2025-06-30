function interpolatedValues = importSpectralProfileFromTXT(filename)
% Open the file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open the file: %s', filename);
end

% Initialize a cell array to store numeric data
numericData = [];

% Read the file line by line
while true
    line = fgetl(fid);
    if ~ischar(line)  % End of file
        break;
    end
    
    % Split the line by comma or tab
    splitLine = regexp(line, '[,\t]', 'split');
    
    % Convert to numeric values
    numericValues = str2double(splitLine);
    
    % Check if the line contains valid numeric data
    if all(~isnan(numericValues))
        numericData = [numericData; numericValues]; % Append as a new row
    end
end

% Close the file
fclose(fid);

% Check if numericData is empty
if isempty(numericData)
    error('No numeric data found in the file.');
end

% Extract wavelengths and corresponding values
wavelengths = numericData(:, 1);
values = numericData(:, 2);

% Normalize the second column (values) to range from 0 to 1
values = (values - min(values)) / (max(values) - min(values));

% Extract values within the wavelength range of 400 to 700
validIndices = wavelengths >= 400 & wavelengths <= 700;
if isempty(validIndices)
    error('Wavelenghts are out of bound. The file should contain spectral data from 400nm to 700nm')
end
filteredWavelengths = wavelengths(validIndices);
filteredValues = values(validIndices);

% Interpolate the data to have 301 values
newWavelengths = linspace(400, 700, 301);
interpolatedValues = interp1(filteredWavelengths, filteredValues, newWavelengths, 'linear', 'extrap');

end
