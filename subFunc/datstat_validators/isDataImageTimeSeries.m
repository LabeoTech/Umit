function tf = isDataImageTimeSeries(dataStruct)
% Validator function for .DATSTAT file.
%
% ISDATAIMAGETIMESERIES checks if all fields in dataCategory are equal to 'Image-time-series'.

tf = all(cellfun(@(x) strcmpi(dataStruct.dataCategory.(x),'Image-time-series'),fieldnames(dataStruct.dataCategory)));
if ~tf
    error('The data should be an image time series!')
end

end