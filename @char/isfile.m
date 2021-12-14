function tf = isfile(fileName)
% Function to overload Matlab's built-in function "isfolder". For backward
% compatobility with versions prior to R2019a.
tf = false;
if exist(fileName, 'file')
    tf = true;
end
end

