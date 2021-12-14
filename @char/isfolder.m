function tf = isfolder(Path)
% Function to overload Matlab's built-in function "isfile". For backward
% compatobility with versions prior to R2019a.
tf = false;
if exist(Path, 'dir')
    tf = true;
end
end

