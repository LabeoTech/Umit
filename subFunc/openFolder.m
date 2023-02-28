function openFolder(folder)
% This function opens the folder "folder" using the system's file explorer
% interface.

assert(isfolder(folder),'The folder doesn''t exist!')
switch computer
    case 'PCWIN64'
        winopen(folder)
    case 'MACI64'
        system(['open ' folder])
    case 'GLNXA64'
        system(['gnome-open ' folder])
    otherwise
        cd(folder)        
end

end
