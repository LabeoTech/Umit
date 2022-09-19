function Umitoolbox_setup
% This setup function adds the toolbox documentation to Matlab and updates
% the toolbox path for the current session.

% Get toolbox folder:
saveDir = fileparts(mfilename('fullpath'));
% Add toolbox to Path for the current Matlab session:
addpath(saveDir);
umIT_path_cleanup(saveDir); % Removes any folders from path containing the toolbox.
% Edit info.xml file to be able to access the documentation through Matlab "doc":
info_file = fileread(fullfile(saveDir,'docs','info.xml'));
path_str = regexp(info_file,'(?<=<help_location>).*(?=</help_location>)', 'match');
% Replace default path string with user's path:
info_file = strrep(info_file,path_str{:}, fullfile(saveDir,'docs'));
% Save info_file:
fid = fopen(fullfile(saveDir,'docs','info.xml'),'w');
fprintf(fid,'%s',info_file);
fclose(fid);
% Build serch database:
folders = dir(fullfile(saveDir, 'docs'));
help_dir = folders([folders.isdir] & startsWith({folders.name}, 'helpsearch-v'));
if ~isempty(help_dir)
    warning('off')
    status = rmdir(fullfile(help_dir.folder, help_dir.name), 's'); %#ok. Erase current database folder.
    warning('on')
end
try
    builddocsearchdb(fullfile(saveDir,'docs'))
catch 
    disp(repmat('-', 1,100))
    warning('Failed to create local documentation! Try to run MATLAB as admin and rerun this function!')    
    disp('If you do not have admin access, the same documentation is available in the project''s wiki page at:')
    disp('https://s-belanger.github.io/Umit/')
    disp(repmat('-', 1,100))
end
disp('Everything is set! You can start using the toolbox now!');
end

